#include "heapsampler.hpp"

// In each 64-bit RegAS in the heap, the top bit (of the reconstructed
// value) is 0, the next 42 bits are the random tag, and the low 21 bits
// are the element value.

#define HEAPSAMPLE_TAG_BITS 42
#define HEAPSAMPLE_ELT_BITS 21

// Make the next random tag
void HeapSampler::make_randtag(MPCTIO &tio, yield_t &yield)
{
    // Make a uniformly random HEAPSAMPLE_TAG_BITS-bit tag. This needs
    // to be RegXS in order for the sum (XOR) of P0 and P1's independent
    // values to be uniform.
    RegXS tagx;
    tagx.randomize(HEAPSAMPLE_TAG_BITS);
    mpc_xs_to_as(tio, yield, randtag, tagx);
}

// Compute the heap size (the smallest power of two strictly greater
// than k) needed to store k elements
static size_t heapsize(size_t k)
{
    size_t ret = 1;
    while (ret <= k) {
        ret <<= 1;
    }
    return ret;
}

// Return a random bit that reconstructs to 1 with probability k/m
static RegBS weighted_coin(MPCTIO &tio, yield_t &yield, size_t k,
    size_t m)
{
    RegAS limit;
    limit.ashare = size_t((__uint128_t(k)<<63)/m) * !tio.player();
    RegXS randxs;
    randxs.randomize(63);
    RegAS randas;
    mpc_xs_to_as(tio, yield, randas, randxs);
    CDPF cdpf = tio.cdpf(yield);
    auto[lt, eq, gt] = cdpf.compare(tio, yield, randas-limit, tio.aes_ops());

    return lt;
}

// Constructor for a HeapSampler that samples k items from a stream
// of abritrary and unknown size, using O(k) memory
HeapSampler::HeapSampler(MPCTIO &tio, yield_t &yield, size_t k)
    : k(k), m(0), heap(tio.player(), heapsize(k))
{
    run_coroutines(tio, [&tio, this](yield_t &yield) {
        heap.init(tio, yield);
    }, [&tio, this](yield_t &yield) {
        make_randtag(tio, yield);
    });
}

// An element has arrived
void HeapSampler::ingest(MPCTIO &tio, yield_t &yield, RegAS elt)
{
    ++m;
    RegAS tagged_elt = (randtag << HEAPSAMPLE_ELT_BITS) + elt;
    RegAS elt_to_insert = tagged_elt;

    if (m > k) {
        RegAS extracted_elt;
        RegBS selection_bit;
        run_coroutines(tio, [&tio, this, &extracted_elt](yield_t &yield) {
            extracted_elt = heap.extract_min(tio, yield);
        }, [&tio, this, &selection_bit](yield_t &yield) {
            selection_bit = weighted_coin(tio, yield, k, m);
        });
        mpc_select(tio, yield, elt_to_insert, selection_bit,
            extracted_elt, tagged_elt);
    }

    run_coroutines(tio, [&tio, this, elt_to_insert](yield_t &yield) {
        heap.insert_optimized(tio, yield, elt_to_insert);
    }, [&tio, this](yield_t &yield) {
        make_randtag(tio, yield);
    });
}

// The stream has ended; output min(k,m) randomly sampled elements.
// After calling this function, the HeapSampler is reset to its
// initial m=0 state.
std::vector<RegAS> HeapSampler::close(MPCTIO &tio, yield_t &yield)
{
    size_t retsize = k;
    if (m < k) {
        retsize = m;
    }
    std::vector<RegAS> ret(retsize);
    for (size_t i=0; i<retsize; ++i) {
        ret[i] = heap.extract_min(tio, yield);
        ret[i] &= ((size_t(1)<<HEAPSAMPLE_ELT_BITS)-1);
    }
    // Compare each output to (size_t(1)<<HEAPSAMPLE_ELT_BITS), since
    // there may be a carry; if the output is greater than or equal to
    // that value, fix the carry
    RegAS limit;
    limit.ashare = (size_t(1)<<HEAPSAMPLE_ELT_BITS) * !tio.player();

    std::vector<coro_t> coroutines;
    for (size_t i=0; i<retsize; ++i) {
        coroutines.emplace_back([&tio, &ret, i, limit](yield_t &yield) {
            CDPF cdpf = tio.cdpf(yield);
            auto[lt, eq, gt] = cdpf.compare(tio, yield,
                ret[i]-limit, tio.aes_ops());
            RegAS fix, zero;
            mpc_select(tio, yield, fix, gt^eq, zero, limit);
            ret[i] -= fix;
        });
    }
    run_coroutines(tio, coroutines);
    heap.init(tio, yield);
    return ret;
}

void heapsampler_test(MPCIO &mpcio, const PRACOptions &opts, char **args)
{
    size_t n = 100;
    size_t k = 10;

    // The number of elements to stream
    if (*args) {
        n = atoi(*args);
        ++args;
    }
    // The size of the random sample
    if (*args) {
        k = atoi(*args);
        ++args;
    }

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);
    run_coroutines(tio, [&mpcio, &tio, n, k] (yield_t &yield) {

        std::cout << "\n===== STREAMING =====\n";

        HeapSampler sampler(tio, yield, k);

        for (size_t i=0; i<n; ++i) {
            // For ease of checking, just have the elements be in a
            // simple sequence
            RegAS elt;
            elt.ashare = (i+1) * (1 + 0xfff*tio.player());
            sampler.ingest(tio, yield, elt);
        }

        std::vector<RegAS> sample = sampler.close(tio, yield);

        tio.sync_lamport();
        mpcio.dump_stats(std::cout);
        mpcio.reset_stats();
        tio.reset_lamport();

        std::cout << "\n===== CHECKING =====\n";
        size_t expected_size = k;
        if (n < k) {
            expected_size = n;
        }
        assert(sample.size() == expected_size);
        std::vector<value_t> reconstructed_sample(expected_size);

        std::vector<coro_t> coroutines;
        for (size_t i=0; i<expected_size; ++i) {
            coroutines.emplace_back(
                [&tio, &sample, i, &reconstructed_sample](yield_t &yield) {
                    reconstructed_sample[i] = mpc_reconstruct(
                        tio, yield, sample[i]);
                });
        }
        run_coroutines(tio, coroutines);
        if (tio.player() == 0) {
            for (size_t i=0; i<expected_size; ++i) {
                printf("%06lx\n", reconstructed_sample[i]);
            }
        }
    });
}

void weighted_coin_test(MPCIO &mpcio, const PRACOptions &opts, char **args)
{
    size_t iters = 100;
    size_t m = 100;
    size_t k = 10;

    // The number of iterations
    if (*args) {
        iters = atoi(*args);
        ++args;
    }
    // The denominator
    if (*args) {
        m = atoi(*args);
        ++args;
    }
    // The numerator
    if (*args) {
        k = atoi(*args);
        ++args;
    }

    MPCTIO tio(mpcio, 0, opts.num_cpu_threads);
    run_coroutines(tio, [&mpcio, &tio, iters, m, k] (yield_t &yield) {

        size_t heads = 0, tails = 0;
        for (size_t i=0; i<iters; ++i) {
            RegBS coin = weighted_coin(tio, yield, k, m);
            bool coin_rec = mpc_reconstruct(tio, yield, coin);
            if (coin_rec) {
                heads++;
            } else {
                tails++;
            }
            printf("%lu flips %lu heads %lu tails\n", i+1, heads, tails);
        }
    });
}
