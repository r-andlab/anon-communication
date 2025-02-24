#ifndef __HEAPSAMPLER_HPP__
#define __HEAPSAMPLER_HPP__

#include <vector>
#include "mpcio.hpp"
#include "coroutine.hpp"
#include "heap.hpp"

// Implement a stream sampler using a MinHeap.  A stream sampler will
// sample a random subset of k elements from an arbitrality long stream
// of elements, while using only O(k) memory.  The number of elements in
// the stream need not be known in advance.  Importantly, no party knows
// which elements ended up in the sample.

// We use the technique from "Path Oblivious Heap" by Elaine Shi
// (IEEE Symsposium on Security and Privacy 2020):
// https://eprint.iacr.org/2019/274.pdf; in particular, Section 5.2
// "Oblivious Streaming Sampler with Applications to Distributed
// Differential Privacy".  We correct the typo that the probability to
// keep the mth item (for m>k) is listed as 1/m, but it should be k/m.
// Note that the Shi paper is in the client-server setting, where the
// client _is_ allowed to know which elements did and did not end up in
// the sample (but the server, who stores the heap, is not). In our MPC
// setting, no party may learn which elements ended up in the sample.

class HeapSampler {
    // The number of items to sample
    size_t k;

    // The number of items that have arrived so far
    size_t m;

    // The MinHeap with O(k) storage that when m>=k stores a
    // uniformly random sample of size k of the m items seen so far
    MinHeap heap;

    // The next random tag to use
    RegAS randtag;

    // Make the next random tag
    void make_randtag(MPCTIO &tio, yield_t &yield);

public:
    // Constructor for a HeapSampler that samples k items from a stream
    // of abritrary and unknown size, using O(k) memory
    HeapSampler(MPCTIO &tio, yield_t &yield, size_t k);

    // An element has arrived
    void ingest(MPCTIO &tio, yield_t &yield, RegAS elt);

    // The stream has ended; output min(k,m) randomly sampled elements.
    // After calling this function, the HeapSampler is reset to its
    // initial m=0 state.
    std::vector<RegAS> close(MPCTIO &tio, yield_t &yield);
};

// A unit test for the HeapSampler

void heapsampler_test(MPCIO &mpcio, const PRACOptions &opts, char **args);

// A unit test for the weighted_coin internal function

void weighted_coin_test(MPCIO &mpcio, const PRACOptions &opts, char **args);

#endif
