#include <bsd/stdlib.h> // arc4random_buf

#include "rdpf.hpp"
#include "bitutils.hpp"

#undef RDPF_MTGEN_TIMING_1

#ifdef RDPF_MTGEN_TIMING_1
// Timing tests for multithreaded generation of RDPFs
// nthreads = 0 to not launch threads at all
// run for num_iters iterations, output the number of millisections
// total for all of the iterations
//
// Results: roughly 50 µs to launch the thread pool with 1 thread, and
// roughly 30 additional µs for each additional thread.  Each iteration
// of the inner loop takes about 4 to 5 ns.  This works out to around
// level 19 where it starts being worth it to multithread, and you
// should use at most sqrt(2^{level}/6000) threads.
static void mtgen_timetest_1(nbits_t level, int nthreads,
    size_t num_iters, const DPFnode *curlevel,
    DPFnode *nextlevel, size_t &aes_ops)
{
    if (num_iters == 0) {
        num_iters = 1;
    }
    size_t prev_aes_ops = aes_ops;
    DPFnode L = _mm_setzero_si128();
    DPFnode R = _mm_setzero_si128();
    // The tweak causes us to compute something slightly different every
    // iteration of the loop, so that the compiler doesn't notice we're
    // doing the same thing num_iters times and optimize it away
    DPFnode tweak = _mm_setzero_si128();
    auto start = boost::chrono::steady_clock::now();
    for(size_t iter=0;iter<num_iters;++iter) {
        tweak += 1;  // This actually adds the 128-bit value whose high
                     // and low 64-bits words are both 1, but that's
                     // fine.
        size_t curlevel_size = size_t(1)<<level;
        if (nthreads == 0) {
            size_t laes_ops = 0;
            for(size_t i=0;i<curlevel_size;++i) {
                DPFnode lchild, rchild;
                prgboth(lchild, rchild, curlevel[i]^tweak, laes_ops);
                L = (L ^ lchild);
                R = (R ^ rchild);
                nextlevel[2*i] = lchild;
                nextlevel[2*i+1] = rchild;
            }
            aes_ops += laes_ops;
        } else {
            DPFnode tL[nthreads];
            DPFnode tR[nthreads];
            size_t taes_ops[nthreads];
            size_t threadstart = 0;
            size_t threadchunk = curlevel_size / nthreads;
            size_t threadextra = curlevel_size % nthreads;
            boost::asio::thread_pool pool(nthreads);
            for (int t=0;t<nthreads;++t) {
                size_t threadsize = threadchunk + (size_t(t) < threadextra);
                size_t threadend = threadstart + threadsize;
                boost::asio::post(pool,
                    [t, &tL, &tR, &taes_ops, threadstart, threadend,
                    &curlevel, &nextlevel, tweak] {
                        DPFnode L = _mm_setzero_si128();
                        DPFnode R = _mm_setzero_si128();
                        size_t aes_ops = 0;
                        for(size_t i=threadstart;i<threadend;++i) {
                            DPFnode lchild, rchild;
                            prgboth(lchild, rchild, curlevel[i]^tweak, aes_ops);
                            L = (L ^ lchild);
                            R = (R ^ rchild);
                            nextlevel[2*i] = lchild;
                            nextlevel[2*i+1] = rchild;
                        }
                        tL[t] = L;
                        tR[t] = R;
                        taes_ops[t] = aes_ops;
                    });
                threadstart = threadend;
            }
            pool.join();
            for (int t=0;t<nthreads;++t) {
                L ^= tL[t];
                R ^= tR[t];
                aes_ops += taes_ops[t];
            }
        }
    }
    auto elapsed =
        boost::chrono::steady_clock::now() - start;
    std::cout << "timetest_1 " << int(level) << " " << nthreads << " "
        << num_iters << " " << boost::chrono::duration_cast
        <boost::chrono::milliseconds>(elapsed) << " " <<
        (aes_ops-prev_aes_ops) << " AES\n";
    dump_node(L);
    dump_node(R);
}

#endif
