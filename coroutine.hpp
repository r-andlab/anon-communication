#ifndef __COROUTINE_HPP__
#define __COROUTINE_HPP__

#include <vector>

#include "corotypes.hpp"
#include "mpcio.hpp"


// The top-level coroutine runner will call run_coroutines with
// a MPCTIO, and we should call its send() method.  Subcoroutines that
// launch their own coroutines (and coroutine running) will call
// run_coroutines with a yield_t instead, which we should just call, in
// order to yield to the next higher level of coroutine runner.
static inline void send_or_yield(MPCTIO &tio) { tio.send(); }
static inline void send_or_yield(yield_t &yield) { yield(); }

// Get and set communication_nthreads for an MPCTIO; for a yield_t, this
// is a no-op.
static inline int getset_communication_nthreads(MPCTIO &tio, int nthreads = 0) {
    return tio.comm_nthreads(nthreads);
}
static inline int getset_communication_nthreads(yield_t &yield, int nthreads = 0) {
    return 0;
}

// Use this version if you have a variable number of coroutines (or a
// larger constant number than is supported below).
template <typename T>
inline void run_coroutines(T &mpctio_or_yield, std::vector<coro_t> &coroutines) {
    // If there's more than one coroutine, at most one of them can have
    // communication_nthreads larger than 1 (see mpcio.hpp for details).
    // For now, we set them _all_ to 1 (if there's more than one of
    // them), and restore communication_nthreads when they're all done.

    int saved_communication_nthreads = 0;
    if (coroutines.size() > 1) {
        saved_communication_nthreads =
            getset_communication_nthreads(mpctio_or_yield, 1);
    }

    // Loop until all the coroutines are finished
    bool finished = false;
    while(!finished) {
        // If this current function is not itself a coroutine (i.e.,
        // this is the top-level function that launches all the
        // coroutines), here's where to call send().  Otherwise, call
        // yield() here to let other coroutines at this level run.
        send_or_yield(mpctio_or_yield);
        finished = true;
        for (auto &c : coroutines) {
            // This tests if coroutine c still has work to do (is not
            // finished)
            if (c) {
                finished = false;
                // Resume coroutine c from the point it yield()ed
                c();
            }
        }
    }

    if (saved_communication_nthreads > 0) {
        getset_communication_nthreads(mpctio_or_yield,
            saved_communication_nthreads);
    }
}

// Use one of these versions if you have a small fixed number of
// coroutines.  You can of course also use the above, but the API for
// this version is simpler.

template <typename T>
inline void run_coroutines(T &mpctio_or_yield, const coro_lambda_t &l1)
{
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(l1);
    run_coroutines(mpctio_or_yield, coroutines);
}

template <typename T>
inline void run_coroutines(T &mpctio_or_yield, const coro_lambda_t &l1,
    const coro_lambda_t &l2)
{
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(l1);
    coroutines.emplace_back(l2);
    run_coroutines(mpctio_or_yield, coroutines);
}

template <typename T>
inline void run_coroutines(T &mpctio_or_yield, const coro_lambda_t &l1,
    const coro_lambda_t &l2, const coro_lambda_t &l3)
{
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(l1);
    coroutines.emplace_back(l2);
    coroutines.emplace_back(l3);
    run_coroutines(mpctio_or_yield, coroutines);
}

template <typename T>
inline void run_coroutines(T &mpctio_or_yield, const coro_lambda_t &l1,
    const coro_lambda_t &l2, const coro_lambda_t &l3,
    const coro_lambda_t &l4)
{
    std::vector<coro_t> coroutines;
    coroutines.emplace_back(l1);
    coroutines.emplace_back(l2);
    coroutines.emplace_back(l3);
    coroutines.emplace_back(l4);
    run_coroutines(mpctio_or_yield, coroutines);
}

#endif
