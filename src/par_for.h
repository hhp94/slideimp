#ifndef PAR_FOR_H
#define PAR_FOR_H

#include <RcppThread.h>

#include <cstddef>

// -----------------------------------------------------------------------------
// par_for - RcppThread::parallelFor with a serial path for a single thread.
// -----------------------------------------------------------------------------
// RcppThread's free parallelFor() reads the global pool's thread count, sets it
// to the requested count, runs the loop, and sets it back
// (RcppThread/parallelFor.hpp). quickpool implements that set as: join every
// worker, then spawn the new count (quickpool.hpp, set_active_threads). The
// global pool starts at hardware_concurrency(), so a one-thread call rebuilds
// the pool twice - down to one worker and back up - and that cost is fixed: it
// does not shrink as the work does.
//
// One thread means there is nothing to schedule, so the body runs inline and
// the pool is never touched. That is not a rare case here: slide_imp() runs
// each window at one core, and group_imp() forces one core per worker under
// mirai, so the whole chain is single-core per call by default.
//
// Interruption and errors are unchanged. The loop body calls
// RcppThread::checkUserInterrupt(), which polls R directly when it is on the
// main thread and throws UserInterruptException; in the parallel path that
// same exception is raised in a worker and rethrown by quickpool in the owner
// thread. Either way it leaves through the same place.
//
// Anything called from `f` must still obey the rule that the R API is off
// limits, since the parallel path is the one that decides what is legal.
template <typename F>
inline void par_for(
    const std::size_t begin,
    const std::size_t end,
    F f,
    const std::size_t n_threads,
    const std::size_t n_batches = 0)
{
    if (n_threads <= 1)
    {
        for (std::size_t i = begin; i < end; ++i)
        {
            f(i);
        }
        return;
    }

    RcppThread::parallelFor(
        static_cast<int>(begin),
        static_cast<int>(end),
        f,
        n_threads,
        n_batches);
}

#endif
