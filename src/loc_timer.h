#ifndef LOC_TIMER_H
#define LOC_TIMER_H
#ifdef LOC_TIMER
#include <Rcpp.h>
#include <rcpptimer.h>
#define LOC_TIMER_OBJ(name) Rcpp::Timer name(#name)
#define LOC_TIMER_CAT_(a, b) a##b
#define LOC_TIMER_CAT(a, b) LOC_TIMER_CAT_(a, b)
#define LOC_TIMER_SCOPED(t, lbl) \
    Rcpp::Timer::ScopedTimer LOC_TIMER_CAT(scpdtmr_, __LINE__)(t, lbl)
#define LOC_TIC(t, lbl) (t).tic(lbl)
#define LOC_TOC(t, lbl) (t).toc(lbl)
#define LOC_TIMER_PARAM(t) , Rcpp::Timer &t
#define LOC_TIMER_ARG(t) , t
// dump raw per-sample (duration_ns, tag) pairs to R's global environment as a
// data.frame. Companion to the auto-written summary; call before the Timer is
// destroyed.
#define LOC_TIMER_DUMP_RAW(t, r_name)                                       \
    do                                                                      \
    {                                                                       \
        Rcpp::Environment _env = Rcpp::Environment::global_env();           \
        _env.assign(r_name, Rcpp::DataFrame::create(                        \
                                Rcpp::Named("duration_ns") = (t).durations, \
                                Rcpp::Named("tag") = (t).tags));            \
    } while (0)
#else
#define LOC_TIMER_OBJ(name) ((void)0)
#define LOC_TIMER_SCOPED(t, lbl) ((void)0)
#define LOC_TIC(t, lbl) ((void)0)
#define LOC_TOC(t, lbl) ((void)0)
#define LOC_TIMER_PARAM(t)
#define LOC_TIMER_ARG(t)
#define LOC_TIMER_DUMP_RAW(t, r_name) ((void)0)
#endif
#endif // LOC_TIMER_H
