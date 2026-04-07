#ifndef LOC_TIMER_H
#define LOC_TIMER_H

#ifdef LOC_TIMER
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
#else
#define LOC_TIMER_OBJ(name) ((void)0)
#define LOC_TIMER_SCOPED(t, lbl) ((void)0)
#define LOC_TIC(t, lbl) ((void)0)
#define LOC_TOC(t, lbl) ((void)0)
#define LOC_TIMER_PARAM(t)
#define LOC_TIMER_ARG(t)
#endif

#endif // LOC_TIMER_H
