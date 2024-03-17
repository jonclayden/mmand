#ifndef PARALLEL_H_
#define PARALLEL_H_

#if defined(HAVE_LIBDISPATCH) && defined(HAVE_BLOCKS)

#include <dispatch/dispatch.h>

#define PARALLEL_LOOP_START(i,n) \
    dispatch_apply(size_t(n), DISPATCH_APPLY_AUTO, ^(size_t i) {
#define PARALLEL_LOOP_END });

#elif _OPENMP

#include <omp.h>

#define PARALLEL_LOOP_START(i,n) \
    _Pragma("omp parallel for")  \
    for (size_t i=0; i<size_t(n); i++) {
#define PARALLEL_LOOP_END }

#else

#define PARALLEL_LOOP_START(i,n) \
    for (size_t i=0; i<size_t(n); i++) {
#define PARALLEL_LOOP_END }

#endif

#endif
