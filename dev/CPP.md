# C++ rules

Rules for code under `src/`. CLAUDE.md names each ban in one line and points
here; this file holds the full statement and the reasons, so a session that
edits `src/` reads it first. Like CLAUDE.md, it holds invariants only.

## No hand-rolled loop unrolling, no OpenMP, no pragma hints

Write the simple loop and let the compiler decide. SIMD is extremely
heterogeneous across the hardware this package will meet, and CRAN owns the
C++ settings - no `-march`, no per-arch flags, no build-time tuning. So the
code has to be transportable: what gets vectorized, and how, is the compiler's
call at whatever flags the target build uses, and the source's job is to be
the form auto-vectorizers recognize. That means contiguous access, plain
single-statement bodies, no per-iteration branches - and NOT manual
multi-accumulator reductions, strip-mining, or duplicated statements, which
hard-code one machine's answer, hide the loop's real shape from the next
compiler, and change floating-point results in the bargain.

Chunking a loop to amortize a CONTROL-FLOW check is not what this bans. An
early-abandon bound tested once per block rather than once per element - as in
`calc_distance_raw` (`src/impute_knn_brute.cpp`), which chunks over a
`GRAIN` constant so the distance loop can quit early without branching on
every row - is an algorithmic optimization that happens to look like a strip
mine. The ban is on chunking whose PURPOSE is to hand-roll vectorization. The
test is what the chunk boundary is for: if removing the block structure would
only change how the compiler vectorizes, it goes; if it would change how much
work the loop does, it stays.

OpenMP is banned as a parallelism mechanism, hint pragmas included (`#pragma
omp parallel`, `for`, `simd`, and friends). Parallelism is spent at the chain
level through RcppThread and mirai, where cores scale near 1:1, and a pragma's
effect is flag-shaped build state that no source file records. The OpenMP
flags in `src/Makevars` are there because the link line wants them, not as a
license to add pragmas. Accepted cost, stated so it is not relitigated per
loop: a strict-fp serial reduction may stay serial; that is the portable
answer, and a measured hotspot earns a BLAS call or a cleaner algorithm, not
an unroll.

## No thread pinning from C++

If threads get pinned, it is from R, via RhpcBLASctl, best-effort. The
threading landscape is too heterogeneous to control from native code -
OpenBLAS, MKL, BLIS, Accelerate and the OpenMP runtimes each have their own
knob, most are not present at compile time, and calling any vendor's API from
package C++ is a portability hazard the CRAN maintainer gets to keep. So: no
`openblas_set_num_threads`, no `mkl_set_num_threads`, no
`omp_set_num_threads`, no vendor headers, no dynamic-lookup tricks. Package
code that wants a pinned BLAS asks from R with
`RhpcBLASctl::blas_set_num_threads()` and accepts that on some machines it is
a no-op.

**Whether a pin actually took is verified by measurement - CPU time over wall
time - never by asking the library, which cannot be trusted to answer
honestly.** A local measurement harness that sets `OPENBLAS_NUM_THREADS` and
friends in a wrapper BEFORE R starts stays legal: it is not package code, and
the env-var route is the one interface every vendor honors.
