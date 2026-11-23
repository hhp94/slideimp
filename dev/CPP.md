# C++ rules

Rules for code under `src/`. `.claude/CLAUDE.md` names each ban in one line and
points here; this file holds the full statement and the reasons. Like
`.claude/CLAUDE.md`, it holds invariants only.

## No hand-rolled loop unrolling, no OpenMP, no pragma hints

Write the simple loop and let the compiler decide. SIMD differs widely across
the hardware this package meets, and CRAN owns the C++ settings - no `-march`,
no per-arch flags, no build-time tuning. What gets vectorized, and how, is the
compiler's call at whatever flags the target build uses, and the source's job
is to be the form auto-vectorizers recognize: contiguous access, plain
single-statement bodies, no per-iteration branches. NOT manual
multi-accumulator reductions, strip-mining, or duplicated statements, which
hard-code one machine's answer, hide the loop's real shape from the next
compiler, and change floating-point results.

Chunking a loop to amortize a CONTROL-FLOW check is not what this bans. An
early-abandon bound tested once per block rather than once per element is an
algorithmic optimization that happens to look like a strip mine.
`calc_distance_raw` (`src/knn_kernel.cpp`) is one: it chunks over a `GRAIN`
constant so the distance loop can quit early without branching on every row. The
test is what the chunk boundary is for: if removing the block structure would
only change how the compiler vectorizes, it goes; if it would change how much
work the loop does, it stays.

OpenMP is banned as a parallelism mechanism, hint pragmas included (`#pragma
omp parallel`, `for`, `simd`, and friends). Parallelism is spent at the chain
level through RcppThread and mirai, and what a pragma does depends on build
flags that no source file records. Neither `src/Makevars` nor
`src/Makevars.win` passes the OpenMP flags, so `_OPENMP` is not defined and no
header compiles its OpenMP path, Eigen's included.

A strict-fp serial reduction may stay serial: that is the portable answer. A
measured hotspot earns a BLAS call or a cleaner algorithm, not an unroll.

## No thread pinning from C++

If threads get pinned, it is from R, via RhpcBLASctl, best-effort. OpenBLAS,
MKL, BLIS, Accelerate and the OpenMP runtimes each have their own knob, most
are not present at compile time, and calling any vendor's API from package C++
is a portability hazard. So: no `openblas_set_num_threads`, no
`mkl_set_num_threads`, no `omp_set_num_threads`, no vendor headers, no
dynamic-lookup tricks. Package code that wants a pinned BLAS asks from R with
`RhpcBLASctl::blas_set_num_threads()` and accepts that on some machines it is
a no-op.

**Whether a pin took is verified by measurement - CPU time over wall time -
never by asking the library.** A local measurement harness that sets
`OPENBLAS_NUM_THREADS` and friends in a wrapper BEFORE R starts stays legal:
it is not package code, and the env-var route is the one interface every
vendor honors.

## All dimnames are handled in R

C++ under `src/` neither sets nor reads dimnames. A matrix goes back to R with
none, and the R function that made the `.Call` names it right after the call,
with the primitive on the list element:

```r
dimnames(res$loc) <- list(NULL, c("row", "col"))
```

So no `Rcpp::rownames(x) =`, `Rcpp::colnames(x) =`, `.attr("dimnames") =` or
`Rf_setAttrib(x, R_DimNamesSymbol, ...)`, and no `Rcpp::rownames(x)` or
`Rcpp::colnames(x)` read either. The Rcpp assignments run through
`Rcpp::internal::DimNameProxy::assign` (`Rcpp/vector/DimNameProxy.h`), which
rchk reports as an unprotected variable: an error in Rcpp's own header, but
reported against slideimp's library and counted in its rchk result. The other
routes are banned with it so that every name a result carries is set in one
place, the R side, where a reader of the R code sees it. A kernel that needs to
know which column is which takes the indices as arguments, never the names.

## No Rcpp::List built from a SEXP

An exported function takes scalars, Rcpp's atomic vectors and matrices, or a
`SEXP` it checks; it never takes `Rcpp::List`. A list from R arrives as
`SEXP`, is refused unless `TYPEOF(x)` is `VECSXP`, and is read with
`Rf_xlength()` and `VECTOR_ELT()`; a list nested in it is read the same way.
An `Rcpp::List` is made only by size, `Rcpp::List x(n)`, or with
`Rcpp::List::create()`, to be filled and returned. The same holds for
`Rcpp::DataFrame`, `Rcpp::ExpressionVector`, `Rcpp::Language`,
`Rcpp::Pairlist` and `Rcpp::Formula` built from a `SEXP`, and for
`Rcpp::Environment::namespace_env()`. No `Storage::set__()` to put a `SEXP`
into an `Rcpp::List` either: it is Rcpp's internal storage policy, not its
API.

Every way to build an `Rcpp::List` from a `SEXP` goes through
`r_cast<VECSXP>`: the `Vector(SEXP)` constructor, a proxy
(`Rcpp::List x = y[i]`), assignment, `Rcpp::as<Rcpp::List>()`, and the
`input_parameter< Rcpp::List >` that `compileAttributes()` writes for an
exported argument. The cast falls back to R's `as.list()` through
`internal::convert_using_rfunction()` (`Rcpp/r_cast.h`), which holds the
result in an `Rcpp::Armor<SEXP>` (`Rcpp/protection/Armor.h`), and rchk
reports that Armor as errors against slideimp's library, as with the dimnames
proxy above. The cast is compiled in whether or not it runs, so a type check
in front of it does not remove it. The casts to `EXPRSXP`, `LANGSXP` and
`LISTSXP` use the same fallback, a `DataFrame` or `Formula` built from a
`SEXP` reaches it too, and `namespace_env()` holds an Armor of its own. Building by size, `create()`,
the copy constructor, `x[i] = value` and returning an `Rcpp::List` do not
reach it.

Code compiled only under `-DLOC_TIMER` is exempt, such as
`LOC_TIMER_DUMP_RAW` in `src/loc_timer.h`, which calls
`Rcpp::DataFrame::create()`: it never reaches a release build.

`dev/check_rcpp_list.R` builds the objects and fails when one holds a symbol
of the fallback; its header says which build it reads and why.
