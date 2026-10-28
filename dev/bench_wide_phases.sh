#!/usr/bin/env bash
# Drive dev/bench_wide_phases.R: one R process per BLAS thread count, with
# OPENBLAS_NUM_THREADS set before R starts. Each process builds with
# load_all1(). R_LIBS_USER must point at the package library (Git Bash does
# not read Documents/.Renviron) and BENCH_OUT at an existing output directory.
# Run from the repository root: bash dev/bench_wide_phases.sh
set -eu
[ -f DESCRIPTION ] || { echo "run from the repository root"; exit 1; }
[ -n "${BENCH_OUT:-}" ] || { echo "set BENCH_OUT"; exit 1; }
R="Rscript dev/bench_wide_phases.R"

clean() { env -u OPENBLAS_NUM_THREADS -u OMP_NUM_THREADS "$@"; }

clean OPENBLAS_NUM_THREADS=1 $R run blas1
clean OPENBLAS_NUM_THREADS=4 $R run blas4
$R compare
