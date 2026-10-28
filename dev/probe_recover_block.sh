#!/usr/bin/env bash
# Drive dev/probe_recover_block.R: one R process per BLAS thread count, with
# OPENBLAS_NUM_THREADS set before R starts. Base R only, so no package
# library is needed. BENCH_OUT must name an existing output directory.
# Run from the repository root: bash dev/probe_recover_block.sh
set -eu
[ -f DESCRIPTION ] || { echo "run from the repository root"; exit 1; }
[ -n "${BENCH_OUT:-}" ] || { echo "set BENCH_OUT"; exit 1; }
R="Rscript dev/probe_recover_block.R"

clean() { env -u OPENBLAS_NUM_THREADS -u OMP_NUM_THREADS "$@"; }

clean OPENBLAS_NUM_THREADS=1 $R run blas1
clean OPENBLAS_NUM_THREADS=4 $R run blas4
$R compare
