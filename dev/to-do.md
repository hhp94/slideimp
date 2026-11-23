# To do

Work that is agreed on but not yet done. An entry stays here until the change
lands, and is then deleted, not marked done: git history records what
happened.

An entry says what to change and why, enough that a session with no memory of
the conversation can pick it up. If an entry turns out to be wrong, delete it
and say so in chat. A question that is settled rather than owed goes in
[DECISIONS.md](DECISIONS.md), not here.

## Kernel cleanup

Each item is a proposal with a measurement, put to the maintainer before it
is made; nothing here licenses a hand-rolled kernel (`dev/CPP.md`) or a call
that R's headers do not declare. `src/Spectra/` is out of scope
(`dev/vendored.md`). The gate of a change is `dev/gate/eigen_p4_gate.sh`: the
reference capture against HEAD and the suite, under CRAN's RcppEigen and the
candidate. The kernels' phase timers profile a release build through
`bash dev/gate/profile_run.sh <scratch> prep|install|run`.

- **LOBPCG slower than exact on `sim_mat(200, 5000)`.** The eigensolver takes
  126 ms a call under LOBPCG against 53 ms under the exact solver, 23 outer
  iterations each; on chr21 it is 14 ms against 29 ms. Not looked into.

## Release condition: pass on the RcppEigen candidate

slideimp 1.2.0 on CRAN does not link to RcppEigen; the sampler adds it to
`LinkingTo`, so the next release is the first that does. RcppEigen is moving
from Eigen 3.4 to Eigen 5, and its maintainer is patching by hand the reverse
dependencies that fail on the candidate
(https://github.com/RcppCore/RcppEigen/issues/151). A release that passes
only on CRAN's RcppEigen would add slideimp to that list.

Before any release that links to RcppEigen, on both CRAN's RcppEigen and the
candidate current on that day:

- The whole suite, `check_manual()` and `devtools::check()`. The candidate
  goes in its own private library, never over the CRAN one:
  `dev/gate/eigen_install_candidate.R` takes the library as an argument,
  installs from the RcppCore r-universe only when the candidate is missing
  or older, and prints the Eigen version its headers define. A run selects
  it with `R_LIBS` set before R starts; `dev/gate/eigen_p4_gate.sh`, stage
  suite, does that for the suite. "Passes" means the reverse-dependency
  run's standard, the suite and the check, not a compile. Record both
  versions tested in the release commit.
- The reference capture on the BLAS of the rhub Linux images, compared
  against one taken here. Until then the default tolerances of
  `dev/gate/ref_compare.R` are a reading from one BLAS, on one machine.
- The sampler's draws under fixed seeds, compared between the two builds:
  worst absolute and worst relative difference, and whether they are
  bitwise identical. If they move beyond what the `test-mipca-*` files
  tolerate, that is a finding to settle before the release.

The three runs are owed on release day, against the candidate current then.
The reference capture on the rhub images is to be taken again with the gate
as it now stands, which runs PCA cases at a fixed iteration count and does
not assert on `criterion_final`. The atlas image does not reproduce its own
capture, so its noise floor is not zero. Nothing is vendored from Eigen and
no version is pinned.

## Sampler

Open, the maintainer's:

- `pca_imp()` could warn when `solver = "lobpcg"` is forced at a high `ncp`.
- Experimental: investigate a new sampler for WGBS/EM-seq data that models
  the binomial C/T read counts at each CpG, in place of the beta values. It
  may never ship.
- The data transformations (`transform`, such as `"logit_biwhiten"`) are
  experimental and may be dropped. If a release keeps them, explain in its
  documentation what each one does.

