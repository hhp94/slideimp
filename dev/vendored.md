# Vendored files

`src/Spectra/` is the one vendored code in slideimp: the Spectra eigensolver
library (https://github.com/yixuan/spectra), MPL-2.0, which `src/mipca.cpp`
includes by a quoted path. `DESCRIPTION` credits it, and its `MANIFEST`
lists the files and the upstream commit they were taken from, with one git
blob hash per file.

Nobody edits it by hand. A re-vendor is the whole directory, written by
`Rscript dev/tool/vendor_spectra.R` from a clone of Spectra at the commit
pinned in that script; the script checks that the include graph from the two
headers `src/mipca.cpp` uses reaches exactly the listed files, and that every
written file hashes to its upstream blob. A patch, if one is ever needed, is
made by the script and named in its header and in the manifest.
