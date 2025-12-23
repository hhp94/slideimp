## Resubmission
This is a resubmission. In this version I have:

* Omitted the single quotes around imputePCA() in my DESCRIPTION
* Removed the .Rd file of the internal only function inject_na() as well as the examples

Other than the changes above that directly address the comments, I have also:

* Updated my DESCRIPTION. Removed references to "intensive longitudinal data" (not applicable). Also applied this to README and vignettes.
* Added group_features(). This helper function tremendously increases user-friendliness of the most important function of my package, group_imp(). This feature added one dependency, the `collapse` package.
* Allow pca_imp() to automatically scale row weights by the number of missing values in the row.
* Removed the development only find_knn_brute() function to reduce compile time.
* Fixed one bug in group_imp() where the row.w argument wasn't included.
* Fixed one bug in qrSVD_cpp() function to warn instead of stop if svd failed to converge.
* Improved examples and docs. Fixed many typos.

Happy Holidays! Thank you very much for what you do. Cheers!

## R CMD check results
0 errors | 0 warnings | 1 note

* New submission.
