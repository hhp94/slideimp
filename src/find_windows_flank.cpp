#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List find_windows_flank(const NumericVector location,
                        const IntegerVector subset,
                        double window_size)
{
  const int n_sub = subset.size(); // number of centers we care about

  // preconditions guaranteed by the R wrapper:
  // - location is finite and sorted ascending
  // - subset is integer-like, has no missing values, and contains valid
  //   R-style 1-based indices into location
  // - window_size is finite and > 0
  //
  // intended behavior:
  // - returns one flanking window per subset value, in subset order
  // - returned start/end indices are R-style 1-based inclusive indices
  // - this function uses closed symmetric flanking windows:
  //   [center - window_size, center + window_size]
  // - if window_size is very small, singleton flank windows are valid output.
  //   R-side feasibility checks decide whether those windows are usable

  if (n_sub == 0)
  {
    return List::create(
        Named("start") = IntegerVector(0),
        Named("end") = IntegerVector(0),
        Named("subset_local") = IntegerVector(0));
  }

  IntegerVector start(n_sub);
  IntegerVector end(n_sub);
  IntegerVector subset_local(n_sub);

  for (int i = 0; i < n_sub; ++i)
  {
    // subset is supplied as R-style 1-based indices.
    // R wrapper guarantees bounds, so convert directly to 0-based.
    const int idx = subset[i];
    const int idx0 = idx - 1;

    const double center = location[idx0];

    const double lo = center - window_size;
    const double hi = center + window_size;

    // core flanking logic using binary search.
    // left flank: first position >= lo
    const int s = static_cast<int>(
        std::lower_bound(location.begin(), location.end(), lo) - location.begin());

    // right flank: last position <= hi.
    // upper_bound gives the first position that is too big,
    // so subtract 1 to get the last included position.
    const int e = static_cast<int>(
                      std::upper_bound(location.begin(), location.end(), hi) - location.begin()) -
                  1;

    // convert back to R's 1-based indexing for return.
    start[i] = s + 1;
    end[i] = e + 1;

    // local 1-based index of the center WITHIN this window.
    //
    // example:
    //  - returned global window is positions 2:4
    //  - center is original position 3
    //  - local index is 2, because 3 is the second element in 2:4.
    //
    // idx is 1-based and s is 0-based, so idx - s gives a 1-based local index.
    subset_local[i] = idx - s;
  }

  return List::create(
      Named("start") = start,
      Named("end") = end,
      Named("subset_local") = subset_local);
}
