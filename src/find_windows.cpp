#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
List find_windows(const NumericVector x, double window, double overlap = 0.0)
{
  const int n = x.size();

  // preconditions guaranteed by the R wrapper:
  // - x is finite and sorted ascending
  // - window is finite and > 0
  // - overlap is finite and satisfies 0 <= overlap < window
  //
  // intended behavior:
  // - windows are coordinate-based, not count-based
  // - returned start/end indices are R-style 1-based inclusive indices
  // - this function uses half-open windows: [x[i], x[i] + window)
  // - if window is very small, singleton windows are valid output.
  //   R-side feasibility checks decide whether those windows are usable

  if (n == 0)
  {
    return List::create(
        Named("start") = IntegerVector(0),
        Named("end") = IntegerVector(0));
  }

  // starts/ends have size at most n
  IntegerVector starts(n);
  IntegerVector ends(n);

  int k = 0;
  int i = 0;

  while (i < n)
  {
    // current threshold is starting value + window
    // half-open interval: [x[i], x[i] + window)
    const double threshold = x[i] + window;

    // last index whose value is < threshold
    // lower_bound gives the first value >= threshold, so subtract 1
    const int j =
        static_cast<int>(
            std::lower_bound(x.begin() + i + 1, x.begin() + n, threshold) - x.begin()) -
        1;

    // note down start/end values, converting to R's 1-based indexing
    starts[k] = i + 1;
    ends[k] = j + 1;
    ++k;

    // if this window already covers the last element, we're done
    if (j == n - 1)
      break;

    // find next start: first index past i whose value >= cutoff.
    //
    // When overlap == 0, cutoff == threshold, so x[m] < threshold for all m <= j,
    // and lower_bound lands on j + 1 -- identical to non-overlapping behavior.
    //
    // When overlap > 0, this starts the next window at the first observed
    // coordinate in the trailing overlap region.
    const double cutoff = threshold - overlap;

    i = static_cast<int>(
        std::lower_bound(x.begin() + i + 1, x.begin() + j + 1, cutoff) - x.begin());
  }

  // k >= 1 here because n > 0
  starts = starts[Range(0, k - 1)];
  ends = ends[Range(0, k - 1)];

  return List::create(
      Named("start") = starts,
      Named("end") = ends);
}
