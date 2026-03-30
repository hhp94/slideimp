#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// find_windows_flank <- function(location, subset, window_size) {
//   n <- length(location)
//   n_sub <- length(subset)
//   start <- integer(n_sub)
//   end   <- integer(n_sub)
//   subset_local <- vector("list", n_sub)
//
//   for (i in seq_len(n_sub)) {
//     center <- location[subset[i]]
//     lo <- center - window_size
//     hi <- center + window_size
//     # binary search via findInterval
//     s <- findInterval(lo, location, left.open = FALSE)
//     s <- max(s, 1L)
//     # findInterval returns the index of the last value <= lo,
//     # so if location[s] < lo we need s+1
//     if (location[s] < lo) s <- s + 1L
//     e <- findInterval(hi, location)
//     if (e < 1L) e <- 1L
//     start[i] <- s
//     end[i]   <- e
//     subset_local[[i]] <- subset[i] - s + 1L
//   }
//
//   list(start = start, end = end, subset_local = subset_local)
// }

// Find the indices in `location` that flank each position in `subset`
// by a symmetric window of size `window_size` on each side.
// Returns:
//   - start / end : 1-based indices of the window in `location`
//   - subset_local: 1-based index of the center *inside* its own window
// [[Rcpp::export]]
List find_windows_flank_cpp(NumericVector location,
                            IntegerVector subset,
                            double window_size)
{
  int n = location.size();   // total number of positions
  int n_sub  = subset.size(); // number of centers we care about

  IntegerVector start(n_sub);
  IntegerVector end(n_sub);
  IntegerVector subset_local(n_sub);

  // location MUST be sorted ascending
  for (int i = 0; i < n_sub; i++)
  {
    // convert to 0-based index
    double center = location[subset[i] - 1];

    double lo = center - window_size;
    double hi = center + window_size;

    // Core flanking logic (binary search)
    // left flank: first position >= lo
    int s = (int)(std::lower_bound(location.begin(), location.end(), lo) - location.begin());
    s = std::clamp(s, 0, n - 1);
    // right flank: last position <= hi. std::upper_bound gives the first thing that is too big,
    // so we want to - 1 it.
    int e = (int)(std::upper_bound(location.begin(), location.end(), hi) - location.begin()) - 1;
    e = std::clamp(e, 0, n - 1);
    // Convert back to R's 1-based indexing for return
    start[i] = s + 1;
    end[i] = e + 1;
    // Local 1-based index of the center WITHIN this window. For example, for a window of 2 to 4,
    // center at 3, the local index is 3 - 2 = 1.
    subset_local[i] = subset[i] - s;
  }

  return List::create(
    Named("start") = start,
    Named("end") = end,
    Named("subset_local") = subset_local
  );
}
