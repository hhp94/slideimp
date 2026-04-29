#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List find_windows(const NumericVector x, double window, double overlap = 0.0)
{
  // x is sorted
  // x, window, overlap are all positive
  // overlap must be < window
  int n = x.size();
  // starts/ends have size at most of n
  IntegerVector starts(n);
  IntegerVector ends(n);
  int k = 0;
  int i = 0;

  while (i < n)
  {
    // current threshold is starting value + window (half-open: [x[i], x[i]+window))
    double threshold = x[i] + window;
    int j = i;
    // scan forward, until reach a value of x >= threshold
    while (j + 1 < n && x[j + 1] < threshold)
    {
      ++j;
    }
    // note down start value and end values and iterate k
    starts[k] = i + 1;
    ends[k] = j + 1;
    ++k;

    // if this window already covers the last element, we're done
    if (j == n - 1)
      break;

    // find next start: first index past i whose value >= cutoff.
    // When overlap == 0, cutoff == threshold, so x[m] < threshold for all m <= j,
    // and we land on j + 1 - identical to non-overlapping behavior.
    double cutoff = threshold - overlap;
    int next_i = i + 1; // guarantee forward progress
    while (next_i <= j && x[next_i] < cutoff)
    {
      ++next_i;
    }
    i = next_i;
  }

  starts = starts[Range(0, k - 1)];
  ends = ends[Range(0, k - 1)];
  return List::create(
      Named("start") = starts,
      Named("end") = ends);
}
