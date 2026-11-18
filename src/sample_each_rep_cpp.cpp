#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include "loc_timer.h"

// One repetition of the constrained NA sampler behind sample_na_loc().
//
// Returns a (total_na x 2) integer matrix of 1-based (row, col) positions, or
// a matrix with ZERO rows when no attempt filled every slot. The caller
// raises the error: a condition the user may want to catch is built in R,
// where it can carry a class and name the arguments to relax.
//
// `pool_idx_in` and `pool_cap` are parallel. `pool_cap[i]` is the most NAs
// column `pool_idx_in[i]` can take: its room under colmax, capped so that two
// observed values stay. sample_na_loc() passes only columns that can take the
// smallest demand and that hold two distinct observed values, so the skips
// below for a short or constant column are a guard, not the policy.
//
// Every uniform draw is R::runif(0, 1), from R's generator, so set.seed()
// fixes the positions.

// Refuses a vector holding a value below 0 (NA_integer_ included) or at or
// above `bound`, so every value can index or count safely as std::size_t.
static void check_range(const Rcpp::IntegerVector &v,
                        const long long bound,
                        const char *name)
{
  for (R_xlen_t i = 0; i < v.size(); ++i)
  {
    if (v[i] < 0 || v[i] >= bound)
    {
      Rcpp::stop(std::string(name) + " holds " +
                 (v[i] == NA_INTEGER ? std::string("NA")
                                     : std::to_string(v[i])) +
                 " at position " + std::to_string(i + 1) +
                 ", outside 0 .. " + std::to_string(bound - 1) + ".");
    }
  }
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix sample_each_rep_cpp(const Rcpp::NumericMatrix &obj,
                                        const Rcpp::IntegerVector &pool_idx_in, // 0-based column indices available for sampling
                                        const Rcpp::IntegerVector &na_per_col,  // how many NAs to place in each needed column, ascending
                                        const Rcpp::IntegerVector &row_room,    // per-row budget: max additional NAs allowed
                                        const Rcpp::IntegerVector &pool_cap,    // per-pool-column budget: max NAs the column can take
                                        const int max_attempts)
{
  LOC_TIMER_OBJ(sample_each_rep);
  LOC_TIC(sample_each_rep, "sample_each_rep_total");

  const std::size_t n_rows = static_cast<std::size_t>(obj.nrow());
  const long long int_max = std::numeric_limits<int>::max() + 1LL;

  // The kernel indexes with these values, so it checks them, not the caller.
  if (row_room.size() != obj.nrow())
  {
    Rcpp::stop("row_room has length " + std::to_string(row_room.size()) +
               "; obj has " + std::to_string(obj.nrow()) + " rows.");
  }
  if (pool_cap.size() != pool_idx_in.size())
  {
    Rcpp::stop("pool_cap and pool_idx_in differ in length.");
  }
  if (max_attempts < 0)
  {
    Rcpp::stop("max_attempts must be 0 or more.");
  }
  check_range(pool_idx_in, obj.ncol(), "pool_idx_in");
  check_range(na_per_col, int_max, "na_per_col");
  check_range(row_room, int_max, "row_room");
  check_range(pool_cap, int_max, "pool_cap");

  const int *pool = pool_idx_in.begin();
  const int *need = na_per_col.begin();
  const int *cap = pool_cap.begin();
  const std::size_t n_cols_needed = static_cast<std::size_t>(na_per_col.size());
  const std::size_t n_pool = static_cast<std::size_t>(pool_idx_in.size());
  const double *x = obj.begin();

  // total NA cells to place: at most one per cell of `obj`, but summed in
  // std::size_t so a caller's error cannot overflow it before it is refused
  std::size_t total_na = 0;
  for (std::size_t i = 0; i < n_cols_needed; ++i)
  {
    total_na += static_cast<std::size_t>(need[i]);
  }
  if (total_na > static_cast<std::size_t>(std::numeric_limits<int>::max()))
  {
    Rcpp::stop("na_per_col sums past the rows an R matrix can hold.");
  }

  // reusable allocations across attempts
  std::vector<std::size_t> left_temp(n_rows);
  std::vector<int> row_out(total_na, 0);
  std::vector<int> col_out(total_na, 0);
  std::vector<int> attempt_row_room(n_rows);

  // Visiting order over the pool, as positions into `pool_idx_in`. It is
  // shuffled lazily, one swap per column visited, and it persists across
  // attempts: a uniform shuffle of any arrangement is uniform. An attempt
  // therefore costs the columns it visits, not the size of the pool.
  std::vector<std::size_t> order(n_pool);
  for (std::size_t i = 0; i < n_pool; ++i)
  {
    order[i] = i;
  }

  for (int attempt = 0; attempt < max_attempts; ++attempt)
  {
    std::copy(row_room.begin(), row_room.end(), attempt_row_room.begin());
    std::size_t written = 0;
    std::size_t slot = 0; // number of columns successfully filled so far
    for (std::size_t i = 0; i < n_pool; ++i)
    {
      if (slot >= n_cols_needed)
      {
        break; // all columns filled, stop scanning available pool
      }

      // Fisher-Yates step: position i takes a uniform draw from [i, n_pool)
      const std::size_t k = std::min(
          n_pool - 1, i + (std::size_t)(R::runif(0, 1) * (n_pool - i)));
      std::swap(order[i], order[k]);
      const std::size_t pos = order[i];
      const int col_idx = pool[pos];

      // NAs required for the current slot
      const std::size_t needed = static_cast<std::size_t>(need[slot]);
      // skip column if its budget can't accommodate the needed count
      if (static_cast<std::size_t>(cap[pos]) < needed)
      {
        continue;
      }

      const double *candidate = x + static_cast<std::size_t>(col_idx) * n_rows;

      // Split the observed rows by row budget.
      //  - A row with budget left can take an NA: it goes to `left_temp`.
      //  - A row at its budget stays observed whatever is drawn, so it
      //    protects the column for free. Track whether those rows already
      //    hold two distinct values.
      LOC_TIC(sample_each_rep, "scan_rows");
      std::size_t L = 0;
      bool fixed_any = false;
      bool fixed_two = false;
      double fixed_val = 0.0;
      for (std::size_t r = 0; r < n_rows; ++r)
      {
        const double v = candidate[r];
        if (!std::isfinite(v))
        {
          continue;
        }
        if (attempt_row_room[r] > 0)
        {
          left_temp[L++] = r;
        }
        else if (!fixed_any)
        {
          fixed_any = true;
          fixed_val = v;
        }
        else if (v != fixed_val)
        {
          fixed_two = true;
        }
      }
      LOC_TOC(sample_each_rep, "scan_rows");

      // Protect two rows with different values, so the column keeps its
      // variance. Rows at their budget are used first, because protecting a
      // row that could take an NA spends a cell the column may need. Rows
      // protected out of `left_temp` are removed from it.
      LOC_TIC(sample_each_rep, "protect_pair");
      if (!fixed_two)
      {
        double v0 = fixed_val;
        if (!fixed_any)
        {
          if (L == 0)
          {
            LOC_TOC(sample_each_rep, "protect_pair");
            continue; // no observed value at all
          }
          // first protected row: uniform over the observed rows
          const std::size_t p0 = std::min(
              L - 1, (std::size_t)(R::runif(0, 1) * L));
          v0 = candidate[left_temp[p0]];
          left_temp[p0] = left_temp[L - 1];
          --L;
        }
        // second protected row: uniform over the rows whose value differs
        std::size_t n_diff = 0;
        for (std::size_t j = 0; j < L; ++j)
        {
          if (candidate[left_temp[j]] != v0)
          {
            ++n_diff;
          }
        }
        if (n_diff == 0)
        {
          LOC_TOC(sample_each_rep, "protect_pair");
          continue; // every observed value is identical
        }
        const std::size_t pick = std::min(
            n_diff - 1, (std::size_t)(R::runif(0, 1) * n_diff));
        std::size_t seen = 0;
        for (std::size_t j = 0; j < L; ++j)
        {
          if (candidate[left_temp[j]] == v0)
          {
            continue;
          }
          if (seen == pick)
          {
            left_temp[j] = left_temp[L - 1];
            --L;
            break;
          }
          ++seen;
        }
      }
      LOC_TOC(sample_each_rep, "protect_pair");

      // randomly select exactly `needed` rows from the eligible set
      // `left_temp` (size L) without replacement (`sample.int(L, needed)` in R)
      // for j = 0, 1, ..., needed-1
      //  1. k <- j + runif{0, ..., L-j-1} // uniform in [j, L-1]
      //  2. swap(left_temp[j], left_temp[k])
      //  3. commit left_temp[j] as the selected row
      if (L < needed)
      {
        continue;
      }
      LOC_TIC(sample_each_rep, "commit");
      for (std::size_t j = 0; j < needed; ++j)
      {
        const std::size_t k2 = std::min(
            L - 1, j + (std::size_t)(R::runif(0, 1) * (L - j)));
        std::swap(left_temp[j], left_temp[k2]);
        const std::size_t r = left_temp[j];
        row_out[written + j] = static_cast<int>(r);
        col_out[written + j] = col_idx;
        attempt_row_room[r] -= 1;
      }
      written += needed;
      ++slot;
      LOC_TOC(sample_each_rep, "commit");
    }

    if (slot == n_cols_needed)
    {
      Rcpp::IntegerMatrix result(Rcpp::no_init_matrix(
          static_cast<int>(total_na), 2));
      int *res = result.begin();
      for (std::size_t j = 0; j < total_na; ++j)
      {
        res[j] = row_out[j] + 1; // 1-based for R
        res[total_na + j] = col_out[j] + 1;
      }
      LOC_TOC(sample_each_rep, "sample_each_rep_total");
      return result;
    }
    // slot < n_cols_needed: reshuffle and retry
  }

  LOC_TOC(sample_each_rep, "sample_each_rep_total");
  return Rcpp::IntegerMatrix(0, 2);
}
