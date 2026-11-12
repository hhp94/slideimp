#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include "loc_timer.h"

// One repetition of the constrained NA sampler behind sample_na_loc().
//
// Returns a (total_na x 2) matrix of 1-based (row, col) positions, or a
// matrix with ZERO rows when no attempt filled every slot. The caller raises
// the error: a condition the user may want to catch is built in R, where it
// can carry a class and name the arguments to relax.
//
// `pool_idx_in` and `pool_cap` are parallel. `pool_cap(i)` is the most NAs
// column `pool_idx_in(i)` can take: its room under colmax, capped so that two
// observed values stay. sample_na_loc() passes only columns that can take the
// smallest demand and that hold two distinct observed values, so the skips
// below for a short or constant column are a guard, not the policy.
// [[Rcpp::export]]
arma::umat sample_each_rep_cpp(const arma::mat &obj,
                               const arma::uvec &pool_idx_in, // 0-based column indices available for sampling
                               const arma::uvec &na_per_col,  // how many NAs to place in each needed column, ascending
                               const arma::uvec &row_room,    // per-row budget: max additional NAs allowed
                               const arma::uvec &pool_cap,    // per-pool-column budget: max NAs the column can take
                               arma::uword max_attempts)
{
  LOC_TIMER_OBJ(sample_each_rep);
  LOC_TIC(sample_each_rep, "sample_each_rep_total");

  const arma::uword n_cols_needed = na_per_col.n_elem; // number of columns to fill
  const arma::uword total_na = arma::accu(na_per_col); // total NA cells to place
  const arma::uword n_rows = obj.n_rows;
  const arma::uword n_pool = pool_idx_in.n_elem;

  // reusable allocations across attempts
  arma::uvec left_temp(n_rows, arma::fill::none);
  arma::uvec row_out(total_na, arma::fill::zeros);
  arma::uvec col_out(total_na, arma::fill::zeros);
  arma::uvec attempt_row_room(row_room.n_elem);

  // Visiting order over the pool, as positions into `pool_idx_in`. It is
  // shuffled lazily, one swap per column visited, and it persists across
  // attempts: a uniform shuffle of any arrangement is uniform. An attempt
  // therefore costs the columns it visits, not the size of the pool.
  arma::uvec order(n_pool, arma::fill::none);
  for (arma::uword i = 0; i < n_pool; ++i)
  {
    order(i) = i;
  }

  for (arma::uword attempt = 0; attempt < max_attempts; ++attempt)
  {
    attempt_row_room = row_room;
    arma::uword written = 0;
    arma::uword slot = 0; // number of columns successfully filled so far
    for (arma::uword i = 0; i < n_pool; ++i)
    {
      if (slot >= n_cols_needed)
      {
        break; // all columns filled, stop scanning available pool
      }

      // Fisher-Yates step: position i takes a uniform draw from [i, n_pool)
      const arma::uword k = std::min(
          n_pool - 1, i + (arma::uword)(arma::randu() * (n_pool - i)));
      std::swap(order(i), order(k));
      const arma::uword pos = order(i);
      const arma::uword col_idx = pool_idx_in(pos);

      // NAs required for the current slot
      const arma::uword needed = na_per_col(slot);
      // skip column if its budget can't accommodate the needed count
      if (pool_cap(pos) < needed)
      {
        continue;
      }

      const double *candidate = obj.colptr(col_idx);

      // Split the observed rows by row budget.
      //  - A row with budget left can take an NA: it goes to `left_temp`.
      //  - A row at its budget stays observed whatever is drawn, so it
      //    protects the column for free. Track whether those rows already
      //    hold two distinct values.
      LOC_TIC(sample_each_rep, "scan_rows");
      arma::uword L = 0;
      bool fixed_any = false;
      bool fixed_two = false;
      double fixed_val = 0.0;
      for (arma::uword r = 0; r < n_rows; ++r)
      {
        const double v = candidate[r];
        if (!std::isfinite(v))
        {
          continue;
        }
        if (attempt_row_room(r) > 0)
        {
          left_temp(L++) = r;
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
          const arma::uword p0 = std::min(
              L - 1, (arma::uword)(arma::randu() * L));
          v0 = candidate[left_temp(p0)];
          left_temp(p0) = left_temp(L - 1);
          --L;
        }
        // second protected row: uniform over the rows whose value differs
        arma::uword n_diff = 0;
        for (arma::uword j = 0; j < L; ++j)
        {
          if (candidate[left_temp(j)] != v0)
          {
            ++n_diff;
          }
        }
        if (n_diff == 0)
        {
          LOC_TOC(sample_each_rep, "protect_pair");
          continue; // every observed value is identical
        }
        const arma::uword pick = std::min(
            n_diff - 1, (arma::uword)(arma::randu() * n_diff));
        arma::uword seen = 0;
        for (arma::uword j = 0; j < L; ++j)
        {
          if (candidate[left_temp(j)] == v0)
          {
            continue;
          }
          if (seen == pick)
          {
            left_temp(j) = left_temp(L - 1);
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
      for (arma::uword j = 0; j < needed; ++j)
      {
        const arma::uword k2 = std::min(
            L - 1, j + (arma::uword)(arma::randu() * (L - j)));
        std::swap(left_temp(j), left_temp(k2));
        const arma::uword r = left_temp(j);
        row_out(written + j) = r;
        col_out(written + j) = col_idx;
        attempt_row_room(r) -= 1;
      }
      written += needed;
      ++slot;
      LOC_TOC(sample_each_rep, "commit");
    }

    if (slot == n_cols_needed)
    {
      arma::umat result(total_na, 2);
      result.col(0) = row_out + 1; // convert to 1-based for R
      result.col(1) = col_out + 1;
      LOC_TOC(sample_each_rep, "sample_each_rep_total");
      return result;
    }
    // slot < n_cols_needed: reshuffle and retry
  }

  LOC_TOC(sample_each_rep, "sample_each_rep_total");
  return arma::umat(0, 2);
}
