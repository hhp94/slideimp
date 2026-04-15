#include <RcppArmadillo.h>
#include "loc_timer.h"

// [[Rcpp::export]]
arma::umat sample_each_rep_cpp(const arma::mat &obj,
                               const arma::uvec &pool_idx_in, // 0-based column indices available for sampling
                               const arma::uvec &na_per_col,  // how many NAs to place in each needed column
                               const arma::uvec &row_room, // per-row budget: max additional NAs allowed
                               const arma::uvec &col_room,    // per-column budget: max NAs allowed
                               arma::uword max_attempts)
{
  LOC_TIMER_OBJ(sample_each_rep);
  LOC_TIC(sample_each_rep, "sample_each_rep_total");

  const arma::uword n_cols_needed = na_per_col.n_elem; // number of columns to fill
  const arma::uword total_na = arma::accu(na_per_col); // total NA cells to place
  const arma::uword n_rows = obj.n_rows;

  // reusable allocations across attempts
  arma::uvec left_temp(n_rows, arma::fill::none);
  arma::uvec row_out(total_na, arma::fill::zeros);
  arma::uvec col_out(total_na, arma::fill::zeros);
  arma::uvec attempt_row_room(row_room.n_elem);

  for (arma::uword attempt = 0; attempt < max_attempts; ++attempt)
  {
    LOC_TIC(sample_each_rep, "attempt_setup");
    // shuffle column order each attempt to escape bad orderings `pool_idx[sample.int(length(pool_idx))]`
    arma::uvec shuffled = arma::shuffle(pool_idx_in);
    attempt_row_room = row_room;
    arma::uword written = 0;
    arma::uword slot = 0; // number of columns successfully filled so far
    LOC_TOC(sample_each_rep, "attempt_setup");
    for (arma::uword i = 0; i < shuffled.n_elem; ++i) // for (col_idx in shuffled)
    {
      if (slot >= n_cols_needed)
      {
        break; // all columns filled, stop scanning available pool
      }

      const arma::uword col_idx = shuffled(i);
      // NAs required for the current slot
      const arma::uword needed = na_per_col(slot);
      // skip column if its budget can't accommodate the needed count
      if (col_room(col_idx) < needed)
      {
        continue;
      }

      const double *candidate = obj.colptr(col_idx);

      LOC_TIC(sample_each_rep, "find_finite");
      arma::uvec obs_row = arma::find_finite(obj.col(col_idx)); // which(!is.na(candidate))
      LOC_TOC(sample_each_rep, "find_finite");

      if (obs_row.n_elem < 2)
      {
        continue;
      }
      // Protect a heterogeneous pair of rows
      // Pick keep0 at random, then keep1 from rows with a different value.
      // Skip the column if all observed values are identical (ZV).
      LOC_TIC(sample_each_rep, "protect_pair");
      LOC_TIC(sample_each_rep, "protect_pair_v0");
      const arma::uword p0_idx = arma::randi<arma::uword>(
          arma::distr_param(0, (arma::sword)obs_row.n_elem - 1));
      const arma::uword keep0 = obs_row(p0_idx);
      const double v0 = candidate[keep0];
      LOC_TOC(sample_each_rep, "protect_pair_v0");

      LOC_TIC(sample_each_rep, "protect_pair_v1");
      // pass 1: count rows whose value differs from v0.
      arma::uword n_diff = 0;
      for (arma::uword j = 0; j < obs_row.n_elem; ++j)
      {
        if (candidate[obs_row(j)] != v0)
        {
          ++n_diff;
        }
      }
      if (n_diff == 0)
      {
        LOC_TOC(sample_each_rep, "protect_pair_v1");
        LOC_TOC(sample_each_rep, "protect_pair");
        continue; // column is ZV among observed values
      }
      // pass 2: pick the pick-th differing row uniformly.
      const arma::uword pick = (arma::uword)(arma::randu() * n_diff);
      arma::uword keep1 = 0;
      arma::uword seen = 0;
      for (arma::uword j = 0; j < obs_row.n_elem; ++j)
      {
        const arma::uword r = obs_row(j);
        if (candidate[r] == v0)
        {
          continue;
        }
        if (seen == pick)
        {
          keep1 = r;
          break;
        }
        ++seen;
      }
      LOC_TOC(sample_each_rep, "protect_pair_v1");
      LOC_TOC(sample_each_rep, "protect_pair");
      // build eligible rows (exclude protected + row-budget exhausted)
      // R does this in two vectorised steps:
      //  1. left_over <- setdiff(obs_row, kept_global)
      //  2. left_over <- left_over[attempt_row_room[left_over] > 0L]
      // Here, we fuse both filters into a single loop, push_back-ing only rows that pass both guards.
      // inside the loop:
      LOC_TIC(sample_each_rep, "build_leftover");
      arma::uword L = 0;
      for (arma::uword j = 0; j < obs_row.n_elem; ++j)
      {
        const arma::uword r = obs_row(j);
        if (r != keep0 && r != keep1 && attempt_row_room(r) > 0)
        {
          left_temp(L++) = r;
        }
      }
      LOC_TOC(sample_each_rep, "build_leftover");
      // randomly select exactly `needed` rows from the eligible set
      // `left_vec` (size L) without replacement (`sample.int(length(left_over), needed)`in R)
      // for j = 0, 1, ..., needed-1
      //  1. k <- j + runif{0, ..., L-j-1} // uniform in [j, L-1]
      //  2. swap(left_vec[j], left_vec[k])
      //  3. commit left_vec[j] as the selected row
      if (L < needed)
      {
        continue;
      }
      LOC_TIC(sample_each_rep, "commit");
      for (arma::uword j = 0; j < needed; ++j)
      {
        const arma::uword k = j + (arma::uword)(arma::randu() * (L - j));
        std::swap(left_temp(j), left_temp(k));
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
    // slot < n_cols_needed: shuffle and retry
  }

  LOC_TOC(sample_each_rep, "sample_each_rep_total");
  Rcpp::stop("Failed to sample NA locations after %d attempts.",
             static_cast<int>(max_attempts));
}
