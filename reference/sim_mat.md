# Simulate Matrix with Metadata

Generates a matrix of random normal data, then optionally scales values
between 0 and 1 column-wise. It also creates corresponding data frames
for feature (column) and sample (row) metadata and can optionally
introduce `NA` values into a specified proportion of rows. A correlation
between columns `rho` (before scaling) can be added.

## Usage

``` r
sim_mat(
  n = 100,
  p = 100,
  rho = 0.5,
  n_col_groups = 2,
  n_row_groups = 1,
  perc_total_na = 0.1,
  perc_col_na = 0.5,
  beta = TRUE
)
```

## Arguments

- n:

  An integer specifying the number of rows (samples). Default is `100`.

- p:

  An integer specifying the number of columns (features). Default is
  `100`.

- rho:

  Columns correlation before scaling (compound symmetry). Default is
  `0.5`.

- n_col_groups:

  An integer for the number of groups to assign to features/columns.
  Default is `2`.

- n_row_groups:

  An integer for the number of groups to assign to samples/rows. Default
  is `1`.

- perc_total_na:

  Proportion of all cells to set to NA. Default is `0.1`.

- perc_col_na:

  Proportion of columns across which those NAs are spread. Default is
  `0.5`.

- beta:

  If `TRUE` (default) scale values between 0 and 1 column wise.

## Value

An object of class `slideimp_sim`. This is a list containing:

- `input`: A numeric matrix of dimension \\n \times p\\ containing the
  simulated values and injected `NA`s.

- `col_group`: A data frame with \$p\$ rows mapping each `feature` to a
  `group`.

- `row_group`: A data frame with \$n\$ rows mapping each `sample` to a
  `group`.

## Examples

``` r
set.seed(123)
sim_data <- sim_mat(n = 50, p = 10, rho = 0.5)
sim_data
#> $col_group (2 column groups)
#>    feature  group
#> 1 feature1 group1
#> 2 feature2 group1
#> 3 feature3 group1
#> 4 feature4 group1
#> 5 feature5 group2
#> 6 feature6 group2
#> 
#> $row_group (1 row groups)
#>    sample  group
#> 1 sample1 group1
#> 2 sample2 group1
#> 3 sample3 group1
#> 4 sample4 group1
#> 5 sample5 group1
#> 6 sample6 group1
#> 
#> $input (50 x 10)
#>          feature1  feature2  feature3  feature4  feature5  feature6
#> sample1 0.3855027 0.2859279        NA 0.7034804 0.3949006 0.3390209
#> sample2 0.3939132 0.5091468 0.5414625 0.6130915 0.4153076 0.3899470
#> sample3 0.7020675 0.7302534 0.7618661 0.6474409 0.6996583 0.6687105
#> sample4 0.6887438 0.4568955 0.3007332 0.5369357 0.5503469 0.3900968
#> sample5 0.4220865 0.3630904 0.4552229        NA 0.7723463 0.5073267
#> sample6 1.0000000 0.7918423 0.6874920 0.6385426 0.7579935 0.9167009
#> # Showing [1:6, 1:6] of full matrix
```
