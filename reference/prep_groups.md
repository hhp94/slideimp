# Group Preparation for Imputation

Normalize and validate a grouping specification for use with
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).

## Usage

``` r
prep_groups(
  obj_cn,
  group,
  subset = NULL,
  min_group_size = 0,
  allow_unmapped = FALSE,
  seed = NULL
)
```

## Arguments

- obj_cn:

  Character. The column names of the data matrix, usually
  `colnames(obj)`.

- group:

  A data frame or a character scalar. The feature groups, in one of
  these forms:

  - A character scalar naming a supported Illumina platform, such as
    `"EPICv2"`. This needs the `slideimp.extra` package.

  - A long-format data frame with columns `group` and `feature`.

  - A list-column data frame with a `feature` list-column. Optional
    list-columns are `aux`, for auxiliary feature names, and
    `parameters`, for group-specific parameter lists.

- subset:

  Optional character vector. The features to impute, by name, or `NULL`
  for every grouped feature. A group's other features become its
  auxiliary columns, and a group left with no feature is dropped with a
  message.

- min_group_size:

  Integer. Minimum number of columns per group, features and auxiliary
  columns together. Groups smaller than this are padded with columns
  sampled at random from the other grouped columns in `obj_cn`. If `0`,
  no padding is performed.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` must appear in
  `group`. If `TRUE`, columns with no group assignment are left
  untouched and are not used as auxiliary columns.

- seed:

  Integer, numeric, or `NULL`. Random seed used when padding small
  groups.

## Value

A data frame of class `slideimp_tbl` containing:

- `group`: original group labels, if provided, or sequential group
  labels.

- `feature`: a list-column of character vectors containing feature
  names.

- `aux`: a list-column of character vectors containing auxiliary names.

- `parameters`: a list-column of per-group configuration lists.

## Details

`prep_groups()` converts a long-format or list-column `group` to a
`slideimp_tbl`, and pads small groups when `min_group_size` is above
`0`.

Let \\A\\ be `obj_cn` and \\B\\ be the union of all feature and
auxiliary names in `group`. When `allow_unmapped = FALSE`,
`prep_groups()` requires \\A \subseteq B\\: every name in `obj_cn` must
appear in `group`.

Names in \\B\\ but not in \\A\\, such as probes dropped in quality
control, are removed from each group. Groups left with no features are
dropped with a message.

## See also

[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
to impute with the prepared groups.

## Examples

``` r
sim <- sim_mat(n = 10, p = 20)
prepped <- prep_groups(colnames(sim$input), sim$col_group)
prepped
#> # slideimp table: 2 x 4
#>   group          feature             aux parameters
#>  group1 <character [12]> <character [0]> <list [0]>
#>  group2  <character [8]> <character [0]> <list [0]>
```
