# Prepare Groups for Imputation

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

  Character vector of column names from the data matrix, usually
  `colnames(obj)`.

- group:

  Specification of how features should be grouped for imputation.
  Accepted formats are:

  - A character scalar naming a supported Illumina platform; see Note.

  - A long-format `data.frame` with columns `group` and `feature`.

  - A list-column `data.frame` with a `feature` list-column. Optional
    list-columns are `aux`, for auxiliary feature names, and
    `parameters`, for group-specific parameter lists.

- subset:

  Optional character vector of feature names to impute. If `NULL`, all
  grouped features are imputed. Features in a group but not in `subset`
  are demoted to auxiliary columns for that group. Groups left with zero
  features after demotion are dropped with a message.

- min_group_size:

  Integer or `NULL`. Minimum total number of columns per group, counting
  both features and auxiliary columns. Groups smaller than this are
  padded with randomly sampled columns from `obj_cn`. If `NULL` or `0`,
  no padding is performed.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` must appear in
  `group`. If `TRUE`, columns with no group assignment are left
  untouched and are not used as auxiliary columns.

- seed:

  Integer, numeric, or `NULL`. Random seed used when padding small
  groups.

## Value

A `data.frame` of class `slideimp_tbl` containing:

- `group`: original group labels, if provided, or sequential group
  labels.

- `feature`: a list-column of character vectors containing feature
  names.

- `aux`: a list-column of character vectors containing auxiliary names.

- `parameters`: a list-column of per-group configuration lists.

## Details

`prep_groups()` converts long-format or list-column group specifications
into a validated `slideimp_tbl`, enforces feature and auxiliary-column
set relationships, prunes dropped columns, and optionally pads small
groups.

Let \\A\\ be `obj_cn` and \\B\\ be the union of all feature and
auxiliary names in `group`. When `allow_unmapped = FALSE`, the function
enforces \\A \subseteq B\\: every matrix column must appear somewhere in
the grouping specification.

Elements in \\B\\ but not in \\A\\, such as previously dropped probes,
are pruned from each group. Groups left with zero features after pruning
or subset demotion are removed with a diagnostic message.

## See also

[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)

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
