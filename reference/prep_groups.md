# Prepare Groups for Imputation

Normalizes and validates a grouping specification for use with
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md).
Converts long-format or canonical list-column input into a validated
`slideimp_tbl`, enforcing set relationships, pruning dropped columns,
and optionally padding small groups.

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

  Character vector of column names from the data matrix (e.g.,
  `colnames(obj)`). Every element must appear in `group$feature` unless
  `allowed_unmapped = TRUE`.

- group:

  Specification of how features should be grouped for imputation.
  Accepts three formats:

  - `character`: (Requires the `{slideimp.extra}` package, available on
    GitHub). A single string naming a supported Illumina platform (e.g.,
    `"EPICv2"`, `"EPICv2_deduped"`). The manifest is fetched
    automatically. For installation instructions, see the package
    README. Supported platforms can be viewed via
    [`slideimp.extra::slideimp_arrays`](https://rdrr.io/pkg/slideimp.extra/man/slideimp_arrays.html).

  - `data.frame` (Long format):

    - `group`: Column identifying the group for each feature.

    - `feature`: Character column of individual feature names.

  - `data.frame` (List-column format):

    - `feature`: List-column of character vectors to impute (e.g., from
      `prep_groups()`).

    - `aux`: (Optional) List-column of auxiliary names used for context.

    - `parameters`: (Optional) List-column of group-specific parameter
      lists.

- subset:

  Character vector of feature names to impute (default `NULL` means
  impute all features). Must be a subset of `obj_cn` (`colnames(obj)`)
  and must appear in at least one group's `feature`. Features in a group
  but not in `subset` are demoted to auxiliary columns for that group.
  Groups left with zero features after demotion are dropped with a
  message.

- min_group_size:

  Integer or `NULL`. Minimum number of columns (features + aux) per
  group. Groups smaller than this are padded with randomly sampled
  columns from `obj`. Passed to `prep_groups()` internally.

- allow_unmapped:

  Logical. If `FALSE`, every column in `colnames(obj)` *must* appear in
  `group`. If `TRUE`, columns that have no group assignment are left
  untouched (neither imputed nor used as auxiliary columns) and a
  message is issued instead of an error.

- seed:

  Numeric or `NULL`. Random seed for reproducibility when padding for
  `min_group_size` and passed to
  [`pca_imp()`](https://hhp94.github.io/slideimp/reference/pca_imp.md).

## Value

A `data.frame` of class `slideimp_tbl` containing:

- `group`: Original group labels (if provided).

- `feature`: A list-column of character vectors (feature names).

- `aux`: A list-column of character vectors (auxiliary names).

- `parameters`: A list-column of per-group configuration lists.

## Details

### Set Validation

Let \\A\\ = `obj_cn` and \\B\\ = the union of all feature and auxiliary
names in `group`. The function enforces \\A \subseteq B\\: every column
in the matrix must appear somewhere in the manifest.

- `Pruning:` Elements in \\B\\ but not in \\A\\ (e.g., QC-dropped
  probes) are silently pruned from each group.

- `Dropping:` Groups left with zero features after pruning are removed
  entirely with a diagnostic message.

## See also

[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)
