# Resolve a group specification to a data.frame

S3 generic for converting various group specifications into the
canonical data.frame form expected by
[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md).
This generic exists only to allow `slideimp.extra` to register the
character method.

## Usage

``` r
slideimp_resolve_group(x)

# S3 method for class 'data.frame'
slideimp_resolve_group(x)

# Default S3 method
slideimp_resolve_group(x)
```

## Arguments

- x:

  A group specification. The base package provides a method for
  `data.frame`. `slideimp.extra` provides `character` for chip-name
  lookup.

## Value

A data.frame with at least a `feature` column, suitable for passing to
[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md).

## Examples

``` r
df <- data.frame(feature = c("cg1", "cg2"), group = c(1, 1))
slideimp_resolve_group(df)
#>   feature group
#> 1     cg1     1
#> 2     cg2     1
```
