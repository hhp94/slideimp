# Resolve a Group Specification to a Data Frame

Convert a group specification to the canonical data-frame form expected
by
[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md).
This S3 generic is exported so that extension packages, such as
`slideimp.extra`, can register additional methods.

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

  A group specification. `slideimp` provides a method for `data.frame`
  objects. The optional `slideimp.extra` package provides a method for
  character platform names.

## Note

This is primarily an extension hook for `slideimp.extra`.

## Examples

``` r
df <- data.frame(feature = c("cg1", "cg2"), group = c(1, 1))
slideimp_resolve_group(df)
#>   feature group
#> 1     cg1     1
#> 2     cg2     1
```
