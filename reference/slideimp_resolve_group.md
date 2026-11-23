# Resolution of a Group Specification

Convert a group specification to the data frame that
[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md)
reads. This S3 generic is exported so that extension packages, such as
`slideimp.extra`, can register methods for other classes.

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

  An object. The group specification: `slideimp` has a method for
  `data.frame`, and the optional `slideimp.extra` package adds one for a
  character platform name.

## Value

A data frame. The `data.frame` method returns `x` with duplicate rows
removed.

## Examples

``` r
df <- data.frame(feature = c("cg1", "cg2"), group = c(1, 1))
slideimp_resolve_group(df)
#>   feature group
#> 1     cg1     1
#> 2     cg2     1
```
