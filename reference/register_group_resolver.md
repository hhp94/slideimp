# Register a Group Resolver

Called by the companion package `slideimp.extra` from their `.onLoad()`
to register a callback that turns character manifest into a `data.frame`
for
[`group_imp()`](https://hhp94.github.io/slideimp/reference/group_imp.md)/[`prep_groups()`](https://hhp94.github.io/slideimp/reference/prep_groups.md).

## Usage

``` r
register_group_resolver(resolver)
```

## Arguments

- resolver:

  A function that takes the manifest character and returns a
  `data.frame`.

## Value

Invisibly returns `NULL`. The function is called only for its side
effect of registering the resolver.

## Examples

``` r
# Typically called from `.onLoad` by `slideimp.extra` package, not by users directly.
dummy_resolver <- function(x) data.frame(feature = character(), group = character())
register_group_resolver(dummy_resolver)
```
