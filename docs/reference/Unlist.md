# Unlist, keeping original list or element names

Unlist, keeping original list or element names

## Usage

``` r
Unlist(L, nested.names = FALSE)
```

## Arguments

- L:

  list to flatten

- nested.names:

  logical; keep nested list names rather than expanding list names.
  Default: FALSE

## Value

a vector of the same length as the combined lengths of list elements.
Names will either be the list names replicated, or, if nested.names is
TRUE, the original list element names will be kept.

## See also

[`setNames`](https://rdrr.io/r/stats/setNames.html)
