# Flip the (nested) elements of a character vector (or list) with its names

A convenience function to switch the names and elements of a named
character vector or a named list of character vectors.

## Usage

``` r
flip(X)
```

## Arguments

- X:

  A named character vector or a named list of character vectors.

## Value

A named character vector or a named list of character vectors. If the
former, will be of the same length as \<X\>. If the latter, will be of
the same length as there are unique elements across all vectors in the
list.

## See also

[`setNames`](https://rdrr.io/r/stats/setNames.html)
