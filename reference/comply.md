# Apply a Function to All Pairs of Elements in `x` or Between `x` and `y`

Apply a Function to All Pairs of Elements in `x` or Between `x` and `y`

## Usage

``` r
comply(x, y = NULL, FUN = intersect, ..., simplify = TRUE)
```

## Arguments

- x:

  a vector or a list of vectors

- y:

  a vector or a list of vectors. If y is NULL, y will be set to x, such
  that if x is a list, the function FUN will be applied to all pairs of
  vectors in the list. Default: NULL

- FUN:

  the function to be applied to every pair of elements in `x`, or to
  every pair of elements between `x` and `y`. The function should take
  two arguments and return a numeric value. Default: intersect

- ...:

  optional arguments to 'FUN'.

- simplify:

  logical; should the result be simplified to a matrix? Default: TRUE

## Value

a matrix of resulting values after applying FUN to all pairs of x and y
(or the one pair if x and y are both atomic vectors). The number of
columns and the number of rows will correspond to the lengths of x and y
respectively.
