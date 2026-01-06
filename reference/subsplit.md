# Split a vector, list, matrix or dataframe into into a list of matrices or nested list according to substring matches.

subsplit() divides the data in 'x' into groups defined by the presence
of shared substrings. If 'x' is a list, the substrings are derived from
'names(x)'. If 'x' is a matrix or data frame, the substrings are instead
derived from 'colnames(x)' or, if 'by.row=FALSE', from 'rownames(x)'.

## Usage

``` r
subsplit(
  x,
  by.row = FALSE,
  sep = "\\.|-|_",
  pattern = NULL,
  pos = 1,
  max.nchar = NULL,
  replace = NULL,
  na.rm = TRUE
)
```

## Arguments

- x:

  A character vector, a named list, a matrix or data frame.

- by.row:

  If 'TRUE' and 'x' is a matrix or data frame, the substrings are
  derived from 'rownames(x)' instead of 'colnames(x)'. Default: FALSE

- sep:

  The pattern to split each element of 'x' by. Default: '\\\|-\|\_'

- pattern:

  The substrings to split by. If 'NULL', substrings will instead be
  defined from the delimiter 'sep' and position index 'pos' with an
  internal call to scalop::substri(). Default: NULL

- pos:

  The substring index(es) to be returned. Default: 1

- max.nchar:

  'NULL', or a maximum number of characters that the substrings should
  be trimmed to contain. Default: NULL

- replace:

  A character vector containing two elements. The first is the pattern
  to be replaced and the second is the pattern that it will be replaced
  with. Alternatively, a list of such character vectors, or 'NULL' if no
  replacement is desired. Default: NULL

- na.rm:

  Remove NA substring positions. Default: TRUE

## Value

A list of the same type as 'x' and of length equal to the number of
unique substrings found.

## See also

[`str_subset`](https://stringr.tidyverse.org/reference/str_subset.html)

## Examples

``` r
if (FALSE) { # \dontrun{
Names = c(rep('dog.cat', 5), rep('dog', 2), rep('cat', 2), 'frog')
x.vec = replicate(expr = sample(100, 5), n = 10, simplify = F)
names(x.vec) = Names
subsplit(x.vec, pattern = c('dog', 'cat', 'frog'))
x.mat = replicate(expr = sample(100, 5), n = 10, simplify = T);
colnames(x.mat) = Names
subsplit(x.mat)
subsplit(x.mat, by.row = T)
} # }
```
