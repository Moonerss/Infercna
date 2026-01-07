# Get unique sample names

Extract unique sample names from cell ids

## Usage

``` r
unique_sample_names(x, sep = "-|_", pos = 1, max.nchar = NULL, replace = NULL)
```

## Arguments

- x:

  character vector of cell ids

- sep:

  separator Default: '-\|\_'

- pos:

  position of sample name given separator. First position is 1. Default:
  1

- max.nchar:

  maximum number of characters to take as sample name (starting from the
  first). Default: NULL

- replace:

  a list of character vectors to replace. Takes the form list(c(old,
  new), c(old, new)). Default: NULL

## Value

character vector of unique sample names
