# Rolling Means

Apply a rolling window mean to a matrix or vector.

## Usage

``` r
runMean(m, k = 100, endrule = "mean", align = "center", verbose = T)
```

## Arguments

- m:

  a numeric vector or matrix. If the latter, each column will be
  processed separately.

- k:

  width of rolling window. Default: 100

- endrule:

  character string indicating how the values at the beginning and the
  end of the data should be treated. One of "mean", "trim", "keep",
  "constant". See caTools::runmean for more details. Default: 'mean'

- align:

  specifies whether result should be centered (default), left-aligned or
  right-aligned. See caTools::runmean for more details. Default:
  'center'

- verbose:

  print progress messages. Default: TRUE

## Value

a numeric vector or matrix of the same size as \<m\>. Only in case of
endrule=trim, the output vectors will be shorter and output matrices
will have fewer rows.

## See also

[`runmean`](https://rdrr.io/pkg/caTools/man/runmean.html)
