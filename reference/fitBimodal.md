# Fit a Bimodal Gaussian Distribution

Fit a bimodal gaussian distribution to a set of observations.

## Usage

``` r
fitBimodal(
  x,
  prob = 0.95,
  coverage = 0.8,
  size = 10,
  assign = FALSE,
  boolean = FALSE,
  verbose = TRUE,
  maxit = 5000,
  maxrestarts = 100,
  bySampling = FALSE,
  nsamp = 2000,
  ...
)
```

## Arguments

- x:

  a named numeric vector of cells/observations or a matrix of genes X
  cells (variables X observations). If the latter, the column means are
  first computed.

- prob:

  a numeric value \>= 0 and \<= 1; the minimum posterior probability
  required for an observation to be assigned to a mode. Default: 0.95

- coverage:

  the fraction of observations that must have a posterior probability
  higher than \<prob\> to one of two modes in order for the distribution
  to qualify as bimodal. Default: 0.8

- size:

  the minimum number of observations that must be assigned to a mode in
  order for the distribution to qualify as bimodal. Default: 10

- assign:

  if set to TRUE, returns a list of length two containing the vector
  names that were assigned to each mode. Default: FALSE

- boolean:

  if set to TRUE, returns a boolean value indicating whether the
  distribution is bimodal. Default: FALSE

- verbose:

  print progress messages. Default: TRUE

- maxit:

  the maximum number of iterations. Default: 5000

- maxrestarts:

  the maximum number of restarts allowed. See `normalmixEM` for details.
  Default: 100

- bySampling:

  logical; if TRUE, the function uses a bootstrapping method to
  subsample values and identify the two modes iteratively. This method
  is more sensitive to differing mode sizes, so will be useful if you
  believe one group to be much smaller than the other. Default: TRUE

- nsamp:

  the number of bootstrap replicates.

- ...:

  Other argument for methods

## Value

The posterior probabilities of each observation to one of two modes. If
boolean = TRUE, return a boolean value indicating whether bimodality was
found. If assign = TRUE, return a list of length two with the
observations (IDs) in each mode.

## See also

[`norMixEM`](https://rdrr.io/pkg/nor1mix/man/norMixFit.html)
