# Find Clones

Assign cells to genetic subclones from their inferred CNA profiles. You
can compute their CNA profiles using
[`Infercna`](https://moonerss.github.io/Infercna/reference/Infercna.md).

## Usage

``` r
findClones(
  cna,
  prob = 0.95,
  coverage = 0.8,
  mode.size = 10,
  clone.size = 3,
  by = "chr",
  bySampling = FALSE,
  nsamp = 2000,
  force.tries = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- cna:

  a matrix of genes X cells (variables X observations) containing CNA
  values.

- prob:

  a numeric value \>= 0 and \<= 1; the minimum posterior probability
  required for an observation to be assigned to a mode. Default: 0.95

- coverage:

  the fraction of observations that must have a posterior probability
  higher than \<prob\> to one of two modes in order for the distribution
  to qualify as bimodal. Default: 0.8

- mode.size:

  the minimum number of observations required to define a mode. Default:
  10

- clone.size:

  the minimum number of cells required to define a clone. Default: 3

- by:

  found clones by different gene position, Default: 'chr'

- bySampling, nsamp, force.tries, ...:

  see
  [`fitBimodal`](https://moonerss.github.io/Infercna/reference/fitBimodal.md)

- verbose:

  show more message
