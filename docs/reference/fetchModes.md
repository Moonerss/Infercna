# Find small clone

Find all small clones by gene position

## Usage

``` r
fetchModes(
  cna,
  prob = 0.95,
  coverage = 0.8,
  size = 10,
  by = "chr",
  minGenes = 50,
  bySampling = FALSE,
  nsamp = 2000,
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

- size:

  the minimum number of observations required to define a mode. Default:
  10

- by:

  found clones by different gene position, Default: `chr`

- minGenes:

  The min number of genes to found clones.

- bySampling, nsamp, ...:

  see
  [`fitBimodal`](https://moonerss.github.io/Infercna/reference/fitBimodal.md)

## Value

return a list of all small clone
