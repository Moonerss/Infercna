# Compute the cell CNA score

Computes the CNA score for each cell which is composed of the CNA signal
(the sum of squares of computed CNAs) and CNA correlation (the pearson's
correlation between the cell CNA profile and the tumor's CNA profile).

## Usage

``` r
cna_compute_scores(
  cna,
  cor.method = "pearson",
  gene.quantile.for.cor = 0.5,
  gene.quantile.for.signal = 0.9,
  refCells = NULL,
  samples = NULL,
  verbose = FALSE
)
```

## Arguments

- cna:

  a matrix of gene rows by cell columns containing CNA values.

- cor.method:

  character string indicating the method to use for the pairwise
  correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'

- gene.quantile.for.cor:

  as above but for CNA correlations specifically. Default: gene.quantile

- gene.quantile.for.signal:

  as above but for CNA signal specifically. Default: gene.quantile

- refCells:

  a character vector of cell ids to exclude from average CNA profile
  that each cell is correlated to. You can pass reference normal cell
  ids to this argument if these are known. Default: NULL

- samples:

  if CNA correlations should be calculated within cell subgroups,
  provide i) a list of cell id groups, ii) a character vector of sample
  names to groups cells by, iii) TRUE to extract sample names from cell
  ids and subsequently groups. Default: NULL

- verbose:

  print progress messages. Default: FALSE

## Value

A `tibble` with three columns:

1.  CellID

2.  Signal

3.  Correlation
