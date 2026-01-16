# Identify Malignant cells with cluster info

Identify Malignant cells by cna signal and cna correlation with cluster
info

## Usage

``` r
cna_classify_cells_with_cluster(
  cna_scores,
  clusters = NULL,
  signal_threshold = 0.05,
  correlation_threshold = 0.05,
  min_cluster_cc_freq = 0.5,
  cna_matrix = NULL,
  verbose = FALSE
)
```

## Arguments

- cna_scores:

  cna scores info get from
  [`cna_compute_scores`](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)
  or
  [`cna_classify_cells`](https://moonerss.github.io/Infercna/reference/cna_classify_cells.md)

- clusters:

  a named vector of cell clusters with cell id, this must be supplied

- signal_threshold:

  The threshold of cna scores, this is used when `cna_scores` have no
  malignant info

- correlation_threshold:

  The threshold of cna correlation, this is used when `cna_scores` have
  no malignant info

- min_cluster_cc_freq:

  The minium cancer cell or normal cell precent in a cluster, Default
  0.5

- cna_matrix:

  The copy number matrix

- verbose:

  print progress messages. Default: FALSE

## Value

return A `tibble` with `Malignant` column show the cell type info

## See also

[`cna_compute_scores`](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)
[`cna_classify_cells`](https://moonerss.github.io/Infercna/reference/cna_classify_cells.md)
