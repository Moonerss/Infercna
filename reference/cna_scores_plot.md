# Visualise of cna signal and correlation

Visualise Malignant and Non-Malignant Subsets of cells. This is achieved
by plotting, for each cell, its CNA signal over its CNA correlation.
Please see
[cna_compute_scores](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md).

## Usage

``` r
cna_scores_plot(
  cna_scores,
  signal_threshold = 0.05,
  correlation_threshold = 0.5
)
```

## Arguments

- cna_scores:

  cna scores info get from
  [`cna_compute_scores`](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)
  or
  [`cna_classify_cells`](https://moonerss.github.io/Infercna/reference/cna_classify_cells.md)

- signal_threshold:

  The threshold of cna scores

- correlation_threshold:

  The threshold of cna correlation

## Value

return a ggplot object

## See also

[`cna_compute_scores`](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)
[`cna_classify_cells`](https://moonerss.github.io/Infercna/reference/cna_classify_cells.md)
