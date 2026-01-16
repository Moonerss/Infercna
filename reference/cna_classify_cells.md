# Classify CNA scores

This function classifies the CNA scores according to the supplied signal
and correlation thresholds.

## Usage

``` r
cna_classify_cells(
  cna_scores,
  signal_threshold = 0.05,
  correlation_threshold = 0.05,
  verbose = FALSE
)
```

## Arguments

- cna_scores:

  cna scores info get from
  [`cna_compute_scores`](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)

- signal_threshold:

  The threshold of cna scores

- correlation_threshold:

  The threshold of cna correlation

- verbose:

  print progress messages. Default: FALSE

## Value

return A `tibble` with five columns:

1.  CellID

2.  Signal

3.  Correlation

4.  CNADetected

5.  Malignant

## See also

[cna_compute_scores](https://moonerss.github.io/Infercna/reference/cna_compute_scores.md)
