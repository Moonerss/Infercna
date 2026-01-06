# Grouped Reorder and Ungroup

Keep desired groups together but reorder members within each group. E.g.
If cells from multiple samples are being reordered, samples can be kept
together and cells reordered within. Option to also reorder the order of
groups or samples, such that groups with similar average profiles are
placed next to one another.

## Usage

``` r
grouped_reorder(
  m,
  groups,
  interorder = FALSE,
  intraorder = TRUE,
  cor.method = "pearson",
  dist.method = "euclidean",
  cluster.method = "average",
  Names = FALSE
)
```

## Arguments

- m:

  matrix to be reordered

- groups:

  groups within which reordering should take place.

- interorder:

  if TRUE, group order itself is reordered such that groups with similar
  profiles are placed near one another. E.g. Samples with similar
  average CNA profiles. Default: FALSE

- intraorder:

  if TRUE, group order itself is reordered such that groups with similar
  profiles are placed near one another. E.g. Samples with similar
  average CNA profiles. Default: TRUE

- cor.method:

  desired correlation metric. Default: 'pearson'

- dist.method:

  desired distance metric to be used on top of correlation matrix.
  Default: 'euclidean'

- cluster.method:

  desired agglomeration method. Default: 'average'

- Names:

  return the vector of ordered IDs instead of the reordered matrix.
  Default: FALSE

## Value

reordered matrix (same dimensions as input) or a character vector of
ordered column names if Names = T.
