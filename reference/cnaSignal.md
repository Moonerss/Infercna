# Calculate the Means of Squared CNA Values

Calculate the Mean of Squared CNA Values

## Usage

``` r
cnaSignal(
  cna,
  gene.quantile = NULL,
  cell.quantile.for.genes = NULL,
  refCells = NULL,
  samples = NULL,
  ...
)
```

## Arguments

- cna:

  a matrix of gene rows by cell columns containing CNA values.

- gene.quantile:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all cells. Value between 0 and 1
  denoting the quantile of genes to include. Default: NULL

- cell.quantile.for.genes:

  calculate CNA measures including only top / "hotspot" cells according
  to their squared CNA values across all genes. Value between 0 and 1
  denoting the quantile of cells to include. Default: NULL

- refCells:

  a character vector of cell ids (e.g. normal reference cell ids) to
  exclude from calculation of cnaHotspotGenes. Only relevant if
  gene.quantile is not NULL. Default: NULL

- samples:

  if cnaHotspotGenes should be calculated within cell subgroups,
  provide i) a list of cell id groups, ii) a character vector of sample
  names to group cells by, iii) TRUE to extract sample names from cell
  ids and subsequently group. Default: NULL

- ...:

  other arguments passed to unique_sample_names if samples = TRUE.

## Value

a numeric vector of CNA signal values or the Mean of Squared CNA values
