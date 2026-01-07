# Cell - Tumour CNA Correlations

Compute the pairwise correlations between individual cells' CNA values
and the average CNA values in their tumour of origin.

## Usage

``` r
cnaCor(
  cna,
  cor.method = "pearson",
  cell.quantile = NULL,
  gene.quantile.for.cells = NULL,
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

- cor.method:

  character string indicating the method to use for the pairwise
  correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'

- cell.quantile:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all cells. Value between 0 and 1
  denoting the quantile of genes to include. Default: NULL

- gene.quantile.for.cells:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all genes Value between 0 and 1
  denoting the quantile of genes to include, This is matched with
  \`cell.quantile\`. Default: NULL

- gene.quantile:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all cells. Value between 0 and 1
  denoting the quantile of genes to include. Default: NULL

- cell.quantile.for.genes:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all genes Value between 0 and 1
  denoting the quantile of genes to include, This is matched with
  \`gene.quantile\`. Default: NULL

- refCells:

  a character vector of cell ids to exclude from average CNA profile
  that each cell is correlated to. You can pass reference normal cell
  ids to this argument if these are known. Default: NULL

- samples:

  if CNA correlations should be calculated within cell subgroups,
  provide i) a list of cell id groups, ii) a character vector of sample
  names to group cells by, iii) TRUE to extract sample names from cell
  ids and subsequently group. Default: NULL

- ...:

  other arguments passed to unique_sample_names if samples = TRUE.

## Value

a numeric vector or list of numeric vectors
