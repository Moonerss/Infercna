# Find the Genes with the highest CNA signal

Find Genes in the top nth quantile for CNA Signal values

## Usage

``` r
cnaHotspotGenes(cna, gene.quantile = NULL, cell.quantile = NULL)
```

## Arguments

- cna:

  a matrix of gene rows by cell columns containing CNA values

- gene.quantile:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all cells. Value between 0 and 1
  denoting the quantile of genes to include.

- cell.quantile:

  calculate CNA measures including only top / "hotspot" cells according
  to their squared CNA values across all genes. Value between 0 and 1
  denoting the quantile of cells to include.

## Value

gene names in the top nth quantile, where n is specified via
\<gene.quantile\>
