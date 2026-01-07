# Find the Cells with the highest CNA signal

Find Cells in the top nth quantile for CNA Signal values

## Usage

``` r
cnaHotspotCells(cna, cell.quantile = NULL, gene.quantile = NULL)
```

## Arguments

- cna:

  a matrix of cell rows by cell columns containing CNA values

- cell.quantile:

  calculate CNA measures including only top / "hotspot" cells according
  to their squared CNA values across all genes. Value between 0 and 1
  denoting the quantile of cells to include.

- gene.quantile:

  calculate CNA measures including only top / "hotspot" genes according
  to their squared CNA values across all cells. Value between 0 and 1
  denoting the quantile of genes to include.

## Value

cell names in the top nth quantile, where n is specified via
\<cell.quantile\>
