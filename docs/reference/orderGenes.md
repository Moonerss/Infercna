# Order Genes by their Genomic Positions

Order Genes by their Genomic Positions

## Usage

``` r
orderGenes(x)
```

## Arguments

- x:

  a matrix with gene names as rows or a character vector of gene names.

## Value

ordered matrix or character vector

## Examples

``` r
if (FALSE) { # \dontrun{
m <- Infercna::useData()
orderGenes(m) %>%
  rownames() %>%
  head()
orderGenes(rownames(m)) %>% head()
} # }
```
