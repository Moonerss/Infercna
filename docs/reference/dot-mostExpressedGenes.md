# Get most expression gene

select the top n number of genes for infer cna analysis

## Usage

``` r
.mostExpressedGenes(m, ngenes)
```

## Arguments

- m:

  a matrix of genes X cells containing scRNA-seq expression data

- ngenes:

  The number or the precent of gene to select

## Value

return a vector of gene name with high expression
