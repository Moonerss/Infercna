# Filter Genes by their Genome Features

Filter genes to keep (the default) or to out according to genomic
features. For example, filter genes that are either on chromosome 7 or
chromosome arm 2p.

## Usage

``` r
filterGenes(x, value, attribute, out = F)
```

## Arguments

- x:

  a character vector of gene names or a matrix with gene row names.

- value:

  a value or a list of values that are the genome features to filter

- attribute:

  the name(s) of the genome feature(s) in question. One of 'chr', 'arm'.

- out:

  boolean value indicating whether the filtered genes should be thrown
  rather than kept. Default: F

## Value

filtered vector or a matrix with filtered rows, depending on class of
\<x\>.
