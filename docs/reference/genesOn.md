# Retrieve Genes by their Genome Features

Retrieve genes by their genome features. e.g. all genes on chromosome 7.

## Usage

``` r
genesOn(value, attribute)
```

## Arguments

- value:

  a value or a list of values that are the genome features to filter

- attribute:

  the name(s) of the genome feature(s) in question. One or more of of
  'chr', 'arm', 'start', 'end'.

## Value

a character vector of gene names
