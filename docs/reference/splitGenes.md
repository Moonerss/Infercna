# Split Genes By Chromosome (Arm)

Split Genes By Chromosome (Arm). Input can be one of a matrix to be
split into several matrices or a character vector of gene names to be
split into several character vectors.

## Usage

``` r
splitGenes(x, by = "chr")
```

## Arguments

- x:

  a matrix with gene names as rows or a character vector of gene names.

- by:

  a string; one of 'chr' or 'arm', determining how the matrix or
  character vector should be split. Default: 'chr'

## Value

if a matrix was provided, a list of matrices. If a character vector was
provided, a list of character vectors. Both lists will be of length \#
of chromosome (arms) in the genome in question. Ie. length 24 for
H.sapiens.

## Examples

``` r
if (FALSE) { # \dontrun{
m <- Infercna::useData()
a <- splitGenes(x = m)
b <- splitGenes(x = rownames(m))
all(sapply(a, nrow) == lengths(b))
names(a)
names(splitGenes(x = m, by = "arm"))
} # }
```
