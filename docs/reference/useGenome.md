# Select a Genome for infercna to use

Select your genome of choice at the start of an analysis. The available
genomes in the current implementation are, for H.sapiens, hg38 (latest)
and hg19 (preceding) and for mouse, mm10. You can see which genomes are
available via availableGenomes().

## Usage

``` r
useGenome(name)
```

## Arguments

- name:

  a character string of genome name. One of 'hg19', 'hg38', 'mm10',
  'mm39'.

## Value

genome variables are set internally. No return.
