# Add your own Genome

Add your own Genome for infercna to use. The

## Usage

``` r
addGenome(genome, name = "userDefined")
```

## Arguments

- genome:

  genome data provided as a dataframe. The dataframe should contain
  columns 'symbol', 'start_position', 'end_position', 'chromosome_name',
  'arm'. The columns 'chromosome_name' and 'arm' should be factors, with
  the chromosome/chromosome arms ordered correctly.

- name:

  a character string; the name of your genome. Default: 'userDefined'

## Value

no return value. The genome variables will be updated internally.
