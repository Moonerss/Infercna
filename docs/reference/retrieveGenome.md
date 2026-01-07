# Retrieve genome data

Returns a tibble dataframe of the current genome in use. If on function
call a character string is supplied that corresponds to a genome in the
package, that genome will instead be returned.

## Usage

``` r
retrieveGenome(name = NULL)
```

## Arguments

- name:

  a genome name. Currently one of 'hg19' (human), 'hg38' (latest human),
  'mm10' (mouse), 'mm39' (latest mouse). Default: NULL

## Value

a tibble

## See also

[`as_tibble`](https://tibble.tidyverse.org/reference/as_tibble.html)
