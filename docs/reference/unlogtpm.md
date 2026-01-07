# Convert TPM to logTPM

Convert TPM to logTPM, i.e. using log2(TPM/10 + 1).

Convert logTPM to TPM, i.e. using 10\*(2^(TPM)-1).

## Usage

``` r
logtpm(m, bulk = F)

unlogtpm(m, bulk = F)
```

## Arguments

- m:

  matrix of logTPM values (gene rows; cell columns)

- bulk:

  if bulk then instead uses 2^(TPM)-1. i.e. no scaling. Default: F

## Value

TPM matrix

TPM matrix

## Details

TPM/10 is used for single cells since 100,000 is a more reasonable
estimate than 1,000,000 for the number of RNA transcripts in a cell.
1,000,000 is reasonable estimate for bulk samples that contain multiple
cells.

TPM/10 is used for single cells since 100,000 is a more reasonable
estimate than 1,000,000 for the number of RNA transcripts in a cell.
1,000,000 is reasonable estimate for bulk samples that contain multiple
cells.
