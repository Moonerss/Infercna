# Convert Relative CNA Values To Absolute

If the identities of normal cells are known, the expected CNA values of
these cells should be 0. Thus, CNA values within the range of normal
cell CNA values are corrected to 0; CNA values below the range have the
minimum substracted; CNA values above the range have the maximum
subtracted.

## Usage

``` r
refCorrect(cna, noise = NULL, isLog = FALSE, ...)
```

## Arguments

- cna:

  a matrix of rows (genes) by columns (cells) of CNA values.

- noise:

  a numeric value, which if given, increases the boundaries within which
  CNA values are considered 0. Increases by \<noise\> i.e. the bounds
  become Minimum - noise and Maximum + noise. Default: NULL

- isLog:

  boolean value indicating whether the input is in log2 form. Default:
  FALSE

- ...:

  normal cell column IDs. Expects at least two character vectors.
