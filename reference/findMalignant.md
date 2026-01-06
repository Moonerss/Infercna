# Find Malignant Subset of Cells

Find the malignant and non-malignant subsets of cells from a
gene-by-cell matrix of CNA values.

## Usage

``` r
findMalignant(
  cna,
  refCells = NULL,
  samples = NULL,
  gene.quantile.for.cor = 0.5,
  gene.quantile.for.signal = 0.9,
  use.bootstraps = TRUE,
  n.bootstraps = 10000,
  prob = 0.95,
  coverage = 0.8,
  verbose = TRUE,
  ...
)
```

## Arguments

- cna:

  a matrix of gene rows by cell columns containing CNA values.

- refCells:

  a character vector or list of character vectors denoting the cell IDs
  of the known non-malignant 'Reference' cells (usually those that were
  used in the call to infercna::infercna to correct the CNA values).
  This is relevant only if these cells are in the CNA matrix. Default:
  NULL

- samples:

  a character vector of unique sample names to group the cells by. This
  is relevant only if multiple tumours are represented in the cna
  matrix, as it allows tumour-specific average CNA profiles to be
  computed in the call to Infercna::cnaCor. See Details section for more
  information. Default: `unique_sample_names(colnames(cna))`

- gene.quantile.for.cor:

  as above but for CNA correlations specifically. Default: gene.quantile

- gene.quantile.for.signal:

  as above but for CNA signal specifically. Default: gene.quantile

- use.bootstraps:

  logical; if TRUE, the function uses a bootstrapping method to
  subsample values and identify the malignant and non-malignant groups
  iteratively. This method is more sensitive to differing group sizes,
  so will be useful if you believe one group to be much smaller than the
  other. Default: TRUE

- n.bootstraps:

  number of bootstrap replicates. Relevant only if \<use.bootstraps\> is
  TRUE. Default: 10000

- prob:

  a numeric value \>= 0 and \<= 1; the minimum posterior probability
  required for an cell to be assigned to a mode. Default: 0.8

- coverage:

  the fraction of cells that must have a posterior probability higher
  than \<prob\> to one of two modes in order for the distribution to
  qualify as bimodal. Default: 0.8

- verbose:

  print progress messages. Default: TRUE

- ...:

  other arguments passed to
  [`fitBimodal`](https://moonerss.github.io/Infercna/reference/fitBimodal.md)

## Value

if bimodality was not found, returns FALSE. If bimodality was found,
returns a list containing the malignant and non-malignant cell IDs.

## Details

This function attempts to fit gaussian bimodal distributions to each of
two parameters that describe the extent of CNAs per cell. The first of
these is CNA correlation
[`cnaCor`](https://moonerss.github.io/Infercna/reference/cnaCor.md),
which measures the pearson correlations between individual cells' CNA
profiles and the average profile of their tumours of origin. The second
measure is CNA signal
[`cnaSignal`](https://moonerss.github.io/Infercna/reference/cnaSignal.md),
which computes the cell averages of squared CNA values. The function
then assigns cells residing in the 2nd (higher-value) mode by both
measures as malignant, cells residing in the 1st (low-value) mode by
both measures as non-malignant, and any remaining cells with conflicting
assignments to modes as unassigned.
