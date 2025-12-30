#' @title Get most expression gene
#' @description
#' select the top n number of genes for infer cna analysis
#' @param m a matrix of genes X cells containing scRNA-seq expression data
#' @param ngenes The number or the precent of gene to select
#' @return return a vector of gene name with high expression
#'
.mostExpressedGenes <- function(m, ngenes) {
  if (ngenes < 0) cli::cli_abort("{.var ngenes} cannot be negative.")
  if (ngenes > nrow(m)) cli::cli_abort("{.var ngenes} cannot be larger than nrow(m).")

  if (ngenes >= 0 & ngenes <= 1) {
    ngenes <- ngenes * nrow(m)
  } else if (ngenes > nrow(m)) {
    ngenes <- nrow(m)
  }
  ## row mean
  rom <- rowMeans(m)
  names(sort(rowMeans(m), decreasing = T)[1:ngenes])
}


#' @title Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @description Infer Copy-Number Alterations From Single-Cell RNA-Seq Data
#' @param expr_mat a normalized expression matrix of genes X cells containing scRNA-seq expression data. It can be TPM, RPKM, CPM etc.
#' @param refCells a list of two or more character vectors, each containg cell IDs of normal cell types (one cell type per list element). Since these cell types are presumed normal cells, their CNA values can be used to correct the remaining CNA values. Note that all cell IDs should correspond to column names in expr_mat. See Infercna::refCells (normal cells) for an example. Default: NULL
#' @param window the size of the window to use for the rolling mean. Units are in number of genes. Default: 100
#' @param range values in expr_mat above and below this range will be set to values in range. See Infercna::clip() for more details. Default: c(-3, 3)
#' @param top_n a numeric value indicating if only top genes should be used in the CNA calculation and how many. if n is a fraction between 0 and 1, the number of genes included will n * nrow(expr_mat). Else the number of genes included will be n. Default: 3000
#' @param noise a numeric value, which if given, increases the boundaries within which CNA values are considered 0. Increases by <noise> i.e. the bounds become Minimum - noise and Maximum + noise. Default: 0.1
#' @param center.method method by which to center the cells after calculating CNA values. One of 'median', 'mean', etc.... Default: 'median'
#' @param verbose print progress messages. Default: FALSE
#'
#'
#' @return a matrix of genes X cells of inferred CNA values. Note that n = (window - 1)/2 genes will be lost from either extremity of the genome (ie. n genes lost at the start of chromosome 1 and at the end of chromosome Y, if the genome in question is H.sapiens.)
#' @details Correction with reference cells' <refCells> CNAs: the boundaries of reference CNA values are the boundaries for what should be considered a CNA of 0. Thus, if the boundary is -0.1 and 0.1, then a cell with CNA = -0.5 will be corrected to -0.4, a cell with CNA value of 1 will be corrected to 0.9 and a cell with CNA value of 0.05 will be corrected to 0.
#' @examples
#' \dontrun{
#' m <- mgh125
#' cna <- Infercna(expr_mat = m, refCells = refCells, top_n = 3000)
#' }
#' @rdname Infercna
#' @export
#'
Infercna <- function(expr_mat, refCells = NULL,
                     window = 100, range = c(-3, 3),
                     top_n = 3000, noise = 0.1,
                     center.method = "median",
                     verbose = FALSE) {
  if (!is.null(top_n)) {
    if (verbose) cli::cli_alert_info("Filtering the expression matrix to include only top {.val {top_n}} genes...")
    expr_mat <- expr_mat[.mostExpressedGenes(expr_mat, ngenes = top_n), ]
  }

  if (verbose) cli::cli_alert_info("Converting {.var expr_mat} to log(2) space...")
  expr_mat <- logtpm(expr_mat, bulk = F)

  if (verbose) cli::cli_alert_info("Performing mean-centering of the genes...")
  expr_mat <- rowcenter(expr_mat, by = "mean")

  if (verbose) cli::cli_alert_info("Ordering the genes by their genomic position...")
  expr_mat <- orderGenes(expr_mat)

  if (verbose) cli::cli_alert_info("Restricting expression matrix values to between {.val {range[[1]]}}  and {.val {range[[2]]}} ...")
  expr_mat <- clip(expr_mat, range = range)

  if (verbose) cli::cli_alert_info("Converting {.var expr_mat} from log(2) space...")
  expr_mat <- unlogtpm(expr_mat, bulk = F)

  if (verbose) cli::cli_alert_info("Preparing to calculate CNA values on each chromosome in turn...")
  ms <- splitGenes(expr_mat, by = "chr")

  if (verbose) cli::cli_alert_info("Calculating rolling means with a window size of {.val {window}} genes...")
  cna <- sapply(ms, function(expr_mat) try(runMean(expr_mat, k = window, verbose = verbose)), simplify = F)
  cna <- cna[sapply(cna, class) != "try-error" | !sapply(class, isFALSE)]
  cna <- Reduce(rbind, cna)
  if (verbose) cli::cli_alert_info("Converting CNA values to log(2) space...")
  cna <- logtpm(cna, bulk = F)
  if (verbose) cli::cli_alert_info("Performing {.val {center.method}}-centering of the cells...")
  cna <- colcenter(cna, by = center.method)

  if (!is.null(refCells)) {
    if (verbose) cli::cli_alert_info("Correcting CNA profiles using CNA values from {.var refCells}...")
    Args <- c(list(cna = cna, noise = noise, isLog = TRUE), refCells)
    cna <- do.call(refCorrect, Args)
  }
  if (verbose) cli::cli_alert_success("Done!")
  return(cna)
}
