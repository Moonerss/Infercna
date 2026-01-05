#' @title Find the Cells with the highest CNA signal
#' @description Find Cells in the top nth quantile for CNA Signal values
#' @param cna a matrix of cell rows by cell columns containing CNA values
#' @param cell.quantile calculate CNA measures including only top / "hotspot" cells according to their squared CNA values across all genes. Value between 0 and 1 denoting the quantile of cells to include.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include.
#' @importFrom stats quantile
#' @return cell names in the top nth quantile, where n is specified via <cell.quantile>
#' @rdname cnaHotspotCells
#' @export
#'
cnaHotspotCells <- function(cna, cell.quantile = NULL, gene.quantile = NULL) {
  cna <- as.matrix(cna)
  stopifnot(is.numeric(cell.quantile))
  if (!is.null(gene.quantile)) {
    stopifnot(is.numeric(gene.quantile))
    cna <- cna[cnaHotspotGenes(cna, gene.quantile = gene.quantile), ]
  }
  msq <- colMeans(cna^2)
  names(msq)[msq >= quantile(msq, cell.quantile)]
}

#' @title Find the Genes with the highest CNA signal
#' @description Find Genes in the top nth quantile for CNA Signal values
#' @param cna a matrix of gene rows by cell columns containing CNA values
#' @param cell.quantile calculate CNA measures including only top / "hotspot" cells according to their squared CNA values across all genes. Value between 0 and 1 denoting the quantile of cells to include.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include.
#' @importFrom stats quantile
#' @return gene names in the top nth quantile, where n is specified via <gene.quantile>
#' @rdname cnaHotspotGenes
#' @export
#'
cnaHotspotGenes <- function(cna, gene.quantile = NULL, cell.quantile = NULL) {
  cna <- as.matrix(cna)
  stopifnot(is.numeric(gene.quantile))
  if (!is.null(cell.quantile)) {
    stopifnot(is.numeric(cell.quantile))
    cna <- cna[, cnaHotspotCells(cna, cell.quantile = cell.quantile)]
  }
  msq <- rowMeans(cna^2)
  names(msq)[msq >= quantile(msq, gene.quantile)]
}


#'
#' @title Calculate the Means of Squared CNA Values
#' @description Calculate the Mean of Squared CNA Values
#'
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param cell.quantile.for.genes calculate CNA measures including only top / "hotspot" cells according to their squared CNA values across all genes. Value between 0 and 1 denoting the quantile of cells to include. Default: NULL
#' @param refCells a character vector of cell ids (e.g. normal reference cell ids) to exclude from calculation of cnaHotspotGenes. Only relevant if gene.quantile is not NULL. Default: NULL
#' @param samples if cnaHotspotGenes should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to group cells by, iii) TRUE to extract sample names from cell ids and subsequently group. Default: NULL
#' @param ... other arguments passed to unique_sample_names if samples = TRUE.
#'
#' @importFrom rlang is_character
#'
#' @return a numeric vector of CNA signal values or the Mean of Squared CNA values
#' @rdname cnaSignal
#' @export
#'
cnaSignal <- function(cna,
                      gene.quantile = NULL,
                      cell.quantile.for.genes = NULL,
                      refCells = NULL,
                      samples = NULL,
                      ...) {
  cna <- as.matrix(cna)

  .cnaSignal <- function(cna) {
    sqmat <- cna^2
    colMeans(sqmat)
  }

  if (is_null(gene.quantile)) {
    return(.cnaSignal(cna))
  }

  if (is_null(refCells)) {
    tmp <- cna
  } else {
    tmp <- cna[, !colnames(cna) %in% unlist(refCells)]
  }

  if (is_null(samples)) {
    genes <- cnaHotspotGenes(tmp, gene.quantile = gene.quantile, cell.quantile = cell.quantile.for.genes)
    return(.cnaSignal(cna[genes, ]))
  }

  if (is_true(samples)) {
    samples <- unique_sample_names(colnames(cna), ...)
    cli::cli_alert_info(str_c('Samples identified:\n', paste0(samples, collapse = "\n")))
  }

  if (is_character(samples)) {
    samples <- subsplit(colnames(cna), pattern = samples)
    samples <- samples[lengths(samples) != 0]
  }

  stopifnot(is_list(samples))

  if (!is_null(gene.quantile)) {
    tmplist <- sapply(samples, function(cells) tmp[, colnames(tmp) %in% cells, drop = FALSE], simplify = F)
    genelist <- lapply(tmplist, cnaHotspotGenes, gene.quantile = gene.quantile, cell.quantile = cell.quantile.for.genes)
    rm(tmplist)
  } else {
    genelist <- replicate(length(samples), rownames(cna), simplify = F)
  }

  cnalist <- Map(function(m, x, y) m[x, y],
                 x = genelist,
                 y = samples,
                 MoreArgs = list(m = cna)
  )

  Unlist(sapply(cnalist, .cnaSignal, simplify = F), nested.names = T)

}

