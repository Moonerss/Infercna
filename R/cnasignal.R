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


#' @title Calculate the Means of Squared CNA Values
#' @description Calculate the Mean of Squared CNA Values
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param refCells a character vector of cell ids (e.g. normal reference cell ids) to exclude from calculation of cnaHotspotGenes. Only relevant if gene.quantile is not NULL. Default: NULL
#' @return a numeric vector of CNA signal values or the Mean of Squared CNA values
#' @rdname cnaSignal
#' @export
#'
cnaSignal <- function(cna,
                      gene.quantile = NULL,
                      cell.quantile.for.genes = NULL,
                      refCells = NULL,
                      samples = NULL) {
  cna <- as.matrix(cna)

  .cnaSignal <- function(cna) {
    sqmat <- cna^2
    colMeans(sqmat)
  }

  if (is.null(gene.quantile)) {
    return(.cnaSignal(cna))
  }

  if (is.null(refCells)) {
    tmp <- cna
  } else {
    tmp <- cna[, !colnames(cna) %in% unlist(refCells)]
  }

  genes <- cnaHotspotGenes(tmp, gene.quantile = gene.quantile, cell.quantile = cell.quantile.for.genes)
  return(.cnaSignal(cna[genes, ]))

}

