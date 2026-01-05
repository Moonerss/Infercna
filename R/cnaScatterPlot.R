#'
#' @title Visualise Malignant and Non-Malignant Subsets
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see `infercna::cnaSignal` and `infercna::cnaCor` for details.
#'
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param gene.quantile calculate CNA measures including only top / "hotspot" genes according to their squared CNA values across all cells. Value between 0 and 1 denoting the quantile of genes to include. Default: NULL
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to groups cells by, iii) TRUE to extract sample names from cell ids and subsequently groups. Default: NULL
#' @param verbose print progress messages. Default: FALSE
#'
#' @importFrom rlang .data
#' @import ggplot2
#'
#' @return return a list of cna signal, cna correlation, and ggplot object of plot.
#' @rdname cnaScatterPlot
#' @export
#'
cnaScatterPlot <- function(cna,
                           cor.method = 'pearson',
                           gene.quantile = NULL,
                           refCells = NULL,
                           samples = NULL,
                           verbose = FALSE
                           ) {

  if (verbose) cli::cli_alert_info('Calculate cna correlation')
  cna_cor <- cnaCor(cna, gene.quantile = gene.quantile, refCells = refCells, samples = samples)

  if (verbose) cli::cli_alert_info('Calculate cna score')
  cna_score <- cnaSignal(cna, gene.quantile = gene.quantile, refCells = refCells, samples = samples)

  plot_data <- data.frame(cna_cor = cna_cor, cna_score = cna_score)

  p <- ggplot(data = plot_data) +
    aes(x = .data[['cna_score']], y = .data[['cna_cor']]) +
    geom_point() +
    labs(x = 'CNA Signal', y = 'CNA Correlation') +
    theme_classic() +
    theme(panel.border = element_rect(colour = 'black'),
          axis.line = element_blank())

  print(p)

  return(invisible(list(cna_cor = cna_cor, cna_signal = cna_score, plot = p)))
}
