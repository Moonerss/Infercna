#'
#' @title Compute the cell CNA score
#'
#' @description Computes the CNA score for each cell which is composed of the CNA
#' signal (the sum of squares of computed CNAs) and CNA correlation (the pearson's
#' correlation between the cell CNA profile and the tumor's CNA profile).
#'
#' @param cna a matrix of gene rows by cell columns containing CNA values.
#' @param cor.method character string indicating the method to use for the pairwise correlations. E.g. 'pearson', 'spearman'. Default: 'pearson'
#' @param gene.quantile.for.cor as above but for CNA correlations specifically. Default: gene.quantile
#' @param gene.quantile.for.signal as above but for CNA signal specifically. Default: gene.quantile
#' @param refCells a character vector of cell ids to exclude from average CNA profile that each cell is correlated to. You can pass reference normal cell ids to this argument if these are known. Default: NULL
#' @param samples if CNA correlations should be calculated within cell subgroups, provide i) a list of cell id groups, ii) a character vector of sample names to groups cells by, iii) TRUE to extract sample names from cell ids and subsequently groups. Default: NULL
#' @param verbose print progress messages. Default: FALSE
#'
#' @importFrom rlang is_character
#' @importFrom tibble tibble
#'
#' @return A \code{tibble} with three columns:
#' \enumerate{
#'   \item CellID
#'   \item Signal
#'   \item Correlation
#' }
#'
#' @export
#'
#'
cna_compute_scores <- function(cna,
                               cor.method = 'pearson',
                               gene.quantile.for.cor = 0.5,
                               gene.quantile.for.signal = 0.9,
                               refCells = NULL,
                               samples = NULL,
                               verbose = FALSE
) {

  if (verbose) cli::cli_alert_info('Calculate cna correlation')
  cna_cor <- cnaCor(cna, gene.quantile = gene.quantile.for.cor, refCells = refCells, samples = samples)

  if (verbose) cli::cli_alert_info('Calculate cna score')
  cna_score <- cnaSignal(cna, gene.quantile = gene.quantile.for.signal, refCells = refCells, samples = samples)

  res <- tibble(CellID = colnames(cna), Signal = cna_score, Correlation = cna_cor)

  return(res)

}

#'
#' @title Classify CNA scores
#'
#' @description This function classifies the CNA scores according to the supplied signal and correlation thresholds.
#'
#' @param cna_scores cna scores info get from \code{\link[Infercna]{cna_compute_scores}}
#' @param signal_threshold The threshold of cna scores
#' @param correlation_threshold The threshold of cna correlation
#' @param verbose print progress messages. Default: FALSE
#'
#' @return return A \code{tibble} with five columns:
#' \enumerate{
#'   \item CellID
#'   \item Signal
#'   \item Correlation
#'   \item CNADetected
#'   \item Malignant
#' }
#'
#' @seealso \link{cna_compute_scores}
#'
#' @export
#'
cna_classify_cells <- function(cna_scores, signal_threshold = 0.05, correlation_threshold = 0.05, verbose = FALSE) {

  cna_detected <- rep("Non-classifiable", nrow(cna_scores))

  cna_detected[which(cna_scores$Signal >= signal_threshold & cna_scores$Correlation >= correlation_threshold)] <- "Detected"
  cna_detected[which(cna_scores$Signal <  signal_threshold & cna_scores$Correlation <  correlation_threshold)] <- "Not detected"
  cna_detected[which(cna_scores$Signal <  signal_threshold & cna_scores$Correlation >= correlation_threshold)] <- "Low signal"
  cna_detected[which(cna_scores$Signal >= signal_threshold & cna_scores$Correlation <  correlation_threshold)] <- "Low correlation"

  cna_scores$CNADetected <- cna_detected

  malignant <- setNames(rep("Unresolved", nrow(cna_scores)), nm = rownames(cna_scores))
  malignant[cna_scores$CNADetected == "Detected"] <- "Malignant"
  malignant[cna_scores$CNADetected == "Not detected"] <- "Nonmalignant"

  cna_scores$Malignant <- malignant


  return (cna_scores)
}


#' @title Visualise of cna signal and correlation
#' @description Visualise Malignant and Non-Malignant Subsets of cells. This is achieved by plotting, for each cell, its CNA signal over its CNA correlation. Please see \link{cna_compute_scores}.
#'
#' @param cna_scores cna scores info get from \code{\link[Infercna]{cna_compute_scores}} or \code{\link[Infercna]{cna_classify_cells}}
#' @param signal_threshold The threshold of cna scores
#' @param correlation_threshold The threshold of cna correlation
#'
#' @return return a ggplot object
#'
#' @seealso \code{\link[Infercna]{cna_compute_scores}} \code{\link[Infercna]{cna_classify_cells}}
#'
#' @export
#'
cna_scores_plot <- function(cna_scores, signal_threshold = 0.05, correlation_threshold = 0.5) {

  if (is.element('CNADetected', colnames(cna_scores)) & is.element('Malignant', colnames(cna_scores))) {
    final_aes <- aes(x = cna_scores$Correlation, y = cna_scores$Signal, color = cna_scores$Malignant)
  } else if (is.element('Malignant', colnames(cna_scores))) {
    final_aes <- aes(x = cna_scores$Correlation, y = cna_scores$Signal, color = cna_scores$Malignant)
  } else if (is.element('CNADetected', colnames(cna_scores))) {
    final_aes <- aes(x = cna_scores$Correlation, y = cna_scores$Signal, color = cna_scores$CNADetected)
  } else {
    final_aes <- aes(x = cna_scores$Correlation, y = cna_scores$Signal)
  }

  p <- ggplot() +
    final_aes +
    geom_vline(xintercept = correlation_threshold) +
    geom_hline(yintercept = signal_threshold) +
    labs(x = 'CNA Correlation', y = 'CNA Signal') +
    theme_classic()

  return (p)
}

#'
#' @title Identify Malignant cells with cluster info
#'
#' @description Identify Malignant cells by cna signal and cna correlation with cluster info
#'
#' @param cna_scores cna scores info get from \code{\link[Infercna]{cna_compute_scores}} or \code{\link[Infercna]{cna_classify_cells}}
#' @param clusters a named vector of cell clusters with cell id, this must be supplied
#' @param signal_threshold The threshold of cna scores, this is used when \code{cna_scores} have no malignant info
#' @param correlation_threshold The threshold of cna correlation, this is used when \code{cna_scores} have no malignant info
#' @param min_cluster_cc_freq The minium cancer cell or normal cell precent in a cluster, Default 0.5
#' @param cna_matrix The copy number matrix
#' @param verbose print progress messages. Default: FALSE
#'
#' @return return A \code{tibble} with \code{Malignant} column show the cell type info
#'
#' @seealso \code{\link[Infercna]{cna_compute_scores}} \code{\link[Infercna]{cna_classify_cells}}
#'
#' @export
#'
cna_classify_cells_with_cluster <- function(cna_scores,
                                            clusters = NULL,
                                            signal_threshold = 0.05,
                                            correlation_threshold = 0.05,
                                            min_cluster_cc_freq = .5,
                                            cna_matrix = NULL,
                                            verbose = FALSE) {

  stopifnot(is.vector(clusters))

  cna_scores$Cluster <- clusters[cna_scores$CellID]

  if (!is.element('Malignant', colnames(cna_scores))) {
    if (verbose) cli::cli_alert_info("Initializing cell Classification")
    cna_scores <- cna_classify_cells(cna_scores, signal_threshold = signal_threshold, correlation_threshold = correlation_threshold)
  }

  cna_scores$CellClass <- cna_scores$Malignant

  if (verbose) cli::cli_alert_info('Classifying clusters')
  cna_scores <- .classify_clusters(cna_scores, min_cluster_cc_freq)

  if (verbose) cli::cli_alert_info("Classifying NiMC (Normal in Malignant Cluster)")
  cna_scores <- .classify_nimc(cna_scores)

  if (verbose) cli::cli_alert_info("Classifying MiNC (Malignant in Normal Cluster)")
  cna_scores <- .classify_minc(cna_scores)

  if (!is.null(cna_matrix)) {
    if (verbose) cli::cli_alert_info("Classifying intermediate scores (low signal/correlation)")
    cna_scores <- .classify_intermediate_scores(cna_scores, cna_matrix)
  }

  if (verbose) cli::cli_alert_info("Classifying malignant cells")
  cna_scores <- .classify_cells(cna_scores)

  return(cna_scores)

}


.classify_clusters <- function(data, min_cluster_cc_freq) {

  cluster_cell_class_freq <- data %>%
    dplyr::group_by(.data$Cluster, .data$CellClass) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(Freq = .data$n / sum(.data$n))

  cluster_class <- setNames(rep("Unresolved", length(unique(data$Cluster))), nm = unique(data$Cluster))

  malignant_clusters <- dplyr::filter(cluster_cell_class_freq, .data$CellClass == "Malignant", .data$Freq >= min_cluster_cc_freq) %>% dplyr::pull(.data$Cluster)
  cluster_class[malignant_clusters] <- "Malignant"

  nonmalignant_clusters <- dplyr::filter(cluster_cell_class_freq, .data$CellClass == "Nonmalignant", .data$Freq >= min_cluster_cc_freq) %>% dplyr::pull(.data$Cluster)
  cluster_class[nonmalignant_clusters] <- "Nonmalignant"

  cluster_class <- tibble(Cluster = names(cluster_class), ClusterClass = cluster_class)

  data <- data %>%
    dplyr::left_join(cluster_class, by = "Cluster")

  return (data)
}

.classify_nimc <- function(data) {

  cell_class <- setNames(data$CellClass, data$CellID)

  nrm_in_mal <- data %>%
    dplyr::filter(.data$CellClass == "Nonmalignant" & .data$ClusterClass == "Malignant")

  if (nrow(nrm_in_mal) > 0)
    cell_class[nrm_in_mal$CellID] <- "NiMC"

  data$CellClass <- cell_class

  return (data)
}

.classify_minc <- function(data) {

  cell_class <- setNames(data$CellClass, data$CellID)

  mal_in_nrm <- data %>%
    dplyr::filter(.data$CellClass == "Malignant" & .data$ClusterClass == "Nonmalignant")

  if (nrow(mal_in_nrm) > 0)
    cell_class[mal_in_nrm$CellID] <- "MiNC"

  data$CellClass <- cell_class

  return (data)
}


#' @importFrom stats cor
#'
.classify_intermediate_scores <- function(data, cna_matrix) {

  cell_class <- setNames(data$CellClass, data$CellID)

  int_scores <- data %>%
    dplyr::group_by(.data$Cluster) %>%
    dplyr::select(.data$CellID, .data$CNADetected, .data$Cluster, .data$ClusterClass, .data$CellClass) %>%
    dplyr::filter(.data$CNADetected == "Low signal" | .data$CNADetected == "Low correlation", .data$CellClass == "Unresolved")

  nonmalignant_clusters <- unique(data %>% dplyr::filter(.data$ClusterClass == "Nonmalignant") %>% dplyr::pull(.data$Cluster))
  malignant_clusters <- unique(data %>% dplyr::filter(.data$ClusterClass == "Malignant") %>% dplyr::pull(.data$Cluster))

  cna_matrix <- t(cna_matrix)

  int_scores$CNACorOwn <- rep(0, nrow(int_scores))
  int_scores$CNACorOth <- rep(0, nrow(int_scores))

  for (i in seq_len(nrow(int_scores))) {
    score_i <- int_scores[i, ]

    cna_i <- cna_matrix[, score_i$CellID]
    cna_c <- cna_matrix[, data %>% dplyr::filter(.data$Cluster == score_i$Cluster) %>% dplyr::pull(.data$CellID)]

    cna_c <- rowMeans(cna_c)

    int_scores$CNACorOwn[i] <- stats::cor(cna_i, cna_c)

    if (score_i$ClusterClass == "Malignant")
      cc <- nonmalignant_clusters
    else if (score_i$ClusterClass == "Nonmalignant")
      cc <- malignant_clusters
    else
      stop("Unclassified cluster")

    int_scores$CNACorOth[i] <- max(sapply(cc, function(x) {
      cna_x <- cna_matrix[, data %>% dplyr::filter(.data$x == data$Cluster) %>% dplyr::pull(.data$CellID)]
      cna_x <- rowMeans(cna_x)
      cor(cna_i, cna_x)
    }))
  }

  cell_class[int_scores %>% dplyr::filter(.data$CNACorOwn > 2*.data$CNACorOth, .data$ClusterClass == "Malignant") %>% dplyr::pull(.data$CellID)] <- "MbCC"
  cell_class[int_scores %>% dplyr::filter(.data$CNACorOwn > 2*.data$CNACorOth, .data$ClusterClass == "Nonmalignant") %>% dplyr::pull(.data$CellID)] <- "NbCC"
  cell_class[int_scores %>% dplyr::filter(.data$CNACorOwn <= 2*.data$CNACorOth) %>% dplyr::pull(.data$CellID)] <- "Unresolved"

  data$CellClass <- cell_class

  return (data)
}

.classify_cells <- function(data) {

  malignant <- setNames(data$CellClass, data$CellID)

  malignant[malignant == "MbCC"] <- "Malignant"
  malignant[malignant == "NbCC"] <- "Nonmalignant"
  malignant[!(malignant %in% c("Malignant", "Nonmalignant"))] <- "Unresolved"

  data$Malignant <- malignant

  return (data)
}
