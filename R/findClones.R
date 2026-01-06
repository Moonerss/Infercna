#'
#' @title Find small clone
#' @description
#' Find all small clones by gene position
#'
#' @param cna a matrix of genes X cells (variables X observations) containing CNA values.
#' @param prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param size the minimum number of observations required to define a mode. Default: 10
#' @param by found clones by different gene position, Default: \code{chr}
#' @param minGenes The min number of genes to found clones.
#' @param bySampling,nsamp,... see \code{\link[Infercna]{fitBimodal}}
#'
#' @return return a list of all small clone
#' @export
#'
fetchModes <- function(cna,
                       prob = 0.95,
                       coverage = 0.8,
                       size = 10,
                       by = 'chr',
                       minGenes = 50,
                       bySampling = FALSE,
                       nsamp = 2000,
                       ...) {

  mats <- splitGenes(cna, by = by)
  Rows <- unlist(sapply(mats, nrow))
  Rows[is.null(Rows)] <- 0
  mats <- mats[Rows >= minGenes]
  modes <- sapply(mats, fitBimodal, prob = prob, coverage = coverage, size = size, assign = T, bySampling = bySampling, nsamp = nsamp, ...)
  modes[!sapply(modes, is_false)]
}


#'
#' @title Find Clones
#' @description
#' Assign cells to genetic subclones from their inferred CNA profiles. You can compute their CNA profiles using \code{\link[Infercna]{Infercna}}.
#'
#' @param cna a matrix of genes X cells (variables X observations) containing CNA values.
#' @param prob a numeric value >= 0 and <= 1; the minimum posterior probability required for an observation to be assigned to a mode. Default: 0.95
#' @param coverage the fraction of observations that must have a posterior probability higher than <prob> to one of two modes in order for the distribution to qualify as bimodal. Default: 0.8
#' @param mode.size the minimum number of observations required to define a mode. Default: 10
#' @param clone.size the minimum number of cells required to define a clone. Default: 3
#' @param by found clones by different gene position, Default: 'chr'
#' @param bySampling,nsamp,force.tries,... see \code{\link[Infercna]{fitBimodal}}
#' @param verbose show more message
#'
#' @export
#'
findClones <- function(cna,
                       prob = 0.95,
                       coverage = 0.8,
                       mode.size = 10,
                       clone.size = 3,
                       by = 'chr',
                       bySampling = FALSE,
                       nsamp = 2000,
                       force.tries = FALSE,
                       verbose = FALSE,
                       ...) {

  if (is_list(as.matrix(cna))) {
    if (verbose) {
      modes <- sapply(cna, fitBimodal, prob = prob, coverage = coverage, size = mode.size, assign = T, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...)
    } else {
      modes <- suppressWarnings(sapply(cna, fitBimodal, prob = prob, coverage = coverage, size = mode.size, assign = T, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...))
    }
  } else {
    if (verbose) {
      modes <- fetchModes(cna, prob = prob, coverage = coverage, size = mode.size, by = by, bySampling = bySampling, nsamp = nsamp,force.tries = force.tries, ...)
    } else {
      modes <- suppressWarnings(fetchModes(cna, prob = prob, coverage = coverage, size = mode.size, by = by, bySampling = bySampling, nsamp = nsamp, force.tries = force.tries,...))
    }
  }

  expandToClones(modes, greaterThan = clone.size)
}


#'
#' @title expand clones
#' @description
#' Order the clones
#'
#' @param modes The result of \code{\link[Infercna]{fetchModes}}
#' @param greaterThan The clone number, Default NULL
#'
#' @export
#'
expandToClones <- function(modes, greaterThan = NULL) {

  modes <- sapply(modes, function(x) {
    stats::setNames(unlist(x), rep(names(x), lengths(x)))},
    simplify = F)

  modes <- unlist(modes)
  clones <- sapply(split(names(modes), modes), paste0, collapse = '--')
  clones <- split(names(clones), clones)
  if (is.null(greaterThan)) return(clones)
  clones[lengths(clones) > greaterThan]
}
