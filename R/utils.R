#' @title check dim
#' @param x the object to check
#' check the object have dim attribute
has_dim <- function(x) {
  if (is.data.frame(x)) x <- as.matrix(x)
  !is.null(attr(x, "dim"))
}

#' @title Convert TPM to logTPM
#' @description Convert TPM to logTPM, i.e. using log2(TPM/10 + 1).
#' @param m matrix of logTPM values (gene rows; cell columns)
#' @param bulk if bulk then instead uses log2(TPM + 1). i.e. no scaling. Default: F
#' @return TPM matrix
#' @details TPM/10 is used for single cells since 100,000 is a more reasonable estimate than 1,000,000 for the number of RNA transcripts in a cell. 1,000,000 is reasonable estimate for bulk samples that contain multiple cells.
#' @rdname unlogtpm
#' @export
logtpm <- function(m, bulk = F) {
  if (has_dim(m)) {
    m <- as.matrix(m)
  }
  if (bulk) {
    x <- 1
  } else {
    x <- 10
  }
  log2((m / x) + 1)
}

#' @title Convert logTPM to TPM
#' @description Convert logTPM to TPM, i.e. using 10*(2^(TPM)-1).
#' @param m matrix of logTPM values (gene rows; cell columns)
#' @param bulk if bulk then instead uses 2^(TPM)-1. i.e. no scaling. Default: F
#' @return TPM matrix
#' @details TPM/10 is used for single cells since 100,000 is a more reasonable estimate than 1,000,000 for the number of RNA transcripts in a cell. 1,000,000 is reasonable estimate for bulk samples that contain multiple cells.
#' @rdname unlogtpm
#' @export
unlogtpm <- function(m, bulk = F) {
  if (has_dim(m)) {
    m <- as.matrix(m)
  }
  if (bulk) {
    x <- 1
  } else {
    x <- 10
  }
  # (2^m) * x - 1
  (2^(m) - 1) * x
}

#' @title Center a matrix row-wise
#' @description Center a matrix row-wise
#' @param m a matrix or Matrix
#' @param by either "mean", "median" or a numeric vector of length equal to the number of rows of ‘m’. Default: "mean"
#' @return row-centered matrix
#' @rdname rowcenter
#' @export
rowcenter <- function(m, by = "mean") {
  m <- as.matrix(m)
  if (by == "mean") {
    by <- rowMeans(m, na.rm = T)
  } else if (by == "median") {
    by <- matrixStats::rowMedians(m, na.rm = T)
  } else {
    stopifnot(is.numeric(by) & length(by) == nrow(m))
  }

  t(scale(t(m), center = by, scale = F))
}

#' @title Center a matrix column-wise
#' @description Center a matrix column-wise
#' @param m a matrix or Matrix
#' @param by either "mean", "median" or a numeric vector of length equal to the number of columns of ‘m’. Default: "mean"
#' @return column-centered matrix
#' @rdname colcenter
#' @export
colcenter <- function(m, by = "mean") {
  m <- as.matrix(m)
  if (by == "mean") {
    by <- colMeans(m, na.rm = T)
  } else if (by == "median") {
    by <- matrixStats::colMedians(m, na.rm = T)
  } else {
    stopifnot(is.numeric(by) & length(by) == ncol(m))
  }
  scale(m, center = by, scale = F)
}

#' @title Squish matrix values into range
#' @description Squish matrix values into range
#' @param m matrix to manipulate
#' @param range numeric vector of length two giving desired output range. Default: c(-3, 3)
#' @rdname clip
#' @export
clip <- function(m, range = c(-3, 3)) {
  m <- as.matrix(m)
  m[m < range[[1]]] <- range[[1]]
  m[m > range[[2]]] <- range[[2]]
  m
}

#' @title Unlist, keeping original list or element names
#' @description Unlist, keeping original list or element names
#' @param L list to flatten
#' @param nested.names logical; keep nested list names rather than expanding list names. Default: FALSE
#' @return a vector of the same length as the combined lengths of list elements. Names will either be the list names replicated, or, if nested.names is TRUE, the original list element names will be kept.
#' @seealso
#'  \code{\link[stats]{setNames}}
#' @rdname Unlist
#' @export
#' @importFrom stats setNames
Unlist <- function(L, nested.names = FALSE) {
  if (nested.names) {
    Names <- unlist(sapply(L, names), use.names = F)
  } else {
    Names <- rep(names(L), lengths(L))
  }
  stats::setNames(unlist(L), Names)
}


#' @title Flip the (nested) elements of a character vector (or list) with its names
#' @description A convenience function to switch the names and elements of a named character vector or a named list of character vectors.
#' @param X A named character vector or a named list of character vectors.
#' @return A named character vector or a named list of character vectors. If the former, will be of the same length as <X>. If the latter, will be of the same length as there are unique elements across all vectors in the list.
#' @seealso
#'  \code{\link[stats]{setNames}}
#' @rdname flip
#' @export
#' @importFrom stats setNames
flip <- function(X) {
  stopifnot(!is.null(names(X)))
  if (is.character(X)) {
    return(stats::setNames(names(X), X))
  } else if (is.list(X) & !has_dim(X)) {
    return(split(rep(names(X), lengths(X)), unlist(X, use.names = F)))
  } else {
    stop("X should be a list or a character vector.")
  }
}


#' @title Rolling Means
#' @description Apply a rolling window mean to a matrix or vector.
#' @param m a numeric vector or matrix. If the latter, each column will be processed separately.
#' @param k width of rolling window. Default: 100
#' @param endrule character string indicating how the values at the beginning and the end of the data should be treated. One of "mean", "trim", "keep", "constant". See caTools::runmean for more details. Default: 'mean'
#' @param align specifies whether result should be centered (default), left-aligned or right-aligned. See caTools::runmean for more details. Default: 'center'
#' @param verbose print progress messages. Default: TRUE
#' @return a numeric vector or matrix of the same size as <m>. Only in case of endrule=trim, the output vectors will be shorter and output matrices will have fewer rows.
#' @seealso
#'  \code{\link[caTools]{runmean}}
#' @rdname runMean
#' @export
#' @importFrom caTools runmean
runMean <- function(m, k = 100, endrule = "mean", align = "center", verbose = T) {
  if (!is.null(dim(m))) {
    m <- as.matrix(m)
  }

  if (is.null(dim(m))) {
    return(FALSE)
  }
  if (nrow(m) == 0) {
    return(FALSE)
  }

  if (nrow(m) < k) {
    k <- nrow(m)
    if (verbose) cli::cli_alert_info("Adjusting window to the max. number of genes in chromosome ( {.val {k}} )")
  }

  mout <- caTools::runmean(m, k = k, endrule = endrule, align = align)

  if (!is.null(dim(m))) {
    colnames(mout) <- colnames(m)
    rownames(mout) <- rownames(m)
    mout <- as.data.frame(mout)
  } else {
    names(mout) <- names(m)
  }

  mout
}
