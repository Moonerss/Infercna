## cor matrix
.hca_cor <- function(m, method) {
  method <- match.arg(method, cor_methods())
  if (method == "none") {
    return(m)
  }
  if (is_cor(m)) {
    cli::cli_alert_warning('Computing correlations over a correlation matrix...\n Set `cor.method = "none"` to skip correlation step.')
  }
  stats::cor(m, method = method)
}

## distance matrix
.hca_dist <- function(m, method, max.dist = NULL) {
  method <- match.arg(method, dist_methods())
  if (!.is_square(m)) m <- t(m)
  if (!is.null(max.dist)) {
    if (.is_square(m) && unique(diag(m)) != max.dist) {
      cli::cli_alert_warning("<max.dist> = {.val max.dist}, but <m> diagonal = {.val {unique(diag(m))}} ...")
    }
    m <- max.dist - m
  }
  if (method == "none") {
    return(stats::as.dist(m))
  }
  stats::dist(m, method = method)
}

## Hierarchical Clustering
.hca_tree <- function(d, method) {
  method <- match.arg(method, cluster_methods())
  stats::hclust(d, method = method)
}

## cut tree
.hca_cutree <- function(tree, k, h, min.size, max.size) {
  groups <- .hca_cutree_as_list(tree = tree, k = k, h = h)
  ncells <- length(tree$labels)
  min.size <- .cluster_size(min.size, ncells)
  max.size <- .cluster_size(max.size, ncells)
  lens <- lengths(groups)
  groups[lens >= min.size & lens <= max.size]
}

.cluster_size <- function(cluster.size, ncells) {
  stopifnot(cluster.size >= 0)
  stopifnot(cluster.size <= ncells)
  if (cluster.size > 1) {
    return(cluster.size)
  }
  cluster.size * ncells
}

.hca_cutree_as_list <- function(tree, h, k) {
  stopifnot(!is.null(h) | !is.null(k))
  groups <- stats::cutree(tree = tree, h = h, k = k)
  if (!has_dim(groups)) {
    return(split(names(groups), groups))
  }

  .clusterNames <- function(cutreeOutput) {
    # change df colnames (tree heights) to have 4 decimal places followed by "_"
    colnames(cutreeOutput) <- paste0(round(as.numeric(colnames(cutreeOutput)), 4), "_")
    # new clusterNames
    names(unlist(apply(cutreeOutput, 2, function(col) 1:length(table(col)))))
  }

  clusterNames <- .clusterNames(cutreeOutput = groups)
  labels <- rownames(groups)
  groups <- as.list(as.data.frame(groups))
  groups <- sapply(groups, function(ID) split(labels, ID), simplify = F)
  groups <- unlist(groups, recursive = F, use.names = F)
  groups <- stats::setNames(groups, clusterNames)
  unique(groups)
}

.hca <- function(x,
                 cor.method = cor_methods(),
                 dist.method = dist_methods(),
                 cluster.method = cluster_methods(),
                 max.dist = 1,
                 h = NULL,
                 k = NULL,
                 min.size = 5,
                 max.size = 0.5,
                 cor.end = F,
                 dist.end = F,
                 hclust.end = F) {
  res <- c()

  if (inherits(x, "matrix")) {
    x <- .hca_cor(x, method = cor.method)
    res <- c(res, list(cr = x))
    if (cor.end) {
      return(res)
    }
    x <- .hca_dist(x, method = dist.method, max.dist = max.dist)
    res <- c(res, list(dist = x))
    if (dist.end) {
      return(res)
    }
  }

  if (inherits(x, "dist")) {
    x <- .hca_tree(x, method = cluster.method)
    res <- c(res, list(tree = x, order = x$labels[x$order]))
    if (hclust.end) {
      return(res)
    }
  }

  if (inherits(x, "hclust")) {
    if (is.null(h)) h <- res$tree$height
    x <- .hca_cutree(x,
                     k = k,
                     h = h,
                     min.size = min.size,
                     max.size = max.size
    )
    res <- c(res, list(groups = x))
  }

  res
}

cor_methods <- function() {
  c("pearson", "kendall", "spearman", "none")
}

dist_methods <- function() {
  c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "none")
}

cluster_methods <- function() {
  c("average", "complete", "single", "ward.D", "ward.D2", "mcquitty", "median", "centroid")
}

# check matrix is a square
.is_square <- function(m) {
  nrow(m) == ncol(m)
}

#' @title Check cor matrix
#' @description Check the matrix whether is a cor matrix
#' @param m a matrix
#' @export
#'
is_cor <- function(m) {
  rg <- range(m)
  if ((.is_square(m)) & (rg[1] >= -1) & (rg[2] <= 1)) {
    dg <- unique(diag(m))
    return(length(dg) == 1 & dg == 1)
  }
  FALSE
}


