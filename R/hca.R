#' @title Hierarchical clustering analysis
#' @description Hierarchical clustering analysis
#' @param x a matrix (features X observations), an object of class `dist` or an object of class `hclust`.
#' @param cor.method option if <x> is a matrix. One of 'pearson' (default), 'kendall', 'spearman' or 'none' (no correlation coefficents computed). Default: cor_methods()
#' @param dist.method option if <x> is a matrix or correlation matrix. One of 'euclidean' (default), 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski' or 'none' (in which case stats::as.dist() will be used). Default: dist_methods()
#' @param max.dist maximum distance between observations. Default: 1
#' @param cluster.method one of "average" (default), "complete", "single", "ward.D", "ward.D2", "mcquitty", "median" or "centroid". Default: cluster_methods()
#' @param h height(s) at which to cut the tree. If NULL, h will be set to all tree heights. Default: NULL
#' @param k number of groups to return from tree. If <h> and <k> are both not NULL, <k> takes precedence. Default: NULL
#' @param min.size minimum allowed cluster/group size. Values between 0 and 1 are interpreted as fractions of total count. Groups smaller than <min.size> are filtered out. Default: 5
#' @param max.size maximum allowed cluster/group size. Values between 0 and 1 are interpreted as fractions of total count. Groups larger than <max.size> are filtered.out. Default: 0.5
#' @return list
#' @rdname hca
#' @export
hca <- function(x,
                cor.method = cor_methods(),
                dist.method = dist_methods(),
                cluster.method = cluster_methods(),
                max.dist = 1,
                h = NULL,
                k = NULL,
                min.size = 5,
                max.size = 0.5) {
  .hca(x,
       cor.method = cor.method,
       dist.method = dist.method,
       cluster.method = cluster.method,
       max.dist = max.dist,
       h = h,
       k = k,
       min.size = min.size,
       max.size = max.size
  )
}


#' @rdname hca
#' @param return.steps whether only return the object of specific step, such as correlation or distance matrix. Default: FALSE
#' @param reorder Whether reorder the result. Default: TRUE
#' @param reorder.col,reorder.row Whether reorder the result. if set `reorder`, it will reorder row and column. Or you set them separately.
#' @param ... Other argument of `.hca`
#' @export
hca_cor <- function(x,
                    return.steps = F,
                    reorder = T,
                    reorder.col = reorder,
                    reorder.row = reorder,
                    ...) {
  if (reorder.col | reorder.row) {
    obj <- .hca(x, hclust.end = T, ...)
    if (reorder.col) obj$cr <- obj$cr[, obj$ord]
    if (reorder.row) obj$cr <- obj$cr[obj$ord, ]
  } else {
    obj <- .hca(x, cor.end = T, ...)
  }
  if (return.steps) {
    return(obj)
  }
  obj$cr
}

#' @rdname hca
#' @export
hca_dist <- function(x, return.steps = F, ...) {
  if (return.steps) {
    return(.hca(x, dist.end = T, ...))
  }
  .hca(x, dist.end = T, ...)$dist
}

#' @rdname hca
#' @export
hca_tree <- function(x, return.steps = F, ...) {
  if (return.steps) {
    return(.hca(x, hclust.end = T, ...))
  }
  .hca(x, hclust.end = T, ...)$tree
}

#' @rdname hca
#' @export
hca_order <- function(x, return.steps = F, cor.method = "pearson", ...) {
  if (return.steps) {
    return(.hca(x, hclust.end = T, cor.method = cor.method, ...))
  }
  .hca(x, hclust.end = T, cor.method = cor.method, ...)$order
}

#' @rdname hca
#' @export
hca_groups <- function(x, return.steps = F, ...) {
  if (return.steps) {
    return(.hca(x, ...))
  }
  .hca(x, ...)$groups
}


#' @rdname hca
#' @export
hca_reorder <- function(x,
                        reorder.col = T,
                        reorder.row = T,
                        cor.method = "none",
                        ...) {
  if (reorder.col) x <- x[, .hca(x = x, hclust.end = T, cor.method = cor.method, ...)$order]
  if (reorder.row) x <- x[.hca(x = t(x), hclust.end = T, cor.method = cor.method, ...)$order, ]
  x
}

#' @title Grouped Reorder and Ungroup
#' @description Keep desired groups together but reorder members within each group. E.g. If cells from multiple samples are being reordered, samples can be kept together and cells reordered within. Option to also reorder the order of groups or samples, such that groups with similar average profiles are placed next to one another.
#' @param m matrix to be reordered
#' @param groups groups within which reordering should take place.
#' @param interorder if TRUE, group order itself is reordered such that groups with similar profiles are placed near one another. E.g. Samples with similar average CNA profiles. Default: FALSE
#' @param intraorder if TRUE, group order itself is reordered such that groups with similar profiles are placed near one another. E.g. Samples with similar average CNA profiles. Default: TRUE
#' @param cor.method desired correlation metric. Default: 'pearson'
#' @param dist.method desired distance metric to be used on top of correlation matrix. Default: 'euclidean'
#' @param cluster.method desired agglomeration method. Default: 'average'
#' @param Names return the vector of ordered IDs instead of the reordered matrix. Default: FALSE
#' @return reordered matrix (same dimensions as input) or a character vector of ordered column names if Names = T.
#' @rdname grouped_reorder
#' @export
grouped_reorder <- function(m,
                            groups,
                            interorder = FALSE,
                            intraorder = TRUE,
                            cor.method = 'pearson',
                            dist.method = 'euclidean',
                            cluster.method = 'average',
                            Names = FALSE) {

  if (interorder) {
    groupAvgs <- sapply(groups, function(group) rowMeans(m[, group, drop = F]))
    groups <- groups[hca_order(rowcenter(groupAvgs))]
  }

  if (intraorder) {
    m.list <- sapply(groups, function(x) m[, x], simplify = F)
    groups <- sapply(m.list, function(m) hca_order(rowcenter(m),
                                                   cor.method = cor.method,
                                                   dist.method = dist.method,
                                                   cluster.method = cluster.method))
  }

  if (Names) {
    return(groups)
  }

  ord <- as.character(unlist(groups, use.names = F))
  m[, ord]
}
