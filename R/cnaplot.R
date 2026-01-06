# prepare the data to plot cna heatmap
# m --> a cna matrix
# reorder --> whether the cells
# reorder.by --> reorder cells by genes on select chromosomes or chromosome arms, e.g., c("7", "10", "1p", "19q"). Default: NULL
# groups --> groups of cells to delineate and label on the plot (e.g., cell types, samples, or subclones). If <reorder> is TRUE, ordering will be performed within groups, Default: NULL
# interorder --> reorder by hierarchical clustering between groups, Default: TRUE
# genome --> set genome to use ('hg19' or 'hg38'), Default: 'hg19'
# dist.method --> distance metric for reordering, Default: 'euclidean'
# cluster.method --> linkage method for reordering, Default: 'complete'

.ggcna_prep <- function(m,
                        reorder = FALSE,
                        reorder.by = NULL,
                        groups = NULL,
                        interorder = TRUE,
                        genome = "hg19",
                        dist.method = "euclidean",
                        cluster.method = "complete") {
  if (!reorder) {
    return(m)
  }

  # null for all gene, not null for select gene
  if (is.null(reorder.by)) {
    genes <- rownames(m)
  } else {
    reorder.by <- as.character(reorder.by)
    useGenome(genome)
    geneList <- c(splitGenes(rownames(m)), splitGenes(rownames(m), by = "arm"))
    genes <- geneList[names(geneList) %in% reorder.by] %>% unlist()
    rm(geneList)
  }

  if (is.null(groups)) {
    cols <- hca_order(m[genes, ], cluster.method = cluster.method, dist.method = dist.method)
  } else {
    cols <- colnames(grouped_reorder(m[genes, ],
                                     groups = groups, interorder = interorder,
                                     cluster.method = cluster.method, dist.method = dist.method
    ))
  }

  m[, cols]
}

.limit <- function(v) {
  m <- 10
  Sign <- sign(v)
  v <- abs(v)
  v.lim <- Sign * (floor(m * v) / m)
  return(v.lim)
}

.limits <- function(v, symmetric = TRUE) {
  stopifnot(length(v) == 2, all(sapply(v, class) == 'numeric'))

  if (!symmetric) {
    return(sapply(v, .limit))
  }

  limit <- .limit(min(abs(v)))
  return(c(-1 * limit, 1 * limit))
}

.axis.spacer <- function(breaks, labels, limits, levels = NULL) {
  if (!is.null(labels) & !is.null(levels)) {
    breaks <- levels %in% labels
    labels <- levels[breaks]
  }
  if (is.null(breaks)) {
    breaks <- seq(limits[[1]], limits[[2]], limits[[2]])
  }
  if (is.null(labels)) {
    labels <- breaks
  }
  return(list(breaks = breaks, labels = labels))
}

# get from ggpubr::rotate_x_text
rotate_x_text <- function(angle = 90, hjust = NULL, vjust = NULL, ...) {
  if (missing(hjust) & angle > 5)
    hjust <- 1
  if (missing(vjust) & angle == 90)
    vjust <- 0.5
  theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust, ...))
}

# get from ggpubr::rotate_y_text
rotate_y_text <- function (angle = 90, hjust = NULL, vjust = NULL, ...) {
  if (missing(hjust) & angle == 90)
    hjust <- 0.5
  else if (missing(hjust) & angle > 5)
    hjust <- 1
  theme(axis.text.y = element_text(angle = angle, hjust = hjust,
                                   vjust = vjust, ...))
}


#' @title Plot A Heatmap with ggplot2
#'
#' @description
#' `ggheatmap()` create a heatmap ggplot2 object.
#'
#' @param dat A \code{data.frame} or a \code{tibble} to plot ggplot.
#' @param x The column name of variable mapping to x axis.
#' @param y The column name of variable mapping to y axis.
#' @param fill The column name of variable mapping to heatmap fill.
#' @param x.num description
#' @param y.num description
#' @param limits A vector to set the min and max value to plot heatmap. Default: -1, 1
#' @param limits.find If used, auto find the value limits to plot heatmap. Default: FALSE
#' @param limits.symmetric Used with \code{limits.find}. If used, set the value limits to be symmetrical. Default: FALSE
#' @param cols a color vector used to plot heatmap. If NULL, automatically set built-in palettes. Default: NULL
#' @param na.value The colot of NA value.
#' @param raster.labels Whether label the value in heatmap. Default: FALSE
#' @param raster.labels.col The color of the label in heatmap
#' @param x.breaks,y.breaks same as \code{breaks} in \code{\link[ggplot2]{scale_x_continuous}} or \code{\link[ggplot2]{scale_y_continuous}}
#' @param x.labels,y.labels same as \code{labels} in \code{\link[ggplot2]{scale_x_continuous}} or \code{\link[ggplot2]{scale_y_continuous}}
#' @param labels.col The color of x axis and y axis label
#' @param x.labels.angle,y.labels.angle The Angle of x axis or y axis label (in \eqn{[0, 360]})
#' @param x.expand,y.expand same as \code{expand} in \code{\link[ggplot2]{scale_x_continuous}} or \code{\link[ggplot2]{scale_y_continuous}}
#' @param x.position,y.position same as \code{position} in \code{\link[ggplot2]{scale_x_continuous}} or \code{\link[ggplot2]{scale_y_continuous}}
#' @param x.title,y.title The title of x axis or y axis
#' @param legend.breaks same as \code{breaks} in \code{\link[ggplot2]{scale_fill_gradientn}}
#' @param legend.labels same as \code{labels} in \code{\link[ggplot2]{scale_fill_gradientn}}
#' @param legend.labels.col The legend label color
#' @param legend.width same as \code{legend.width} in \code{\link[ggplot2]{theme}}
#' @param legend.height same as \code{legend.height} in \code{\link[ggplot2]{theme}}
#' @param legend.margin same as \code{legend.margin} in \code{\link[ggplot2]{theme}}
#' @param legend.title The name of legend
#' @param legend.position same as \code{legend.position} in \code{\link[ggplot2]{theme}}
#' @param legend.justification same as \code{legend.justification} in \code{\link[ggplot2]{theme}}
#' @param title same as \code{title} in \code{\link[ggplot2]{labs}}
#' @param subtitle same as \code{subtitle} in \code{\link[ggplot2]{labs}}
#' @param caption same as \code{caption} in \code{\link[ggplot2]{labs}}
#' @param text.size The base size of text in the plot
#' @param plot.margin same as \code{plot.margin} in \code{\link[ggplot2]{theme}}
#'
#' @return return a ggplot2 object
#'
#' @export
#' @rdname cnaplot
#'
#' @import ggplot2
#' @importFrom rlang :=
#'
ggheatmap <- function(dat,
                      x = Var1,
                      y = Var2,
                      fill = value,
                      x.num = TRUE,
                      y.num = TRUE,
                      limits = c(-1, 1),
                      limits.find = FALSE,
                      limits.symmetric = FALSE,
                      cols = NULL, na.value = NULL,
                      raster.labels = FALSE,
                      raster.labels.col = 'black',
                      x.breaks = waiver(),
                      y.breaks = waiver(),
                      x.labels = waiver(),
                      y.labels = waiver(),
                      labels.col = 'black',
                      x.labels.angle = 0,
                      y.labels.angle = 0,
                      x.expand = c(0, 0),
                      y.expand = c(0, 0),
                      x.position = 'bottom',
                      y.position = 'left',
                      x.title = NULL,
                      y.title = NULL,
                      legend.breaks = NULL,
                      legend.labels = NULL,
                      legend.labels.col = 'black',
                      legend.width = unit(0.3, 'cm'),
                      legend.height = unit(0.8, 'cm'),
                      legend.margin = ggplot2::margin(t = -0.02, unit='npc'),
                      legend.title = 'Scaled Expression',
                      legend.position = 'bottom',
                      legend.justification = 'right',
                      title = waiver(),
                      subtitle = waiver(),
                      caption = waiver(),
                      text.size = 12,
                      plot.margin = ggplot2::margin(0.15, 0.15, 0.15, 0.15,'cm')
                      ) {
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  fill <- rlang::enquo(fill)
  xname <- rlang::quo_name(x)
  yname <- rlang::quo_name(y)

  ## axis breaks and labels
  if (inherits(x.breaks, 'waiver')) x.breaks <- x.breaks
  if (inherits(x.labels, 'waiver')) x.labels <- x.labels
  if (inherits(y.breaks, 'waiver')) y.breaks <- y.breaks
  if (inherits(y.labels, 'waiver')) y.labels <- y.labels


  if (x.num | inherits(dat %>% dplyr::pull(!!x), 'numeric')) {
    dat <- dat %>% dplyr::mutate(!!xname := as.numeric(!!x))
    x.scale.FUN <- ggplot2::scale_x_continuous
  } else {
    x.scale.FUN <- ggplot2::scale_x_discrete
  }

  if (y.num | inherits(dat %>% dplyr::pull(!!y), 'numeric')) {
    dat <- dat %>% dplyr::mutate(!!yname := as.numeric(!!y))
    y.scale.FUN <- ggplot2::scale_y_continuous
  } else {
    y.scale.FUN <- ggplot2::scale_y_discrete
  }

  ## show number
  if (is_true(raster.labels)) {
    dat <- dat %>% dplyr::mutate(label = round(!!fill, 1))
  } else if (!is_null(raster.labels) & length(raster.labels) == nrow(dat)) {
    dat$label <- raster.labels
  }
  ## set limit
  if (limits.find) {
    v <- dat %>% dplyr::select(!!fill) %>% range
    limits <- .limits(v = v, symmetric = limits.symmetric)
  }

  ## set color palettes
  if (is_null(cols) & limits[[1]] >= 0) {
    # The same: cols <- RColorBrewer::brewer.pal(9, 'YlOrRd')
    cols <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A",
              "#E31A1C", "#BD0026", "#800026")
  } else if (is_null(cols)) {
    cols <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0",
              "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
  }

  ## set legend breaks
  legend <- .axis.spacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
  legend.breaks <- legend$breaks
  legend.labels <- legend$labels

  ## plot heatmap
  p <- ggplot(data = dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
    geom_raster() +
    scale_fill_gradientn(colors = cols,
                         limits = limits,
                         oob = scales::squish,
                         breaks = legend.breaks,
                         labels = legend.breaks,
                         name = legend.title,
                         na.value = na.value,
                         guide = ggplot2::guide_colorbar(frame.colour = 'black', ticks.colour='black',
                                                         title.position = 'top', title.hjust = 0.5)
                         ) +
    eval(quote(x.scale.FUN(breaks = x.breaks, labels = x.labels, expand = x.expand, position = x.position))) +
    eval(quote(y.scale.FUN(breaks = y.breaks, labels = y.labels, expand = y.expand, position = y.position)))

  ## add raster label
  if ("label" %in% colnames(dat)) {
    p <-  p + ggplot2::geom_text(aes(label = label), colour = raster.labels.col)
  }

  p <- p + labs(x = x.title, y = y.title, title = title, subtitle = subtitle, caption = caption) +
    theme_bw(base_size = text.size) +
    ## set axis text angle
    rotate_x_text(angle = x.labels.angle) +
    rotate_y_text(angle = y.labels.angle) +
    theme(axis.text = element_text(color = labels.col),
          legend.position = legend.position,
          legend.justification = legend.justification,
          legend.text = ggplot2::element_text(colour = legend.labels.col),
          legend.margin = legend.margin,
          legend.key.height = legend.height,
          legend.key.width = legend.width,
          plot.margin = plot.margin)

  return(p)
}




#' @title plot a CNA heatmap
#' @description
#' Uses `ggplot::geom_raster` to generate a plot of copy number aberration (CNA) values. Delineates chromosomes and chromosome arms, and provides options for clustering cells or groups of cells (e.g., cell types, samples, or subclones) prior to plotting.
#'
#' @param cna CNA matrix (genes by cells). It can be generated using \code{\link[Infercna]{Infercna}}.
#' @param genome The genome annotation to order the genes.
#' @param reorder reorder cells using hierarchical clustering, Default: FALSE
#' @param reorder.by reorder cells by genes on select chromosomes or chromosome arms, e.g., c("7", "10", "1p", "19q"). Default: NULL
#' @param groups groups of cells to delineate and label on the plot (e.g., cell types, samples, or subclones). If <reorder> is TRUE, ordering will be performed within groups, Default: NULL
#' @param interorder reorder by hierarchical clustering betweeen groups, Default: TRUE
#' @param dist.method distance metric for reordering, Default: 'euclidean'
#' @param cluster.method linkage method for reordering, Default: 'complete'
#' @param hide hide x-axis labels for specific chromosomes (e.g., small chromosomes whose labels would overlap flanking labels), Default: c("21", "Y")
#' @param limits for the color key; replaces out of bounds values with the nearest limit. Default: c(-1, 1)
#' @param cols custom color palette (character vector), Default: NULL
#' @param legend.title legend title, Default: NULL
#' @param title plot title, Default: 'Copy-number aberrations'
#' @param ... Other argument used in \code{\link[Infercna]{ggheatmap}}
#'
#' @return return a ggplot object
#'
#' @export
#' @rdname cnaplot
#'
#'
cnaplot <- function(cna,
                    genome = 'hg19',
                    reorder = FALSE,
                    reorder.by = NULL,
                    groups = NULL,
                    interorder = TRUE,
                    dist.method = 'euclidean',
                    cluster.method = 'complete',
                    hide = c('21','Y'),
                    limits = c(-1, 1),
                    cols = NULL,
                    legend.title = 'Scaled CNA values',
                    title = 'Copy-number aberrations',
                    ...
                    ) {
  if (is_null(cols)) {
    cols <- c(
      "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
      "white", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  }

  ## set geom vlines,breaks,labels for chr and chr arms
  chr.size <-  0.1
  arm.size <- 0.05
  genes <- rownames(cna)
  useGenome(genome)
  m <- orderGenes(cna)
  L <- splitGenes(genes)
  xints <- cumsum(lengths(L)) + chr.size
  L2 <- splitGenes(genes, by = 'arm')
  xints2 <- cumsum(lengths(L2)) + arm.size
  xints2 <- xints2[str_detect(names(xints2), "p")]
  x.breaks <- (xints - chr.size) - lengths(L)/2
  x.labels <- names(xints)
  x.labels[names(xints) %in% hide] <- ''

  m <- .ggcna_prep(m,
                   reorder = reorder,
                   reorder.by = reorder.by,
                   groups = groups,
                   genome = genome,
                   interorder = interorder,
                   dist.method = dist.method,
                   cluster.method = cluster.method)
  if (!is_null(groups)) {
    line.size <- 0.1
    ord <- colnames(m)
    groupNames <- flip(Unlist(groups))[ord]
    groupNames <- groupNames[!duplicated(groupNames)]
    yints <- cumsum(lengths(groups)) + line.size
    y.breaks <- (yints-line.size) - lengths(groups)/2
    y.labels <- names(yints)
  } else {
    y.labels <- ggplot2::waiver()
    y.breaks <- ggplot2::waiver()
  }

  d <- reshape2::melt(as.matrix(m)) %>%
    dplyr::mutate(Var1 = as.numeric(factor(as.character(Var1), levels = rownames(m))),
                  Var2 = as.numeric(factor(as.character(Var2), levels = colnames(m))))

  p <- ggheatmap(d, x = Var1, y = Var2, fill = value, title = title,
                 limits = limits, legend.title = legend.title, cols = cols,
                 x.breaks = x.breaks, x.labels = x.labels,
                 y.breaks = y.breaks, y.labels = y.labels,
                 ...) +
    geom_vline(xintercept = xints, color = 'grey20', linetype = 1, linewidth = chr.size) +
    geom_vline(xintercept = xints2, color = 'grey20', linetype = 2, linewidth = arm.size) +
    # scalop::theme_scalop(legend.text.size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
                   axis.ticks.x = ggplot2::element_blank())
  return(p)
}
