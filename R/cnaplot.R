.ggcna_prep <- function(m, reorder = FALSE, reorder.by = NULL, groups = NULL, interorder = TRUE,
                        genome = "hg38", dist.method = "euclidean", cluster.method = "complete") {
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


## re change the value of vector
.limit <- function(v) {
  m <- 10
  Sign <- sign(v)
  v <- abs(v)
  v.lim <- Sign * (floor(m * v) / m)
  return(v.lim)
}


.limits <- function(v, symmetric = TRUE) {
  stopifnot(length(v) == 2, all(sapply(v, class) == "numeric"))

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

rotate_x_text <- function(angle = 90, hjust = NULL, vjust = NULL, ...) {
  if (missing(hjust) & angle > 5)
    hjust <- 1
  if (missing(vjust) & angle == 90)
    vjust <- 0.5
  theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust, ...))
}

#'
#' Plot A Heatmap with ggplot2
#'
#'
#' `ggheatmap()` create a heatmap ggplot2 object.
#'
#'
#' @param dat A \code{data.frame} or a \code{tibble} to plot ggplot.
#' @param x The column name of variable mapping to x axis.
#' @param y The column name of variable mapping to y axis.
#' @param fill The column name of variable mapping to heatmap fill.
#' @param limits A vector to set the max and min value to plot heatmap.
#' @param raster.labels Whether label the value in heatmap. Defualt: 45.
#' @param raster.label.col The color of raster.labels in heatmap. Default: black.
#' @param lim.find If you want to auto set limit, set this argument as TRUE.
#' @param lim.sym Set with \code{lim.find}, whether use the same upper or lower threshold
#' @param title Heatmap title.
#' @param subtitle Heatmap subtitle.
#' @param caption Heatmap caption title.
#' @param x.name The x axis name.
#'
#' @param y.name The y axis name.
#'
#' @param x.angle The angle to change of x axis text. Default: TRUE
#'
#' @param axis.rel The relative size of axis name. Default: 1
#'
#' @param title.rel The relative size of title name. Default: 1.1
#'
#' @param text.size base font size, given in pts.
#'
#' @param aspect.ratio aspect ratio of the panel
#'
#' @param cols a vector of color to use in heatmap
#'
#' @param na.value Colour to use for missing values
#'
#' @param plot.margin margin around entire plot \code{unit} with the sizes of the `top`, `right`, `bottom`, and `left` margins
#'
#' @param legend.position the default position of legends `none`, `left`, `right`, `bottom`, `top`, `inside`
#'
#' @param legend.justification anchor point for positioning legend inside plot `center` or two-element numeric vector or the justification according to the plot area when positioned outside the plot
#'
#' @param legend.margin the margin around each legend \code{\link[ggplot2]{element}}
#'
#' @param legend.title The legend title name
#'
#' @param legend.title.position placement of legend title relative to the main legend. `top`, `right`, `bottom` or `left`.
#'
#' @param legend.title.rel The relative size of legend title. Default: 0.8
#'
#' @param legend.key.height,legend.key.width size of legend keys \code{unit}
#'
#' @param legend.text.rel The relative size of legend text. Default: 0.8
#'
#' @param legend.text.color The color of legend text. Default: black
#'
#' @param legend.breaks The breaks of legend. One of:
#'
#' - NULL for no breaks
#'
#' - \code{waiver()} for the default breaks computed by the \code{\link[scales]{new_transform}}
#'
#' - A numeric vector of positions
#'
#' - A function that takes the limits as input and returns breaks as output e.g., a function returned by \code{\link[scales]{extended_breaks}}. Note that for position scales, limits are provided after scale expansion. Also accepts rlang \code{\link[rlang]{as_function}} function notation.
#'
#'  @param legend.labels The labels of legend. One of:
#'
#' - NULL for no labels
#'
#' - \code{waiver()} for the default labels computed by the transformation object
#'
#' - A character vector giving labels, must be same length as \code{breaks}
#'
#' - An expression vector, must be the same length as breaks. See ?plotmath for details.
#'
#' - A function that takes the breaks as input and returns labels as output. Also accepts rlang \code{\link[rlang]{as_function}} function notation.
#'
#' @param ticks.linewidth The ticks width. Default: 0.5
#' @param x.axis.position,y.axis.position For position scales, The position of the axis. `top` or `bottom` for x axes, \code{left} or \code{right} for y axes.
#'
#' @param breaks,x.breaks,y.breaks The breaks of x axis or y axis. If you set \code{breaks}, it will be use in \code{x.breaks} and \code{y.breaks}. Or you set them separately. One of:
#'
#' - NULL for no breaks.
#'
#' - \code{waiver()} for the default breaks computed by the \code{\link[scales]{new_transform}}
#'
#' - A numeric vector of positions
#'
#' - A function that takes the limits as input and returns breaks as output e.g., a function returned by \code{\link[scales]{extended_breaks}}. Note that for position scales, limits are provided after scale expansion. Also accepts rlang \code{\link[rlang]{as_function}} function notation.
#'
#' @param labels,x.labels,y.labels The labels of x axis or y axis. If you set \code{breaks}, it will be use in \code{x.breaks} and \code{y.breaks}. Or you set them separately. One of:
#'
#' - NULL for no breaks
#'
#' - \code{waiver()} for the default labels computed by the \code{\link[scales]{new_transform}}
#'
#' - A character vector giving labels, must be same length as \code{breaks}
#'
#' - An expression vector, must be the same length as breaks. See ?plotmath for details.
#'
#' - A function that takes the breaks as input and returns labels as output. Also accepts rlang \code{\link[rlang]{as_function}} function notation.
#'
#' @param num,x.num,y.num Whether the set x axis or y axis as numeric, if you set \code{num}, it will be use in \code{x.num} and \code{y.num}. Or you set them separately. Default: TRUE
#'
#' @param expand For position scales, a vector of range expansion constants used to add some padding around the data to ensure that they are placed some distance away from the axes.
#' Use the convenience function \code{expansion()} to generate the values for the expand argument. The defaults are to expand the scale by 5\% on each side for continuous variables,
#' and by 0.6 units on each side for discrete variables.
#'
#' @importFrom rlang enquo quo_name :=
#' @importFrom scales squish
#' @importFrom ggplot2 aes waiver
#'
#' @return return a ggplot object
#'
ggheatmap <- function(dat, x = Var1, y = Var2, fill = value, limits = c(-0.5, 0.5),
                      raster.labels = FALSE,
                      raster.label.col = "black",
                      lim.find = F, lim.sym = T,
                      title = waiver(), subtitle = waiver(), caption = waiver(),
                      x.name = NULL, y.name = NULL, x.angle = 45,
                      axis.rel = 1, title.rel = 1.1, text.size = 12,
                      aspect.ratio = NULL,
                      cols = NULL, na.value = "white",
                      plot.margin = margin(0, 0.15, 0, 0, "cm"),
                      legend.position = "bottom",
                      legend.justification = "right",
                      legend.margin = margin(0, 0.15, 0, 0, "cm"),
                      legend.title = NULL,
                      legend.title.position = "top",
                      legend.title.rel = 0.8,
                      legend.key.height = grid::unit(0.4, "cm"),
                      legend.key.width = grid::unit(0.6, "cm"),
                      legend.text.rel = 0.8,
                      legend.text.color = "black",
                      legend.breaks = NULL,
                      legend.labels = NULL,
                      ticks.linewidth = 0.5,
                      x.axis.position = "bottom",
                      y.axis.position = "left",
                      breaks = waiver(),
                      x.breaks = waiver(),
                      y.breaks = waiver(),
                      labels = waiver(),
                      x.labels = waiver(),
                      y.labels = waiver(),
                      num = T,
                      y.num = num,
                      x.num = num,
                      expand = c(0, 0)) {
  ## get x y
  x <- enquo(x)
  y <- enquo(y)
  fill <- enquo(fill)
  xname <- quo_name(x)
  yname <- quo_name(y)

  if (!inherits(breaks, "waiver")) {
    x.breaks <- breaks
    y.breaks <- breaks
  }

  if (!inherits(labels, "waiver")) {
    x.labels <- labels
    y.labels <- labels
  }

  ## set scale
  if (x.num | num | inherits(dat %>% dplyr::pull(!!x), "numeric")) {
    dat <- dat %>% dplyr::mutate(!!xname := as.numeric(!!x))
    x.scale.FUN <- ggplot2::scale_x_continuous
  } else {
    x.scale.FUN <- ggplot2::scale_x_discrete
  }
  if (y.num | num | inherits(dat %>% dplyr::pull(!!y), "numeric")) {
    dat <- dat %>% dplyr::mutate(!!yname := as.numeric(!!y))
    y.scale.FUN <- ggplot2::scale_y_continuous
  } else {
    y.scale.FUN <- ggplot2::scale_y_discrete
  }
  x.scale <- quote(x.scale.FUN(
    expand = expand,
    breaks = x.breaks,
    labels = x.labels,
    position = x.axis.position
  ))
  y.scale <- quote(y.scale.FUN(
    expand = expand,
    breaks = y.breaks,
    labels = y.labels,
    position = y.axis.position
  ))

  ## label value
  if (isTRUE(raster.labels)) {
    dat <- dat %>% dplyr::mutate(label = round(!!fill, 1))
  } else if (!is.null(raster.labels) & length(raster.labels) == nrow(dat)) {
    dat$label <- raster.labels
  }

  ## set value limit
  if (lim.find) {
    v <- dat %>%
      dplyr::select(!!fill) %>%
      range()
    limits <- .limits(v = v, symmetric = lim.sym)
  }

  ## set color
  if (is.null(cols) & limits[[1]] >= 0) {
    cols <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C",
              "#FC4E2A","#E31A1C", "#BD0026", "#800026")
  } else if (is.null(cols)) {
    cols <- c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7",
              "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f")
  }


  legend <- .axis.spacer(breaks = legend.breaks, labels = legend.labels, limits = limits)
  legend.breaks <- legend$breaks
  legend.labels <- legend$labels

  ## plot heatmap
  G <- ggplot2::ggplot(dat, aes(x = !!x, y = !!y, fill = !!fill, group = 1)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_gradientn(
      colors = cols,
      limits = limits,
      expand = expand,
      oob = squish,
      breaks = legend.breaks,
      labels = legend.breaks,
      name = legend.title,
      na.value = na.value,
      guide = ggplot2::guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = 3,
        barheight = 0.75
      )
    ) +
    eval(x.scale) +
    eval(y.scale) +
    ggplot2::labs(x = x.name, y = y.name, title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_bw(base_size = text.size) +
    rotate_x_text(angle = x.angle, vjust = 1)
  ggplot2::theme(
    aspect.ratio = aspect.ratio,
    title = ggplot2::element_text(size = ggplot2::rel(title.rel)),
    axis.title = ggplot2::element_text(size = ggplot2::rel(axis.rel)),
    legend.position = legend.position,
    legend.justification = legend.justification,
    legend.text = ggplot2::element_text(size = ggplot2::rel(legend.text.rel), color = legend.text.color, hjust = 0.5),
    legend.title = ggplot2::element_text(size = ggplot2::rel(legend.title.rel)),
    legend.margin = ggplot2::margin(t = -0.5, unit = "cm"),
    plot.margin = plot.margin,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width
  )


  if ("label" %in% colnames(dat)) {
    G <- G +
      ggplot2::geom_text(aes(label = label), colour = raster.label.col)
  }

  return(G)
}




#' @title Plot a CNA heatmap
#' @description Uses `ggplot::geom_raster` to generate a plot of copy number aberration (CNA) values. Delineates chromosomes and chromosome arms, and provides options for clustering cells or groups of cells (e.g., cell types, samples, or subclones) prior to plotting.
#' @param m CNA matrix (genes by cells). <m> can be generated using `Infercna::Infercna`.
#' @param reorder reorder cells using hierarchical clustering, Default: FALSE
#' @param reorder.by reorder cells by genes on select chromosomes or chromosome arms, e.g., c("7", "10", "1p", "19q"). Default: NULL
#' @param groups groups of cells to delineate and label on the plot (e.g., cell types, samples, or subclones). If <reorder> is TRUE, ordering will be performed within groups, Default: NULL
#' @param interorder reorder by hierarchical clustering betweeen groups, Default: TRUE
#' @param genome set genome to use ('hg19' or 'hg38'), Default: 'hg19'
#' @param dist.method distance metric for reordering, Default: 'euclidean'
#' @param cluster.method linkage method for reordering, Default: 'complete'
#' @param hide hide x-axis labels for specific chromosomes (e.g., small chromosomes whose labels would overlap flanking labels), Default: c("21", "Y")
#' @param limits for the colour key; replaces out of bounds values with the nearest limit. Default: c(-1, 1)
#' @param y.angle y-axis labels' angle, Default: 90
#' @param legend.title legend title, Default: NULL
#' @param title plot title, Default: 'Copy-number aberrations'
#' @param axis.text.size x and y axes label size, Default: 12
#' @param legend.height legend bar height, Default: 0.3
#' @param legend.width legend bar width, Default: 0.5
#' @param cols custom colour palette (character vector), Default: NULL
#' @param ... other arguments to pass to `ggheatmap`.
#'
#' @importFrom stringr str_detect
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @return a `ggplot` object
#' @examples
#' \dontrun{
#' m <- Infercna(mgh125, isLog = TRUE, refCells = refCells)
#' malCells <- list(Malignant = setdiff(colnames(mgh125), unlist(refCells)))
#' groups <- c(malCells, refCells)
#' p <- ggcna(m, reorder = T, groups = refCells)
#' }
#' @rdname ggcna
#' @export
cnaplot <- function(m, reorder = FALSE, reorder.by = NULL, groups = NULL, interorder = TRUE,
                    genome = "hg19", dist.method = "euclidean", cluster.method = "complete",
                    hide = c("21", "Y"), limits = c(-1, 1), y.angle = 90,
                    legend.title = NULL, title = "Copy-number aberrations",
                    axis.text.size = 12, legend.height = grid::unit(0.3, "cm"), legend.width = grid::unit(0.5, "cm"),
                    cols = NULL, ...) {
  ## set default color
  if (is.null(cols)) {
    cols <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
              "white", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
  }
  # geom vlines,breaks,labels for chr and chr arms
  chr.size <- 0.1
  arm.size <- 0.05
  genes <- rownames(m)
  useGenome(genome)
  m <- orderGenes(m)
  ## chromosome split
  L <- splitGenes(genes, by = "chr")
  xints <- cumsum(lengths(L)) + chr.size

  ## bind split
  L2 <- splitGenes(genes, by = "arm")
  xints2 <- cumsum(lengths(L2)) + arm.size
  xints2 <- xints2[str_detect(names(xints2), "p")]

  breaks <- (xints - chr.size) - lengths(L) / 2

  labels <- names(xints)
  labels[names(xints) %in% hide] <- ""

  ## change gene order by cluster
  m <- .ggcna_prep(m,
                   reorder = reorder,
                   reorder.by = reorder.by,
                   groups = groups,
                   genome = genome,
                   interorder = interorder,
                   dist.method = dist.method,
                   cluster.method = cluster.method
  )

  if (!is.null(groups)) {
    line.size <- 0.1
    ord <- colnames(m)
    groupNames <- flip(Unlist(groups))[ord]
    groupNames <- groupNames[!duplicated(groupNames)]
    yints <- cumsum(lengths(groups)) + line.size
    ybreaks <- (yints - line.size) - lengths(groups) / 2
    ylabels <- names(yints)
  } else {
    ylabels <- ggplot2::waiver()
    ybreaks <- ggplot2::waiver()
  }

  d <- as.data.frame(m) %>%
    rownames_to_column(var = "Var1") %>%
    pivot_longer(!Var1, names_to = "Var2", values_to = "value") %>%
    dplyr::mutate(
      Var1 = as.numeric(factor(as.character(Var1), levels = rownames(m))),
      Var2 = as.numeric(factor(as.character(Var2), levels = colnames(m)))
    )

  ## plot
  G <- ggheatmap(d,
                 x = Var1, y = Var2, fill = value, limits = limits,
                 legend.title = legend.title, legend.key.height = legend.height,
                 legend.key.width = legend.width, ...
  ) +
    ggplot2::geom_vline(xintercept = xints, col = "grey20", linetype = 1, size = chr.size) +
    ggplot2::geom_vline(xintercept = xints2, col = "grey20", linetype = 2, size = arm.size) +
    ggplot2::scale_x_continuous(breaks = breaks, labels = labels, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = ybreaks, labels = ylabels, expand = c(0, 0)) +
    theme_cna(legend.text.size = 14) +
    ggplot2::theme(
      plot.margin = margin(.2, .2, .2, .2, "cm"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
      axis.text = ggplot2::element_text(size = axis.text.size),
      axis.ticks.x = ggplot2::element_blank(),
      legend.key.height = legend.height,
      legend.key.width = legend.width
    ) +
    ggplot2::scale_fill_gradientn(
      name = legend.title,
      oob = scales::squish,
      colors = cols,
      breaks = c(limits[1], 0, limits[2]),
      limits = limits,
      guide = ggplot2::guide_colorbar(frame.colour = "black", ticks.colour = "black")
    ) +
    ggplot2::labs(title = title)

  if (!is.null(groups)) {
    G <- G +
      ggplot2::geom_hline(yintercept = yints, col = "grey20", size = line.size, linetype = 1) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(angle = y.angle, hjust = 0.5),
        axis.ticks.y = ggplot2::element_blank()
      )
  }

  return(G)
}
