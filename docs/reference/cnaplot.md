# Plot A Heatmap with ggplot2

\`ggheatmap()\` create a heatmap ggplot2 object.

Uses \`ggplot::geom_raster\` to generate a plot of copy number
aberration (CNA) values. Delineates chromosomes and chromosome arms, and
provides options for clustering cells or groups of cells (e.g., cell
types, samples, or subclones) prior to plotting.

## Usage

``` r
ggheatmap(
  dat,
  x = Var1,
  y = Var2,
  fill = value,
  x.num = TRUE,
  y.num = TRUE,
  limits = c(-1, 1),
  limits.find = FALSE,
  limits.symmetric = FALSE,
  cols = NULL,
  na.value = NULL,
  raster.labels = FALSE,
  raster.labels.col = "black",
  x.breaks = waiver(),
  y.breaks = waiver(),
  x.labels = waiver(),
  y.labels = waiver(),
  labels.col = "black",
  x.labels.angle = 0,
  y.labels.angle = 0,
  x.expand = c(0, 0),
  y.expand = c(0, 0),
  x.position = "bottom",
  y.position = "left",
  x.title = NULL,
  y.title = NULL,
  legend.breaks = NULL,
  legend.labels = NULL,
  legend.labels.col = "black",
  legend.width = unit(0.3, "cm"),
  legend.height = unit(0.8, "cm"),
  legend.margin = ggplot2::margin(t = -0.02, unit = "npc"),
  legend.title = "Scaled Expression",
  legend.position = "bottom",
  legend.justification = "right",
  title = waiver(),
  subtitle = waiver(),
  caption = waiver(),
  text.size = 12,
  plot.margin = ggplot2::margin(0.15, 0.15, 0.15, 0.15, "cm")
)

cnaplot(
  cna,
  genome = "hg19",
  reorder = FALSE,
  reorder.by = NULL,
  groups = NULL,
  interorder = TRUE,
  dist.method = "euclidean",
  cluster.method = "complete",
  hide = c("21", "Y"),
  limits = c(-1, 1),
  cols = NULL,
  legend.title = "Scaled CNA values",
  title = "Copy-number aberrations",
  ...
)
```

## Arguments

- dat:

  A `data.frame` or a `tibble` to plot ggplot.

- x:

  The column name of variable mapping to x axis.

- y:

  The column name of variable mapping to y axis.

- fill:

  The column name of variable mapping to heatmap fill.

- x.num:

  description

- y.num:

  description

- limits:

  for the color key; replaces out of bounds values with the nearest
  limit. Default: c(-1, 1)

- limits.find:

  If used, auto find the value limits to plot heatmap. Default: FALSE

- limits.symmetric:

  Used with `limits.find`. If used, set the value limits to be
  symmetrical. Default: FALSE

- cols:

  custom color palette (character vector), Default: NULL

- na.value:

  The colot of NA value.

- raster.labels:

  Whether label the value in heatmap. Default: FALSE

- raster.labels.col:

  The color of the label in heatmap

- x.breaks, y.breaks:

  same as `breaks` in
  [`scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  or
  [`scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)

- x.labels, y.labels:

  same as `labels` in
  [`scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  or
  [`scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)

- labels.col:

  The color of x axis and y axis label

- x.labels.angle, y.labels.angle:

  The Angle of x axis or y axis label (in \\\[0, 360\]\\)

- x.expand, y.expand:

  same as `expand` in
  [`scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  or
  [`scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)

- x.position, y.position:

  same as `position` in
  [`scale_x_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  or
  [`scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)

- x.title, y.title:

  The title of x axis or y axis

- legend.breaks:

  same as `breaks` in
  [`scale_fill_gradientn`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)

- legend.labels:

  same as `labels` in
  [`scale_fill_gradientn`](https://ggplot2.tidyverse.org/reference/scale_gradient.html)

- legend.labels.col:

  The legend label color

- legend.width:

  same as `legend.width` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- legend.height:

  same as `legend.height` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- legend.margin:

  same as `legend.margin` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- legend.title:

  legend title, Default: NULL

- legend.position:

  same as `legend.position` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- legend.justification:

  same as `legend.justification` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- title:

  plot title, Default: 'Copy-number aberrations'

- subtitle:

  same as `subtitle` in
  [`labs`](https://ggplot2.tidyverse.org/reference/labs.html)

- caption:

  same as `caption` in
  [`labs`](https://ggplot2.tidyverse.org/reference/labs.html)

- text.size:

  The base size of text in the plot

- plot.margin:

  same as `plot.margin` in
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html)

- cna:

  CNA matrix (genes by cells). It can be generated using
  [`Infercna`](https://moonerss.github.io/Infercna/reference/Infercna.md).

- genome:

  The genome annotation to order the genes.

- reorder:

  reorder cells using hierarchical clustering, Default: FALSE

- reorder.by:

  reorder cells by genes on select chromosomes or chromosome arms, e.g.,
  c("7", "10", "1p", "19q"). Default: NULL

- groups:

  groups of cells to delineate and label on the plot (e.g., cell types,
  samples, or subclones). If \<reorder\> is TRUE, ordering will be
  performed within groups, Default: NULL

- interorder:

  reorder by hierarchical clustering betweeen groups, Default: TRUE

- dist.method:

  distance metric for reordering, Default: 'euclidean'

- cluster.method:

  linkage method for reordering, Default: 'complete'

- hide:

  hide x-axis labels for specific chromosomes (e.g., small chromosomes
  whose labels would overlap flanking labels), Default: c("21", "Y")

- ...:

  Other argument used in `ggheatmap`

## Value

return a ggplot2 object

return a ggplot object
