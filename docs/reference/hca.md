# Hierarchical clustering analysis

Hierarchical clustering analysis

## Usage

``` r
hca(
  x,
  cor.method = cor_methods(),
  dist.method = dist_methods(),
  cluster.method = cluster_methods(),
  max.dist = 1,
  h = NULL,
  k = NULL,
  min.size = 5,
  max.size = 0.5
)

hca_cor(
  x,
  return.steps = F,
  reorder = T,
  reorder.col = reorder,
  reorder.row = reorder,
  ...
)

hca_dist(x, return.steps = F, ...)

hca_tree(x, return.steps = F, ...)

hca_order(x, return.steps = F, cor.method = "pearson", ...)

hca_groups(x, return.steps = F, ...)

hca_reorder(x, reorder.col = T, reorder.row = T, cor.method = "none", ...)
```

## Arguments

- x:

  a matrix (features X observations), an object of class \`dist\` or an
  object of class \`hclust\`.

- cor.method:

  option if \<x\> is a matrix. One of 'pearson' (default), 'kendall',
  'spearman' or 'none' (no correlation coefficents computed). Default:
  cor_methods()

- dist.method:

  option if \<x\> is a matrix or correlation matrix. One of 'euclidean'
  (default), 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'
  or 'none' (in which case stats::as.dist() will be used). Default:
  dist_methods()

- cluster.method:

  one of "average" (default), "complete", "single", "ward.D", "ward.D2",
  "mcquitty", "median" or "centroid". Default: cluster_methods()

- max.dist:

  maximum distance between observations. Default: 1

- h:

  height(s) at which to cut the tree. If NULL, h will be set to all tree
  heights. Default: NULL

- k:

  number of groups to return from tree. If \<h\> and \<k\> are both not
  NULL, \<k\> takes precedence. Default: NULL

- min.size:

  minimum allowed cluster/group size. Values between 0 and 1 are
  interpreted as fractions of total count. Groups smaller than
  \<min.size\> are filtered out. Default: 5

- max.size:

  maximum allowed cluster/group size. Values between 0 and 1 are
  interpreted as fractions of total count. Groups larger than
  \<max.size\> are filtered.out. Default: 0.5

- return.steps:

  whether only return the object of specific step, such as correlation
  or distance matrix. Default: FALSE

- reorder:

  Whether reorder the result. Default: TRUE

- reorder.col, reorder.row:

  Whether reorder the result. if set \`reorder\`, it will reorder row
  and column. Or you set them separately.

- ...:

  Other argument of \`.hca\`

## Value

list
