#' hg19 reference genome gene position annotation.
#'
#' A data.frame with gene position inforamtion get from hg19 reference genome.
#'
#' @format A data.frame with gene annotation information.
#' @source \url{https://grch37.ensembl.org}
"hg19"


#' hg38 reference genome gene position annotation.
#'
#' A data.frame with gene position inforamtion get from hg38 reference genome.
#'
#' @format A data.frame with gene annotation information.
#' @source \url{https://useast.ensembl.org}
"hg38"


#' mm10 reference genome gene position annotation.
#'
#' A data.frame with gene position inforamtion get from mm10 reference genome.
#'
#' @format A data.frame with gene annotation information.
#' @source \url{https://feb2021.archive.ensembl.org}
"mm10"


#' mm39 reference genome gene position annotation.
#'
#' A data.frame with gene position inforamtion get from mm39 reference genome.
#'
#' @format A data.frame with gene annotation information.
#' @source \url{https://useast.ensembl.org}
"mm39"

#' Normal Cell IDs in a cohort of 28 Glioblastoma samples.
#'
#' A list of two character vectors containing cell IDs of the normal cells that were found in the GBM study.
#'
#' @format A list of length two containing cell IDs of two normal cell types.
#' \describe{
#'   \item{oligodendrocytes}{219 normal oligodendrocyte cells in the GBM cohort}
#'   \item{macrophages}{707 normal macrophage cells in the GBM cohort}
#' }
#' @source \url{portals.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma}
"refCells"

#' scRNA-seq data from a patient with Glioblastoma
#'
#' A single-cell RNA-sequencing dataset containing 8556 HQ genes and 1265 high quality cells. The data was generated using the Smart-Seq2 protocol and was analysed in the form TPM.
#'
#' @format A matrix with 8556 genes (rows) and 1265 cells (columns):
#' @import Matrix
#' @source \url{https://portals.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma}
"bt771"

#' scRNA-seq data from a patient with Glioblastoma
#'
#' A single-cell RNA-sequencing dataset containing 8556 HQ genes and 1266 high quality cells. The data was generated using the Smart-Seq2 protocol and was analysed in the form TPM.
#'
#' @format A matrix with 8556 genes (rows) and 1266 cells (columns):
#' @import Matrix
#' @source \url{https://portals.broadinstitute.org/single_cell/study/SCP393/single-cell-rna-seq-of-adult-and-pediatric-glioblastoma}
"mgh125"
