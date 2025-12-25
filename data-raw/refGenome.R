## code to prepare `refGenome` dataset goes here
quickmart <- function(genome = "hg19", host = "https://www.ensembl.org") {
  if (genome == "hg19") {
    dataset <- "hsapiens_gene_ensembl"
    host <- "https://grch37.ensembl.org"
  } else if (genome == "hg38") {
    dataset <- "hsapiens_gene_ensembl"
    host <- "https://useast.ensembl.org"
  } else if (genome == "mm39") {
    dataset <- "mmusculus_gene_ensembl"
    host <- "https://useast.ensembl.org"
  } else if (genome == "mm10") {
    dataset <- "mmusculus_gene_ensembl"
    host <- 'https://feb2021.archive.ensembl.org'
  } else {
    stop('genome "', genome, '" not found')
  }
  # Choose which species to use and server to download from
  a <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
            dataset = dataset,
            host = host
  )
}

getGenome <- function(mart, add_attribute = "hgnc_symbol") {
  filters <- c("chromosome_name")
  ## must load genes
  values <- list(c(1:22, "X", "Y"))
  attributes <- c(
    "chromosome_name",
    "band",
    "strand",
    "start_position",
    "end_position",
    "ensembl_gene_id"
  )
  attributes <- unique(c(attributes, add_attribute))
  result <- biomaRt::getBM(
    attributes = attributes,
    mart = mart,
    filters = filters,
    values = values,
    uniqueRows = T
  )
  if ("hgnc_symbol" %in% add_attribute) {
    result <- dplyr::filter(result, hgnc_symbol != "")
  }
  if ("mgi_symbol" %in% add_attribute) {
    result <- dplyr::filter(result, mgi_symbol != "")
  }
  result
}

settingchrlevels <- function(dat) {
  x <- unique(dat$chromosome_name)
  xnum <- sort(as.numeric(x[x %in% as.character(1:100)]))
  xchar <- sort(x[!x %in% xnum])
  chrlev <- c(xnum, xchar)
  dat <- dplyr::mutate(dat, chromosome_name = factor(as.character(chromosome_name), levels = chrlev))
  dat <- dplyr::arrange(dat, chromosome_name, start_position)
  armlev <- unique(dat$arm)
  dat <- dplyr::mutate(dat, arm = factor(as.character(arm), levels = armlev))
  dat
}

## run & get data
hg19 <- quickmart("hg19")
hg38 <- quickmart("hg38")
mm10 <- quickmart("mm10")
mm39 <- quickmart("mm39")

marts <- list(hg19, hg38, mm10, mm39)
add_attribute <- c("hgnc_symbol", "hgnc_symbol", "mgi_symbol", "mgi_symbol")
dats <- Map(getGenome, mart = marts, add_attribute = add_attribute)

dats <- lapply(dats, function(x) {
  colnames(x)[7] <- "symbol"
  x
})

dats[1:2] <- lapply(dats[1:2], function(d) {
  d %>% dplyr::mutate(arm = paste0(chromosome_name, stringr::str_extract(band, "p|q")))
})

dats[3:4] <- lapply(dats[3:4] , function(d) {
  d %>% dplyr::mutate(arm = chromosome_name)
})

dats <- lapply(dats, function(dat) settingchrlevels(dat))

dats <- lapply(dats, function(dat) {
  dat[, c(7, 4, 5, 1, 8, 2, 3, 6)]
})

hg19 <- dats[[1]]
hg38 <- dats[[2]]
mm10 <- dats[[3]]
mm39 <- dats[[4]]


usethis::use_data(hg19, hg38, mm10, mm39, overwrite = TRUE)
