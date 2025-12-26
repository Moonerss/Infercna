#' @title List of available genomes
#' @description See which genomes are available for use
#' @return character vector
#' @details more genomes will be added in future dev patches.
#' @rdname availableGenomes
#' @export
availableGenomes <- function() {
  avai <- paste(c("hg19", "hg38", "mm10", "mm39"), collapse = "\n")
  message("Available genomes:\n", avai)
}


#' @title Select a Genome for infercna to use
#' @description Select your genome of choice at the start of an analysis. The available genomes in the current implementation are, for H.sapiens, hg38 (latest) and hg19 (preceding) and for mouse, mm10. You can see which genomes are available via availableGenomes().
#' @param name a character string of genome name. One of 'hg19', 'hg38', 'mm10', 'mm39'.
#' @return genome variables are set internally. No return.
#' @rdname useGenome
#' @export
useGenome <- function(name) {
  message("Genome has been set to ", name)
  return(.configureGenome(name = name))
}

#' @title List the current genome name
#' @description E.g. "hg19" if the genome being used is hg19
#' @return string
#' @details Default genome is "hg19"
#' @rdname currentGenome
#' @export
currentGenome <- function() {
  message("Genome: ", Genv$name)
}

#' @title Retrieve genome data
#' @description Returns a tibble dataframe of the current genome in use. If on function call a character string is supplied that corresponds to a genome in the package, that genome will instead be returned.
#' @param name a genome name. Currently one of 'hg19' (human), 'hg38' (latest human), 'mm10' (mouse), 'mm39' (latest mouse). Default: NULL
#' @return a tibble
#' @seealso
#'  \code{\link[tibble]{as_tibble}}
#' @rdname retrieveGenome
#' @export
#' @importFrom tibble as_tibble
#' @importFrom cli cli_alert_info
#'
retrieveGenome <- function(name = NULL) {
  if (is.null(name)) {
    name <- Genv$name
    data <- Genv$data
  } else {
    data <- .fetchGenome(name)
  }

  cli::cli_alert_info("Retrieving: {.val {name}}")
  tibble::as_tibble(data)
}


#' @title Add your own Genome
#' @description Add your own Genome for infercna to use. The
#' @param genome genome data provided as a dataframe. The dataframe should contain columns 'symbol', 'start_position', 'end_position', 'chromosome_name', 'arm'. The columns 'chromosome_name' and 'arm' should be factors, with the chromosome/chromosome arms ordered correctly.
#' @param name a character string; the name of your genome. Default: 'userDefined'
#' @return no return value. The genome variables will be updated internally.
#' @rdname addGenome
#' @export
#' @importFrom dplyr arrange
#' @importFrom stats setNames
#' @importFrom cli cli_abort cli_alert_info
addGenome <- function(genome, name = "userDefined") {
  columns <- c("symbol", "start_position", "chromosome_name", "arm")

  if (!all(columns %in% colnames(genome))) {
    cli::cli_abort("Columns must include: {.val {columns}}")
  }

  if (is.null(levels(genome$chromosome_name))) {
    cli::cli_abort("Please add levels to the chromosome_name column")
  }

  if (is.null(levels(genome$arm))) {
    cli::cli_abort("Please add levels to the arm column")
  }

  genome <- dplyr::arrange(genome, chromosome_name, start_position)
  .configureGenome(data = genome, name = name)
  cli::cli_alert_info("Genome has been set to {.val {name}}")
}


.fetchGenome <- function(name) {
  stopifnot(is.character(name))
  stopifnot(length(name) == 1)
  result <- try(utils::getFromNamespace(x = name, ns = "Infercna"))
  if (inherits(result, "try-error")) {
    cli::cli_abort("Genome name {.val {name}} does not exist in package data.")
  }
  result
}

.parseGenome <- function(data, name) {
  vars <- as.list(as.data.frame(data))
  nam <- c("symbol", "chromosome_name", "arm", "start_position", "end_position")
  nam2 <- c("gene", "chr", "arm", "start", "end")
  vars <- stats::setNames(vars[nam], nam2)
  vars[2:5] <- sapply(vars[2:5], stats::setNames, vars[["gene"]], simplify = F)
  c(list(name = name, data = data), vars)
}

.configureGenome <- function(name = NULL, data = NULL) {
  if (is.null(data) & is.null(name)) {
    cli::cli_abort("Please provide {.var name} or/and {.var data}.")
  } else if (is.null(data)) {
    data <- .fetchGenome(name)
  } else if (is.null(name)) {
    name <- "userDefined"
  }
  vars <- .parseGenome(data = data, name = name)
  Genv <- new.env(parent = .GlobalEnv)
  Genv <- list2env(vars, envir = Genv)
  assign("Genv", Genv, envir = .GlobalEnv)
}


