#species.R
#' Check whether a species is supported by the fishFingers model
#'
#' @param species character; scientific or common name
#' @param list logical; if TRUE, return all supported species
#'
#' @return Logical (TRUE/FALSE) or a data.frame of supported species
#' @export
check_species <- function(species = NULL, list = FALSE) {

  species_index <- utils::read.csv(
    system.file("extdata", "species_index.csv", package = "fishFingers"),
    stringsAsFactors = FALSE
  )

  if (isTRUE(list)) {
    return(species_index)
  }

  if (is.null(species) || !nzchar(species)) {
    stop("Please provide a species name, or set list = TRUE", call. = FALSE)
  }

  species <- trimws(tolower(species))

  sci_match <- tolower(species_index$scientific_name) == species

  has_common <- !is.na(species_index$common_name)
  com_match <- rep(FALSE, nrow(species_index))
  com_match[has_common] <-
    tolower(species_index$common_name[has_common]) == species

  in_domain <- any(sci_match | com_match)

  if (!in_domain) {
    warning(
      "Species '", species,
      "' is not in the fishFingers applicability domain.\n",
      "Use check_species(list = TRUE) to see supported species.",
      call. = FALSE
    )
  }

  return(in_domain)
}
