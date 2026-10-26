#' appdomain.R

#' Calculate the applicability domain of predicted substances based on Tanimoto similarity
#'
#' @param fingerprints matrix or data.frame; binary structural fingerprints (typically from predict_bcf())
#' @return numeric vector; maximum Tanimoto similarity (0-1) for each observation compared to training set chemicals
#' @details
#' Calculates the maximum Tanimoto similarity between each observation's fingerprints and all chemicals
#' in the training set. Tanimoto similarity ranges from 0 (completely dissimilar) to 1 (identical).
#' Higher values indicate that an observation is more similar to known training chemicals.
#' @export
appdomain <- function(fingerprints) {
  
  if (missing(fingerprints)) {
    stop("Argument 'fingerprints' must be provided.", call. = FALSE)
  }

  # Load training data
  rawdata_path <- system.file(
    "extdata",
    "train_df.csv",
    package = "fishFingers"
  )

  rawdata <- read.csv(rawdata_path, check.names = FALSE)

  fpFinal <- read.csv(system.file("extdata", "final_fps.csv", package = "fishFingers"), check.names = FALSE)
  fp_names <- fpFinal$fpName

  # Extract training fingerprints
  train_fingerprints <- as.matrix(rawdata[, fp_names])

  # Ensure input fingerprints are numeric matrix
  if (!is.matrix(fingerprints)) {
    fingerprints <- as.matrix(fingerprints)
  }

  storage.mode(fingerprints) <- "numeric"
  storage.mode(train_fingerprints) <- "numeric"

  # Calculate max Tanimoto for each observation
  max_similarities <- apply(fingerprints, 1, calc_max_tanimoto, training_matrix = train_fingerprints)

  return(max_similarities)
}