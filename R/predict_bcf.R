# predict_bcf.R
#' Predict bioconcentration factor (BCF)
#'
#' \code{predict_bcf} is a top-level wrapper to predict BCF from either SMILES strings
#' or a SIRIUS project directory. Internally generates fingerprints, applies the trained
#' XGBoost model, and optionally performs Monte Carlo thresholding. \code{predict_bcf_fps}
#' is a lower-level function that accepts fingerprints, and is used by \code{predict_bcf}
#' to perform the calculations.
#'
#' @param x Character vector of SMILES or a file path to a SIRIUS project.
#' @param input Character string, either "smiles" or "sirius".
#' @param species Character string specifying the species.
#' @param threshold Character string specifying thresholding approach:
#'   "basic" or "mc".
#' @param topMost logical; if input = "sirius", then select if the top ranked or all candidates are predicted.
#' @param N integer; if threshold = "mc" this number specifies the number of Monte Carlo simulations performed. Higher iterations take longer.
#' @param cutoff double; if threshold = "mc" this ratio specifies the cutoff after MC simulation.
#' @param fingerprints A \code{data.frame} with fingerprints for prediction.
#'
#' @return A data.frame containing the original input and predicted BCF.
#' @export
#'
#' @examples
#' predict_bcf(
#'   x = c("CCO", "CC(=O)O"),
#'   input = "smiles",
#'   species = "Cyprinus carpio"
#' )
predict_bcf <- function(
    x,
    input = c("smiles", "sirius"),
    species = "Cyprinus carpio",
    threshold = c("mc"),
    topMost = TRUE,
    N = 10000,
    cutoff = 0.5
) {

  ## ---- argument validation --------------------------------------------------
  input <- match.arg(input)
  threshold <- match.arg(threshold)


  if (missing(x) || length(x) == 0) {
    stop("Argument 'x' must not be empty.", call. = FALSE)
  }


  ## ---- fingerprint generation ----------------------------------------------
  if (input == "smiles") {

    if (!is.character(x)) {
      stop("For input = 'smiles', x must be a character vector.",
           call. = FALSE)
    }

    valid <- vapply(
      x,
      webchem::is.smiles,
      logical(1)
    )

    if (!all(valid)) {
      stop(
        "Invalid SMILES detected at positions: ",
        paste(which(!valid), collapse = ", "),
        call. = FALSE
      )
    }

    fingerprints <- generate_fingerprints(smiles = x)

  } else if (input == "sirius") {

    if (!is.character(x) || length(x) != 1L) {
      stop("For input = 'sirius', x must be a single file path.",
           call. = FALSE)
    }

    if (!dir.exists(x)) {
      stop("SIRIUS project directory does not exist: ", x,
           call. = FALSE)
    }

    fingerprints <- read_sirius_fingerprints(
      sirius_project_dir = x,
      topMost = topMost
    )
  }
  
  return(predict_bcf_fps(
    input = input,
    fingerprints = fingerprints,
    species = species,
    threshold = threshold,
    topMost = topMost,
    N = N,
    cutoff = cutoff
  ))
}

#' @rdname predict_bcf
#' @export
predict_bcf_fps <- function(
    input = c("smiles", "sirius"),
    fingerprints,
    species = "Cyprinus carpio",
    threshold = c("mc"),
    topMost = TRUE,
    N = 10000,
    cutoff = 0.5
) {

  ## ---- argument validation --------------------------------------------------
  input <- match.arg(input)
  threshold <- match.arg(threshold)

  if (missing(species)) {
    warning("Argument 'species' not provided. Defaulting to Cyprinus carpio.", call. = FALSE)
  }

  if (input == "sirius" && missing(threshold)) {
    warning("Argument 'threshold' not provided. Defaulting to 'mc'.", call. = FALSE)
  }

  check_species(species)

  ## ---- model loading --------------------------------------------------------
  model_path <- system.file(
    "extdata",
    "fishFingers.json",
    package = "fishFingers"
  )

  #meta_path <- system.file(
  #  "extdata",
  #  "fishFingers_metadata.rds",
  #  package = "fishFingers"
  #)

  if (model_path == "") {
    stop("Model file 'fishFingers.json' not found in inst/extdata. Try reinstalling fishFingers.",
         call. = FALSE)
  }

  model <- xgboost::xgb.load(model_path)
  #meta <- readRDS(meta_path)

  ## ---- fingerprint generation ----------------------------------------------
  if (input == "smiles") {

    input_df <- data.frame(
      input = x,
      stringsAsFactors = FALSE
    )

  } else if (input == "sirius") {

    input_df <- data.frame(
      fingerprints[,1:2]
    )
  }

  ## ---- thresholding ---------------------------------------------------------
  if (input == "sirius") {

    if (threshold == "mc") {

      fingerprints <- fishFingers_mc(
        x = fingerprints,
        N = N,
        threshold = cutoff
      )

    } else if (threshold == "basic") {
      # no transformation; keep raw posterior probabilities
      fingerprints <- as.numeric(fingerprints)
    }

  }


  ## ---- combine species ------------------------------------------------------

  fish <- make_species_matrix(species)

  fish_rep <- fish[rep(1, nrow(fingerprints)), , drop = FALSE]

  fishFinger <- cbind(fingerprints, fish_rep)

  ## ---- prediction -----------------------------------------------------------
  new_x <- as.matrix(fishFinger)
  storage.mode(new_x) <- "numeric"

  dmat <- xgboost::xgb.DMatrix(new_x)

  pred <- predict(
    model,
    newdata = dmat
  )

  ## ---- output ---------------------------------------------------------------
  out <- cbind(
    input_df,
    bcf_pred = pred
  )

  rownames(out) <- NULL
  return(out)
}
