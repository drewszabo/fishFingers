# predict_bcf.R
#' Predict bioconcentration factor (BCF)
#'
#' Top-level wrapper to predict BCF from either SMILES strings or a SIRIUS
#' project directory. Internally generates fingerprints, applies the trained
#' XGBoost model, and optionally performs Monte Carlo thresholding.
#'
#' @param x Character vector of SMILES (when input="smiles") or path to SIRIUS project directory (when input="sirius"). If NULL and input="sirius", selects an open SIRIUS project interactively.
#' @param input Character string, either "smiles" or "sirius".
#' @param species Character string specifying the species.
#' @param threshold Character string specifying thresholding approach:
#'   "basic" or "mc".
#' @param topMost logical; if input = "sirius", then select if the top ranked or all candidates are predicted.
#' @param N integer; if threshold = "mc" this number specifies the number of Monte Carlo simulations performed. Higher iterations take longer.
#' @param cutoff double; if threshold = "mc" this ratio specifies the cutoff after MC simulation.
#'.@param charge integer; if input = "sirius", this specifies the charge state of the features to be predicted. Default is 1.
#' 
#' @return A data.frame containing the original input, predicted BCF, and applicability domain.
#' @export
#'
#' @examples
#' # predict BCF from SMILES
#' predict_bcf(
#'   x = c("CCNC1=NC(=NC(=N1)Cl)NC(C)C", "CCNC1=NC(=NC(=N1)Cl)NCC"), # atrazine and simazine
#'   input = "smiles",
#'   species = "Cyprinus carpio"
#' )
predict_bcf <- function(
  x = NULL,
  input = c("smiles", "sirius"),
  species = "Cyprinus carpio",
  threshold = "mc",
  topMost = c("formula", "compound", "all"),
  N = 10000,
  cutoff = 0.5,
  charge = 1
) {

  ## ---- argument validation --------------------------------------------------
  input <- match.arg(input)
  topMost <- match.arg(topMost)

  if (missing(species)) {
    warning("Argument 'species' not provided. Defaulting to Cyprinus carpio.", call. = FALSE)
  }

  if (input == "sirius" && missing(threshold)) {
    warning("Argument 'threshold' not provided. Defaulting to 'mc'.", call. = FALSE)
  }

  check_species(species)

  if (is.null(x) && input == "smiles") {
    stop("Argument 'x' must not be empty.", call. = FALSE)
  }

  if (missing(input)) {
    stop("Argument 'input' must be specified as either 'smiles' or 'sirius'.", call. = FALSE)
  }

  ## ---- model loading --------------------------------------------------------
  model_path <- system.file(
    "extdata",
    "fishFingers.json",
    package = "fishFingers"
  )

  if (model_path == "") {
    stop("Model file 'fishFingers.json' not found in inst/extdata. Try reinstalling fishFingers.",
         call. = FALSE)
  }

  model <- xgboost::xgb.load(model_path)
  #meta <- readRDS(meta_path)

  # Load fpIndex to filter fingerprints to model features
  fp_index_path <- system.file(
    "extdata",
    "fpIndex_v2.0.csv",
    package = "fishFingers"
  )
  fpIndex <- read.csv(fp_index_path, check.names = FALSE)
  fp_names <- fpIndex$fpName[fpIndex$fpType != "ecfp6"]

  species_path <- system.file(
    "extdata",
    "species_index.csv",
    package = "fishFingers"
  )

  species_index <- read.csv(species_path, check.names = FALSE)

  #meta_path <- system.file(
  #  "extdata",
  #  "fishFingers_metadata.rds",
  #  package = "fishFingers"
  #)

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

    input_df <- data.frame(
      input = x,
      stringsAsFactors = FALSE
    )

    fingerprints <- generate_fingerprints(smiles = x)

  } else if (input == "sirius") {

    post_prob_matrix <- read_sirius_fingerprints(
      path = x,
      topMost = topMost,
      charge = charge
    )

    # Convert to data.frame (not data.table) to avoid data.table issues
    post_prob_matrix <- as.data.frame(post_prob_matrix)
    
    # Find missing columns and add them with 0 values
    missing_cols <- setdiff(fp_names, names(post_prob_matrix))
    
    for (col in missing_cols) {
      post_prob_matrix[[col]] <- 0L
    }
    
    # Reorder to metadata columns + fp_names order
    metadata_cols <- names(post_prob_matrix)[1:5]
    post_prob_matrix <- post_prob_matrix[, c(metadata_cols, fp_names), drop = FALSE]

    input_df <- data.frame(
      post_prob_matrix[, 1:5]
    )
  }

  ## ---- thresholding ---------------------------------------------------------
  if (input == "sirius") {

    if (threshold == "mc") {

      fingerprints <- fishFingers_mc(
        x = post_prob_matrix,
        N = N,
        threshold = cutoff
      )

    } else if (threshold == "basic") {
      post_prob_matrix[, fp_names] <- lapply(post_prob_matrix[, fp_names], function(x) {
        round(as.numeric(as.character(x)))
      })
      fingerprints <- as.matrix(post_prob_matrix[, fp_names])
    }

  }


  ## ---- calculate similarity -------------------------------------------------

  tansim <- appdomain(fingerprints)

  interp <- case_when(
    tansim == 1 ~ "match",
    tansim >= 0.8 ~ "good",
    tansim >= 0.6 ~ "fair",
    TRUE ~ "low"
  )

    ## ---- calculate fish n -------------------------------------------------

  train_n <- species_index$n[match(species, species_index$scientific_name)]

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
    bcf_pred = pred,
    conf = tansim,
    conf_interp = interp,
    species_train_n = train_n
  )

  rownames(out) <- NULL
  return(out)
}
