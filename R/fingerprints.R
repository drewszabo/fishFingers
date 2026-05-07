#' R/fingerprints.R
#'
#' Calculate fingerprints from SMILES
#'
#' @param smiles character; vector of SMILES
#' @return data.frame of SMILES and 2691 structural fingerprints
#' @export
generate_fingerprints <- function(smiles) {
  

  fp_index_path <- system.file( # find path to fpIndex
    "extdata",
    "fpIndex_v2.0.csv",
    package = "fishFingers"
  )

  fpIndex <- read.csv(fp_index_path, check.names = FALSE)

  # get mols from SMILES
  mols <- parse.smiles(smiles) # parse SMILES to S4 object
  if (is.null(mols) || length(mols) == 0) {
    stop("Failed to parse SMILES")
  }

  # OpenBabel fingerprints
  openbabel <- as.character(fpIndex$smarts[fpIndex$fpType == "openbabel"]) # read SMARTS for custom fp # nolint
  colnames1 <- as.character(fpIndex$fpName[fpIndex$fpType == "openbabel"]) # read names for custom fp # nolint
  custom.fp1 <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = openbabel) # calculate fp
  custom.fp1 <- fp.to.matrix(custom.fp1) # convert fp list to matrix
  colnames(custom.fp1) <- colnames1 # apply fp names to columns

  # CDK (substructure) fingerprints
  substructure <- lapply(mols, get.fingerprint, type = "substructure")
  substructure <- fp.to.matrix(substructure)
  colnames(substructure) <- paste0("Un", seq(55, 54 + ncol(substructure)))
  substructure <- substructure[,fpIndex$fpName[fpIndex$fpType == "substructure"], drop = FALSE]

  # MACCS
  maccs <- lapply(mols, get.fingerprint, type = "maccs")
  maccs <- fp.to.matrix(maccs)
  colnames(maccs) <- paste0("Un", seq(362, 361 + ncol(maccs)))
  maccs <- maccs[,fpIndex$fpName[fpIndex$fpType == "maccs"], drop = FALSE]

  # PubChem CACTVS
  cactvs <- lapply(mols, get.fingerprint, type = "pubchem")
  cactvs <- fp.to.matrix(cactvs)
  colnames(cactvs) <- paste0("Un", seq(528, 527 + ncol(cactvs)))
  cactvs <- cactvs[,fpIndex$fpName[fpIndex$fpType == "cactvs"], drop = FALSE]

  # Klekota-Roth
  kroth <- lapply(mols, get.fingerprint, type = "kr")
  kroth <- fp.to.matrix(kroth)
  colnames(kroth) <- paste0("Un", seq(1409, 1408 + ncol(kroth)))
  kroth <- kroth[,fpIndex$fpName[fpIndex$fpType == "kroth"], drop = FALSE]

  # BioSMARTS
  biosmarts <- as.character(fpIndex$smarts[fpIndex$fpType == "biosmarts"])
  colnames2 <- as.character(fpIndex$fpName[fpIndex$fpType == "biosmarts"])
  custom.fp2 <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = biosmarts)
  custom.fp2 <- fp.to.matrix(custom.fp2)
  colnames(custom.fp2) <- colnames2

  # Ring
  ring.smarts <- as.character(fpIndex$smarts[fpIndex$fpType == "ring"])
  colnames3 <- as.character(fpIndex$fpName[fpIndex$fpType == "ring"])
  ring <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = ring.smarts)
  ring <- fp.to.matrix(ring)
  colnames(ring) <- colnames3

  # In Silico
  insilico.smarts <- as.character(fpIndex$smarts[fpIndex$fpType == "insilico"])
  colnames4 <- as.character(fpIndex$fpName[fpIndex$fpType == "insilico"])
  insilico <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = insilico.smarts)
  insilico <- fp.to.matrix(insilico)
  colnames(insilico) <- colnames4



  # build model data
 fp <- custom.fp1 %>%
    bind_cols(substructure) %>%
    bind_cols(maccs) %>%
    bind_cols(cactvs) %>%
    bind_cols(kroth) %>%
    bind_cols(custom.fp2) %>%
    bind_cols(ring) %>%
    bind_cols(insilico)

  return(fp)

}





#' Read SIRIUS v5.x fingerprints from a project directory
#'
#' This function imports posterior probabilities from files in the SIRIUS CSI:FingerID project directory.
#'
#' @description This function is deprecated and will be removed in a future version. Please use `read_sirius6_fingerprints()` instead for SIRIUS v6.x projects.
#' @details Deprecated in v0.2.0 - SIRIUS v5.x is no longer supported. Please use `read_sirius6_fingerprints()` for SIRIUS v6.x projects.
#' @param sirius_project_dir character; path to the SIRIUS project folder
#' @param topMost logical; if TRUE, then only import the top 1 ranked candidate from each result
#' @return a data frame with file, molecular_formula and posterior probability vectors (named by fingerprint column)
#' @keywords internal
#' @export
read_sirius5_fingerprints <- function(sirius_project_dir, topMost = TRUE) {

  tmp <- make_temp_dir()
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fp_index_path <- system.file( # find path to fpIndex
    "extdata",
    "fpIndex.csv",
    package = "fishFingers"
  )

  fpIndex <- read.csv(fp_index_path, check.names = FALSE)

  sirius_folders <- list.dirs(
    sirius_project_dir,
    full.names = TRUE,
    recursive = FALSE
  )

  result_list <- vector("list", length(sirius_folders))
  out_i <- 1L

  for (feature in sirius_folders) {

    feature_id <- basename(feature)

    fpZip <- list.files(feature, "fingerprint", full.names = TRUE)
    fpZip <- fpZip[!file.info(fpZip)$isdir]
    if (length(fpZip) == 0) {
      next
    }

    scoresZip <- list.files(feature, "scores", full.names = TRUE)
    scoresZip <- scoresZip[!file.info(scoresZip)$isdir]
    if (length(scoresZip) == 0) {
      next
    }

    feat_tmp <- file.path(tmp, feature_id)
    dir.create(feat_tmp, showWarnings = FALSE)

    unzip(scoresZip, exdir = feat_tmp, overwrite = TRUE)

    scoreFiles <- list.files(feat_tmp, pattern = "\\.info$", full.names = TRUE)
    scoreFiles <- scoreFiles[!grepl("compound\\.info$", scoreFiles)]

    if (length(scoreFiles) == 0) {
      next
    }

    scoreTable <- data.table::rbindlist(
      lapply(scoreFiles, function(file) {
        name <- tools::file_path_sans_ext(basename(file))
        dt <- data.table::fread(file, header = FALSE)
        dt <- dt[V1 == "sirius.scores.SiriusScore"]
        dt[, candidate := name]
        dt
      }),
      use.names = TRUE,
      fill = TRUE
    )

    if (topMost && nrow(scoreTable) == 0) {
      next
    }

    if (topMost) {
      top_name <- scoreTable$candidate[which.max(scoreTable$V2)]
      unzip(
        fpZip,
        files = paste0(top_name, ".fpt"),
        exdir = feat_tmp,
        overwrite = TRUE
      )
    } else {
      unzip(fpZip, exdir = feat_tmp, overwrite = TRUE)
    }

    fpFiles <- list.files(feat_tmp, pattern = "\\.fpt$", full.names = TRUE)
    if (length(fpFiles) == 0) {
      next
    }

    pos_files <- fpFiles[grepl("\\]\\+\\.fpt$", fpFiles)]
    neg_files <- fpFiles[grepl("\\]\\-\\.fpt$", fpFiles)]

    pos_dt <- if (length(pos_files)) {
      data.table::rbindlist(
        lapply(pos_files, read_fpt, fpIndex = fpIndex),
        use.names = TRUE
      )
    }

    neg_dt <- if (length(neg_files)) {
      data.table::rbindlist(
        lapply(neg_files, read_fpt, fpIndex = fpIndex),
        use.names = TRUE
      )
    }

    feature_dt <- data.table::rbindlist(
      list(pos_dt, neg_dt),
      use.names = TRUE,
      fill = TRUE
    )

    if (!is.null(feature_dt)) {
      feature_dt[, feature_id := feature_id]
      data.table::setcolorder(feature_dt, "feature_id")
      result_list[[out_i]] <- feature_dt
      out_i <- out_i + 1L
    }
  }

  if (out_i == 1L) {
    return(data.table::data.table())
  }

  result_dt <- data.table::rbindlist(result_list[seq_len(out_i - 1L)], use.names = TRUE, fill = TRUE)
  fp_names <- fpIndex$fpName[fpIndex$fpType != "ecfp6"]

  # Add missing columns explicitly
  missing_cols <- setdiff(fp_names, names(result_dt))

  for (col in missing_cols) {
    result_dt[, (col) := 0L]
  }

  # Now enforce full column set and order
  result_dt <- result_dt[, c("feature_id", "compound", fp_names), with = FALSE]

  # Replace NA with 0
  for (col in fp_names) {
    set(result_dt, which(is.na(result_dt[[col]])), col, 0L)
  }

  return(result_dt)
}




#' Read SIRIUS fingerprints from active REST API
#'
#' This function imports posterior probabilities from HTTPS API started by SIRIUS v6.x. The API must be active and accessible at the specified URL.
#'
#' @param path character; path of .sirius project file or directory (if NULL, will attempt to connect to any open SIRIUS project)
#' @param topMost logical; if TRUE, then only import the top 1 ranked candidate from each result
#' @return a data frame with file, molecular_formula and posterior probability vectors (named by fingerprint column)
#' @export
read_sirius_fingerprints <- function(path = NULL, topMost = TRUE, charge = 1) {

  # Create a local SDK object each call.
  sdk <- SiriusSDK$new()
  api <- NULL

  # Try to attach to an existing SIRIUS instance first.
  api <- tryCatch(
    sdk$attach_to_sirius(sirius_major_version = 6),
    error = function(e) NULL
  )

  if (!is.null(api)) {
    healthy <- tryCatch(
      api$actuator_api$Health()$status == "UP",
      error = function(e) FALSE
    )
    if (!healthy) {
      api <- NULL
    }
  }

  # If not attached or unhealthy, start or attach-or-start.
  if (is.null(api)) {
    api <- tryCatch(
      sdk$attach_or_start_sirius(),
      error = function(e) NULL
    )
    if (is.null(api)) {
      stop("Unable to connect to SIRIUS. Ensure SIRIUS is running or available in PATH; set SIRIUS_EXE if needed.")
    }
  }

  project_id <- sirius_init(path, api)

  # Get aligned features from the API
  features <- api$features_api$GetAlignedFeatures(project_id)
  cat("Found", length(features), "aligned features\n")
  if (length(features) == 0) {
    stop()
  }

  # Get the reference fp index from the project
  # Its unclear if this changes depending on the presence/absence of fp in the data set
  fingerid_data <- api$projects_api$GetFingerIdData(project_id, charge = charge)
  fingerid_data$fpName <- paste0("Un", fingerid_data$absoluteIndex)

  if (topMost == "formula") {
    cat("Importing fingerprints for top-ranked formulas only\n")
    top_fp <- extract_fingerprints(features, api, fingerid_data, project_id, topMost = "formula")
    return(top_fp)
} else if (topMost == "compound") {
    cat("Importing fingerprints for top-ranked compounds only\n")
    top_fp <- extract_fingerprints(features, api, fingerid_data, project_id, topMost = "compound")
    return(top_fp)
} else {
    cat("Importing fingerprints for all candidates\n")
    all_fp <- extract_fingerprints(features, api, fingerid_data, project_id, topMost = FALSE)
    return(all_fp)
    
  }
}