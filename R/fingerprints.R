# R/fingerprints.R
#'
#' Calculate fingerprints from SMILES
#'
#' @param smiles character; vector of SMILES
#' @return data.frame of SMILES and 2691 structural fingerprints
#' @export
generate_fingerprints <- function(smiles) {

  fp_index_path <- system.file( # find path to fpIndex
    "extdata",
    "fpIndex.csv",
    package = "fishFingers"
  )

  fpIndex <- read.csv(fp_index_path, check.names = FALSE)

  # get mols from SMILES
  mols <- parse.smiles(smiles) # parse SMILES to S4 object
  if (is.null(mols) || length(mols) == 0) {
    stop("Failed to parse SMILES")
  }

  # Un0 - Un54 fingerprints in SIRIUS custom
  smarts1 <- as.character(fpIndex$smarts[fpIndex$fpType == "custom1"]) # read SMARTS for custom fp
  colnames1 <- as.character(fpIndex$fpName[fpIndex$fpType == "custom1"]) # read names for custom fp
  custom.fp1 <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = smarts1) # calculate fp
  custom.fp1 <- fp.to.matrix(custom.fp1) # convert fp list to matrix
  colnames(custom.fp1) <- colnames1 # apply fp names to columns

  # CDK (substructure) fingerprints
  substructure <- lapply(mols, get.fingerprint, type = "substructure")
  substructure <- fp.to.matrix(substructure)
  colnames(substructure) <- paste0("Un", seq(55, 54 + ncol(substructure)))
  substructure_filter <- as.integer(fpIndex$typeIndex[fpIndex$fpType == "substructure"])
  substructure <- substructure[, substructure_filter, drop = FALSE]

  # MACCS
  maccs <- lapply(mols, get.fingerprint, type = "maccs")
  maccs <- fp.to.matrix(maccs)
  colnames(maccs) <- paste0("Un", seq(362, 361 + ncol(maccs)))
  maccs_filter <- as.integer(fpIndex$typeIndex[fpIndex$fpType == "maccs"])
  maccs <- maccs[, maccs_filter, drop = FALSE]

  # PubChem CACTVS
  cactvs <- lapply(mols, get.fingerprint, type = "pubchem")
  cactvs <- fp.to.matrix(cactvs)
  colnames(cactvs) <- paste0("Un", seq(528, 527 + ncol(cactvs)))
  cactvs_filter <- as.integer(fpIndex$typeIndex[fpIndex$fpType == "cactvs"])
  cactvs <- cactvs[, cactvs_filter, drop = FALSE]

  # Klekota-Roth
  kroth <- lapply(mols, get.fingerprint, type = "kr")
  kroth <- fp.to.matrix(kroth)
  colnames(kroth) <- paste0("Un", seq(1409, 1408 + ncol(kroth)))
  kroth_filter <- as.integer(fpIndex$typeIndex[fpIndex$fpType == "kroth"])
  kroth <- kroth[, kroth_filter, drop = FALSE]

  # Custom Fingerprints 2
  smarts2 <- as.character(fpIndex$smarts[fpIndex$fpType == "custom2"])
  colnames2 <- as.character(fpIndex$fpName[fpIndex$fpType == "custom2"])

  custom.fp2 <- lapply(mols, get.fingerprint, type = "substructure", substructure.pattern = smarts2)
  custom.fp2 <- fp.to.matrix(custom.fp2)
  colnames(custom.fp2) <- colnames2

  # build model data

  fp <- custom.fp1 %>%
    bind_cols(substructure) %>%
    bind_cols(maccs) %>%
    bind_cols(cactvs) %>%
    bind_cols(kroth) %>%
    bind_cols(custom.fp2)

  return(fp)

}





#' Read SIRIUS fingerprints from a project directory
#'
#' This function imports posterior probabilities from files in the SIRIUS CSI:FingerID project directory.
#'
#' @param sirius_project_dir character; path to the SIRIUS project folder
#' @param topMost logical; if TRUE, then only import the top 1 ranked candidate from each result
#' @return a data frame with file, molecular_formula and posterior probability vectors (named by fingerprint column)
#' @export
read_sirius_fingerprints <- function(sirius_project_dir, topMost = TRUE) {

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
    if (length(fpZip) == 0) {
      next
    }

    scoresZip <- list.files(feature, "scores", full.names = TRUE)
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
  fp_names_present <- intersect(fp_names, names(result_dt))
  result_dt <- result_dt[, c("feature_id", "compound", fp_names_present), with = FALSE]
  for (col in fp_names_present) {
    set(result_dt, which(is.na(result_dt[[col]])), col, 0L)
  }

  return(result_dt)
}


