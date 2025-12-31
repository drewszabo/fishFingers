# R/utils.R
make_temp_dir <- function() {
  tmp <- tempfile("fishFingers_")
  dir.create(tmp)
  tmp
}

check_model <- function(model_path = system.file("extdata", "xgb_model.rds", package = "fishFingers")) {
  if (!file.exists(model_path)) stop("xgb_model.rds not found at: ", model_path, ". Please reinstall fishFingers from Github.")
  model <- readRDS(model_path)
  model
}

read_fpt <- function(fp, fpIndex) {

  dt <- data.table::fread(fp, header = FALSE)
  dt <- data.table::transpose(dt)

  polarity <- if (grepl("\\]\\+\\.fpt$", fp)) {
    "pos"
  } else if (grepl("\\]\\-\\.fpt$", fp)) {
    "neg"
  } else {
    stop("Cannot determine polarity from filename: ", fp)
  }

  fp_names <- if (polarity == "pos") {
    fpIndex$fpName[fpIndex$pos]
  } else {
    fpIndex$fpName[fpIndex$neg]
  }

  if (ncol(dt) != length(fp_names)) {
    stop(
      "Fingerprint length mismatch for ", basename(fp),
      ": expected ", length(fp_names),
      ", got ", ncol(dt)
    )
  }

  name <- tools::file_path_sans_ext(basename(fp))
  data.table::setnames(dt, fp_names)
  dt[, compound := name]
  data.table::setcolorder(dt, "compound")
  dt
}


make_species_matrix <- function(species) {

  if (!check_species(species)) {
    stop(
      "Please check spelling of species and try again.",
      call. = FALSE
    )
  }

  species <- tolower(trimws(species))
  species <- str_replace(species, " ", "_")

  species_factor <- readRDS(
    system.file("extdata", "species.rds", package = "fishFingers")
  )

  species_levels <- tolower(levels(species_factor))

  sp_vec <- as.integer(species_levels == species)

  matrix(
    sp_vec,
    nrow = 1,
    dimnames = list(NULL, levels(species_factor))
  )
}
