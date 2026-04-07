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


sirius_init <- function(path = NULL, api) {

  # Check if path is provided and valid
  if (!is.null(path)) {
    if (!dir.exists(path)) {
      stop("Provided path does not exist: ", path)
    } else if (length(list.files(path, pattern = "\\.sirius$")) == 0) {
      stop("Provided path does not appear to be a valid SIRIUS project directory: ", path)
    }
  }

  # Get and select project information
  if(is.null(path)) {
    projects <- api$projects_api$GetProjects()
  } else {
    projects <- api$projects_api$GetProject(path)
  }

  if (length(projects) == 0) {
    stop("No projects found on localhost.")
  } else if (length(projects) == 1) {
    project_id <- projects[[1]][["projectId"]]
    cat("Found 1 open project\n")
  } else {
    projects_df <- map_df(projects, function(p) {
      tibble(id  = p$projectId, loc = p$location)
    })
    cat("\nAvailable Projects:\n")
    for (i in seq_len(nrow(projects_df))) {
      cat(sprintf("[%d] %s (%s)\n", i, projects_df$id[i],
                  projects_df$loc[i]))
    }
    selection <- as.integer(readline(prompt = "\nSelect project number: "))
    if (is.na(selection) || selection < 1 ||
        selection > nrow(projects_df)) {
      stop("Invalid selection.")
    }
    project_id <- projects_df$id[selection]
    cat("Selected project:", project_id, "\n")
  }

  if (length(project_id) == 0) {
    stop("No project selected.")
  }

  project_id
}


extract_fingerprints <- function(features, api, fingerid_data, project_id, topMost = TRUE) {

  aligned_id <- as.character()
  for (i in seq_along(1:length(features))) {
    aligned_id <- append(aligned_id, features[[i]]$alignedFeatureId)
  }
  
  # Retrieve all formula candidates
  df <- map_df(aligned_id, function(fid) {
    # Get formula candidates for this aligned feature
    formulas <- api$features_api$GetFormulaCandidates(project_id, fid)
    
    if (length(formulas) == 0)
      return(NULL)
    
    # Convert formulas to a tibble
    map_df(formulas, function(f) {
      tibble(
        feature_id = features[[which(aligned_id == fid)]][["name"]],
        aligned_id = fid,
        formula_id = f$formulaId,
        mol_form   = f$molecularFormula,
        rank       = as.numeric(unlist(f$rank %||% NA_real_)),
        score_norm = as.numeric(unlist(f$siriusScoreNormalized %||% NA_real_)),
        score      = as.numeric(unlist(f$siriusScore %||% NA_real_))
      )
    })
  })
  
  if (topMost == TRUE) {
    # Select top‑ranked formula by rank, excluding NA ranks
    df <- df %>%
      group_by(aligned_id) %>%
      filter(rank == 1)
  }

  # Retrieve fingerprint prediction for selected formulas
  failed_count <- 0
  all_fps <- map_df(seq_len(nrow(df)), function(i) {
    # Skip samples without any calculated fps
    res <- tryCatch({
      api$features_api$GetFingerprintPrediction(project_id, df$aligned_id[i], df$formula_id[i])
    }, error = function(e) {
      # Increment failure count instead of warning immediately
      failed_count <<- failed_count + 1
      return(NULL)
    })
    
    # If API failed, skip this iteration
    if (is.null(res) || length(res) == 0) {
      failed_count <<- failed_count + 1
      return(NULL)
    }
    
    # Convert probabilities to characters
    probs <- as.character(res)
    
    # Make a single-row tibble (1 row, N columns)
    fp_df <- as_tibble(t(probs))
    names(fp_df) <- fingerid_data$fpName  # 5000 column names
    
    # Add identifiers
    fp_df %>%
      mutate(
        feature_id = df$feature_id[i],
        aligned_id = df$aligned_id[i],
        formula_id = df$formula_id[i],
        molecularFormula = df$mol_form[i]
      ) %>%
      select(feature_id, aligned_id, formula_id, molecularFormula, everything())
  })
  
  # Warn once with total failed count
  if (failed_count > 0) {
    warning(sprintf("%d feature(s) were omitted due to API failures or missing fingerprints.", failed_count))
  }
  
  return(all_fps)
  
}

# Function to calculate max Tanimoto similarity for each query against all training
calc_max_tanimoto <- function(query, training_matrix) {

  intersection <- as.numeric(training_matrix %*% query)
  ones_query <- sum(query)
  ones_ref <- rowSums(training_matrix)
  union <- ones_ref + ones_query - intersection
  
  tanimoto <- intersection / union
  max(tanimoto, na.rm = TRUE)
}
