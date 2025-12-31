# MassBank database

massbank_msp <- patRoon::loadMSLibrary("MassBank_NIST.msp", "msp")

library(data.table)
library(fastcluster)

#  SMILES list
valid_smiles <- df$SMILES
valid_inchikey <- df$InChIKey

# Identify which records to keep
keep_idx <- (massbank_msp@records$SMILES %in% valid_smiles) |
  (massbank_msp@records$InChIKey %in% valid_inchikey)

# Subset both the records table and the corresponding spectra
filtered_massbank <- massbank_msp
filtered_massbank@records <- massbank_msp@records[keep_idx, ]
filtered_massbank@spectra <- massbank_msp@spectra[keep_idx]

# (optional) Check how many records are left
nrow(filtered_massbank@records)

makeConsensus <- function(massbank, bin_size = 0.01, mgf_file = "consensus.mgf") {
  # Extract metadata
  records <- as.data.table(massbank@records)
  records[, contributor := tstrsplit(DB_ID, "-", keep = 2)]
  records[, group_id := paste(InChIKey, contributor, Instrument_type, Ion_mode, Precursor_type, sep = "_")]
  
  # Flatten spectra into long table
  all_specs <- rbindlist(lapply(names(massbank@spectra), function(id) {
    dt <- as.data.table(massbank@spectra[[id]])
    dt[, DB_ID := id]
    dt
  }))
  
  specs <- merge(
    all_specs,
    records[, .(DB_ID, group_id, Name, PrecursorMZ, Precursor_type, Ion_mode, Collision_energy)],
    by = "DB_ID",
    all.x = TRUE
  )
  
  # Helper: bin spectrum -> returns data.table with columns bin (int), mz (centroid), intensity (sum)
  binSpectrum <- function(mz, intensity, bin_size) {
    bin_index <- as.integer(round(mz / bin_size))
    dt <- data.table(bin = bin_index, mz_orig = mz, intensity = intensity)
    dt <- dt[, .(
      mz = sum(mz_orig * intensity) / sum(intensity),   # intensity-weighted centroid within that spectrum
      intensity = sum(intensity)
    ), by = bin]
    return(dt)
  }
  
  # Process each group
  consensus_list <- lapply(split(specs, specs$group_id), function(g) {
    # binned: named list of data.tables for each original spectrum (names = DB_ID)
    binned <- lapply(split(g, g$DB_ID), function(sp) binSpectrum(sp$mz, sp$intensity, bin_size))
    names(binned) <- unique(g$DB_ID)
    
    # list of unique bin indices across spectra
    all_bins <- sort(unique(unlist(lapply(binned, function(x) x$bin))))
    if (length(all_bins) == 0) return(NULL)
    
    # build intensity matrix: rows = spectra, cols = bin indices
    mat <- do.call(rbind, lapply(binned, function(x) {
      row <- numeric(length(all_bins))
      row[match(x$bin, all_bins)] <- x$intensity
      row
    }))
    rownames(mat) <- names(binned)
    colnames(mat) <- as.character(all_bins)   # keep as strings for names
    
    # Cluster (or skip if single)
    if (nrow(mat) > 1) {
      hc <- hclust(dist(mat, method = "euclidean"), method = "average")
      cl <- cutree(hc, k = 1)  # one cluster -> consensus; can be extended later
      spec_vec <- colMeans(mat[cl == 1, , drop = FALSE])
    } else {
      spec_vec <- mat[1, ]
    }
    # spec_vec is a numeric vector named by bin indices (as strings)
    
    # Convert bin -> true centroid m/z using intensity-weighted mean across spectra
    consensus_bins <- as.integer(names(spec_vec))
    consensus_mz <- numeric(length(consensus_bins))
    for (i in seq_along(consensus_bins)) {
      b <- consensus_bins[i]
      mz_vals <- numeric(0)
      int_vals <- numeric(0)
      # collect centroid mz and intensity from each spectrum that had this bin
      for (sp in binned) {
        idx <- which(sp$bin == b)
        if (length(idx) > 0) {
          mz_vals <- c(mz_vals, sp$mz[idx])
          int_vals <- c(int_vals, sp$intensity[idx])
        }
      }
      if (length(mz_vals) == 0) {
        consensus_mz[i] <- NA_real_
      } else {
        consensus_mz[i] <- sum(mz_vals * int_vals) / sum(int_vals)   # weighted mean across spectra
      }
    }
    
    # Build consensus table and drop any NA mz (shouldn't usually happen)
    consensus_dt <- data.table(
      bin = consensus_bins,
      mz = consensus_mz,
      intensity = as.numeric(spec_vec),
      group_id = unique(g$group_id),
      precursor_mz = unique(g$PrecursorMZ),
      precursor_type = unique(g$Precursor_type),
      ion_mode = unique(g$Ion_mode),
      name = unique(g$Name),
      collision_energy = paste(unique(na.omit(g$Collision_energy)), collapse = ";")
    )
    consensus_dt <- consensus_dt[!is.na(mz) & intensity > 0]
    
    # order peaks by decreasing intensity (optional but common)
    setorder(consensus_dt, -intensity)
    
    return(consensus_dt)
  })
  
  # Combine all consensus
  consensus <- rbindlist(consensus_list, fill = TRUE)
  if (nrow(consensus) == 0) {
    message("No consensus spectra produced.")
    return(consensus)
  }
  
  # --- Export to MGF ---
  con <- file(mgf_file, "w")
  for (grp in unique(consensus$group_id)) {
    sub <- consensus[group_id == grp]
    if (nrow(sub) == 0) next
    # TITLE: include collision energy string if present (clean digits)
    ce_str <- if (nzchar(sub$collision_energy[1])) {
      paste0("_CE", gsub("[^0-9]", "", sub$collision_energy[1]))
    } else {
      ""
    }
    title_line <- paste0(sub$group_id[1], ce_str)
    writeLines("BEGIN IONS", con)
    writeLines(paste0("TITLE=", title_line), con)
    writeLines(paste0("PEPMASS=", sub$precursor_mz[1]), con)
    writeLines(paste0("CHARGE=", ifelse(toupper(sub$ion_mode[1]) == "POSITIVE", "1+", "1-")), con)
    writeLines(paste0("ION=", sub$precursor_type[1]), con)
    
    # write each peak: mz to 4 dp, intensity rounded to integer
    for (r in seq_len(nrow(sub))) {
      mz_val <- sprintf("%.4f", sub$mz[r])
      int_val <- sprintf("%.0f", sub$intensity[r])
      writeLines(paste(mz_val, int_val), con)
    }
    
    writeLines("END IONS\n", con)
  }
  close(con)
  
  message("Consensus spectra written to: ", mgf_file)
  return(consensus)
}

consensus <- makeConsensus(filtered_massbank, bin_size = 0.01, mgf_file = "consensus.mgf")
