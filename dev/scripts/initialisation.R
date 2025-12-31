library(tidyverse)
library(caret)

# Parameters
input_path <- "bcfModelData.csv"
seed <- 42
train_prop <- 0.80
test_prop <- 0.10
valid_prop <- 0.10

set.seed(seed)

# --------------- Load data ----------------------------------------------------------------
df <- read_csv(input_path, show_col_types = FALSE, na = c("", "NA", "N/R"))

# remove NA species
df <- filter(df, !is.na(species))
# clean species name
df <- df %>%
  mutate(species = gsub(" ", "_", df$species),
         species = str_trim(species, side = "both"))

massbank <- read.csv("massbank_filtered.csv", check.names = FALSE) %>%
  filter(grepl("ESI", Instrument_type)) # only soft ionisation

# --------------- Basic checks ----------------------------------------------------------------
# canonical column names handling
colnames(df) <- make.names(colnames(df))

# Identify obvious columns
smiles_col <- "SMILES"

species_col <- "species"

target_col <- "log_bcf"

# Fingerprint columns: beginning with "Un"
fp_cols <- colnames(df)[grep("Un", colnames(df))]

# Check binary nature of a subset (fast)
check_n <- length(fp_cols)
if (check_n > 0) {
  samp_fp <- fp_cols[1:check_n]
  is_binary <- sapply(samp_fp, function(c) all(df[[c]] %in% c(0,1)))
  cat(sum(is_binary), "/", check_n, "checked fingerprint columns are binary (0/1).\n")
} else {
  cat("One or more fingerprint columns are not in binary format.\n")
}

# Check for missing values
missing_total <- sum(is.na(df[,fp_cols]))
cat("Total missing cells in dataset:", missing_total, "\n")
missing_per_col <- colSums(is.na(df))
if (any(missing_per_col>0)) print(sort(missing_per_col[missing_per_col>0], decreasing = TRUE))

# Basic target distribution
cat("Target (", target_col, ") summary:\n")
print(summary(df[[target_col]]))

# Per-species counts
cat("Top species counts:\n")
print(df %>% count(!!sym(species_col)) %>% arrange(desc(n)) %>% head(20))

# Unique SMILES
n_smiles <- df[[smiles_col]] %>% unique() %>% length()
unique_smiles_vec <- unique(df$SMILES)
cat("Unique SMILES:", n_smiles, "Total observations:", nrow(df), "\n")

# --------------- EDA visuals (optional) ----------------------------------------------------
# We'll create a few simple plots and save them to PNG files. The user can open them later.

# 1) Target histogram + density
png("target_hist_density.png", width = 800, height = 600)
qplot(df[[target_col]], geom = c("histogram", "density"), bins = 30, main = paste("Distribution of", target_col), xlab = target_col)
dev.off()

# 2) Boxplot of target by top species (top 15 species)
top_species <- df %>% count(!!sym(species_col)) %>% arrange(desc(n)) %>% head(5) %>% pull(!!sym(species_col))
png("target_by_species_boxplot.png", width = 1000, height = 600)
  df %>% filter(species %in% top_species) %>%
    ggplot(aes(x = species, y = log_bcf)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Target by top", length(top_species), "species"))
dev.off()


# 3) Fingerprint sparsity summary (fraction of ones per fingerprint) - show top 20 most/least frequent
if (length(fp_cols)>0) {
  fp_sums <- colSums(df[fp_cols]==1)
  fp_freq <- sort(fp_sums / nrow(df), decreasing = TRUE)
  fp_freq_df <- tibble(feature = names(fp_freq), freq = as.numeric(fp_freq))
  write_csv(fp_freq_df, "fp_frequency_summary.csv")
  cat("Fingerprint frequency summary saved to fp_frequency_summary.csv\n")
}

# --------------- SMILES-level 80/10/10 split with MassBank-prioritized validation -----------

# Identify SMILES that are present in MassBank
massbank_smiles <- unique(all_sirius$SMILES) # perform SIRIUS calculations and validation first to identify valid features
massbank_inchikey <- unique(all_sirius$InChIKey)
valid_candidates <- unique(df$SMILES[df$SMILES %in% massbank_smiles | df$InChIKey %in% massbank_inchikey])

n_valid_smiles <- floor(valid_prop * n_smiles)
n_train_smiles <- floor(train_prop * n_smiles)
n_test_smiles  <- n_smiles - n_valid_smiles - n_train_smiles



# Select validation SMILES (prioritized)
if (length(valid_candidates) < n_valid_smiles) {
  stop("Not enough MassBank matches to fill validation set.")
}
valid_smiles <- sample(valid_candidates, n_valid_smiles)

# Remaining SMILES (excluding validation)
remaining_smiles <- setdiff(unique_smiles_vec, valid_smiles)

# Shuffle remaining for train/test
shuffled_remaining <- sample(remaining_smiles, length(remaining_smiles))

# Split remaining into train and test
train_smiles <- shuffled_remaining[1:n_train_smiles]
test_smiles  <- shuffled_remaining[(n_train_smiles + 1):(n_train_smiles + n_test_smiles)]

cat("SMILES counts -> train:", length(train_smiles), 
    " test:", length(test_smiles), 
    " valid:", length(valid_smiles), "\n")

train_df <- df %>% filter(!!sym(smiles_col) %in% train_smiles)
test_df  <- df %>% filter(!!sym(smiles_col) %in% test_smiles)
valid_df <- df %>% filter(!!sym(smiles_col) %in% valid_smiles)

write.csv(train_df, "train_df.csv", row.names = FALSE)
write.csv(test_df, "test_df.csv", row.names = FALSE)
write.csv(valid_df, "valid_df.csv", row.names = FALSE)

cat("Observation counts -> train:", nrow(train_df), 
    " test:", nrow(test_df), 
    " valid:", nrow(valid_df), "\n")

# ---------------- Preprocess -------------------------------------------------
# Ensure species is a factor with levels set across ALL data
all_species_levels <- unique(df$species)
train_df$species <- factor(train_df$species, levels = all_species_levels)
test_df$species  <- factor(test_df$species, levels = all_species_levels)
valid_df$species <- factor(valid_df$species, levels = all_species_levels)

# Now one-hot encode consistently
options(na.action='na.pass')
train_species_mat <- model.matrix(~ species - 1, data = train_df)
colnames(train_species_mat) <- sub("^species", "", colnames(train_species_mat))
test_species_mat  <- model.matrix(~ species - 1, data = test_df)
colnames(test_species_mat) <- sub("^species", "", colnames(test_species_mat))
valid_species_mat <- model.matrix(~ species - 1, data = valid_df)
colnames(valid_species_mat) <- sub("^species", "", colnames(valid_species_mat))

# Combine
train_x <- cbind(as.data.frame(select(train_df, fp_cols)), as.data.frame(train_species_mat))
test_x  <- cbind(as.data.frame(select(test_df, fp_cols)),  as.data.frame(test_species_mat))
valid_x <- cbind(as.data.frame(select(valid_df, fp_cols)), as.data.frame(valid_species_mat))

# Targets
train_y <- train_df[[target_col]]
test_y  <- test_df[[target_col]]
valid_y <- valid_df[[target_col]]


