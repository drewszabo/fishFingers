#' SIRIUS Validation
#' 

library(tidyverse)

# Import fpIndex
fpIndex <- read_delim("fpIndex.csv", delim = "\t")

# Locate all score files and extract
scoreZip <- list.files("sirius/", pattern = "scores", recursive = TRUE, full.names = TRUE)

# Unzip scores compressed files
for(scores in scoreZip){
  dir <- dirname(scores)
  unzip(scores, exdir = dir, overwrite = TRUE)
}

# Determine annotation with best score
siriusFolders <- list.dirs("sirius/", full.names = TRUE, recursive = FALSE)

for (folder in siriusFolders) {
  scoreFiles <- list.files(folder, pattern = "\\.info$", full.names = TRUE)
  scoreFiles <- scoreFiles[!grepl(scoreFiles, pattern = "compound.info")]
  
  i = 1
  for(file in scoreFiles){
    name <- tools::file_path_sans_ext(basename(file))
    score <- as_tibble(read.delim(file, header = FALSE)) %>%
      filter(V1 == "sirius.scores.SiriusScore") %>%
      mutate(V1 = name)
    
    if(i == 1){
      scoreTable <- score
    } else {
      scoreTable <- bind_rows(scoreTable, score)
    }
    i = i + 1
  }
    # get max score and extract fp file
    scoreTable <- scoreTable %>%
      arrange(desc(V2)) %>%
      slice_head()
    
    # extract fp from top ranked sirius feature
    fpFile <- list.files(folder, pattern = "fingerprint", full.names = TRUE)
    fpFile <- fpFile[file_test("-f", fpFile)]
    
    # skip folders without fingerprints file
    if (length(fpFile) == 0) {
      message("No fingerprint zip found in folder: ", folder)
      next
    }
    
    unzip(zipfile = fpFile, files = paste0(scoreTable$V1[1], ".fpt"), exdir = folder, overwrite = TRUE)
  
}

# Locate all fingerprint zip files and extract
fpFiles <- list.files("sirius/", pattern = ".fpt", recursive = TRUE, full.names = TRUE)


i = 1
j = 1
for(fp in fpFiles) {
  
  # positive mode
  if(grepl("\\]\\+\\.fpt$", fp)) {
    
    # read in .fpt file
    dp <- read.delim(fp, header = FALSE)
    dp <- as.tibble(t(dp)) %>%
      mutate(path = fp, .before = 1)
    
    if(i == 1) {
      all_pos <- dp
    } else {
      all_pos <- bind_rows(all_pos, dp)
    }
    
    i = i + 1
  }
  
  # negative mode
  if(grepl("\\]\\-\\.fpt$", fp)) {
    
    # read in .fpt file
    dp <- read.delim(fp, header = FALSE)
    dp <- as.tibble(t(dp)) %>%
      mutate(path = fp, .before = 1)
    
    if(j == 1) {
      all_neg <- dp
    } else {
      all_neg <- bind_rows(all_neg, dp)
    }
    
    j = j + 1
  }
  
}

# get fingerprint names for pos and neg
fpName_pos <- fpIndex %>%
  filter(pos == TRUE) %>%
  select(fpName)

fpName_neg <- fpIndex %>%
  filter(neg == TRUE) %>%
  select(fpName)

# rename columns to match
colnames(all_pos) <- c("path", fpName_pos$fpName)
colnames(all_neg) <- c("path", fpName_neg$fpName)

all_pos <- all_pos %>%
  mutate(InChIKey = str_split_i(path, "_", i = 3),
         Lab = str_split_i(path, "_", i = 4),
         Inst = str_split_i(str_extract(path, "(LC-|GC-|ESI-|APCI-)[^/]+"), "_", 1),
         .after = 1)

all_neg <- all_neg %>%
  mutate(InChIKey = str_split_i(path, "_", i = 3),
         Lab = str_split_i(path, "_", i = 4),
         Inst = str_split_i(str_extract(path, "(LC-|GC-|ESI-|APCI-)[^/]+"), "_", 1),
         .after = 1)

# add missing fingerprints P = 0
for (f in fp_cols) {
  if (!f %in% colnames(all_pos)) {
    all_pos[[f]] <- 0
  }
}
all_pos <- all_pos[, c("InChIKey", "Lab", "Inst", fp_cols), drop = FALSE]

# add missing fingerprints P = 0
for (f in fp_cols) {
  if (!f %in% colnames(all_neg)) {
    all_neg[[f]] <- 0
  }
}
all_neg <- all_neg[, c("InChIKey", "Lab", "Inst", fp_cols), drop = FALSE]

# join pos and neg SIRIUS fp prob
all_sirius <- rbind(all_pos, all_neg)

# sort then select unique InChIKey (QFT>IFIT>QTOF LC>GC>Infusion)
write.csv(all_sirius, "all_sirius.csv", row.names = FALSE)
all_sirius <- read.csv("all_sirius.csv", check.names = FALSE)

# join with SMILES
filtered_massbank_table <- patRoon::as.data.table(filtered_massbank)
all_sirius <- filtered_massbank_table %>%
  select(SMILES, InChIKey) %>%
  distinct(InChIKey, .keep_all = TRUE) %>%
  left_join(all_sirius, by = "InChIKey") %>%
  drop_na(Lab) %>%
  select(-Lab, -Inst)

# amend fp in valid_df
valid_df <- valid_df %>%
  select(-c(fp_cols, SMILES)) %>%
  left_join(valid_x_mc, by = "InChIKey")

# Evaluate model
evaluate_model(xgb_model, valid_x, valid_y, cv = FALSE)

# Monte Carlo simulation (parallelized) for SIRIUS+CSI:FingerID posterior probabilities

# Load required libraries
library(future.apply)

# ----- PARAMETERS -----
N <- 10000  # number of Monte Carlo samples

# ----- FUNCTIONS -----

generate_mc_matrix <- function(prob_matrix) {
  prob_matrix <- as.matrix(prob_matrix)
  
  nr <- nrow(prob_matrix)
  nc <- ncol(prob_matrix)
  
  mc_bin <- matrix(
    rbinom(nr * nc, size = 1, prob = prob_matrix),
    nrow = nr,
    ncol = nc,
    dimnames = dimnames(prob_matrix)
  )
  
  return(mc_bin)
}

generate_best_mc_binary <- function(prob_matrix, N = 10000, threshold = 0.5) {
  
  prob_matrix <- as.matrix(prob_matrix)
  nr <- nrow(prob_matrix)
  nc <- ncol(prob_matrix)
  
  # Storage for running average
  avg_matrix <- matrix(0, nrow = nr, ncol = nc)
  
  for (i in seq_len(N)) {
    avg_matrix <- avg_matrix + generate_mc_matrix(prob_matrix)
  }
  
  # Convert counts to probabilities
  avg_matrix <- avg_matrix / N
  
  # Threshold → final deterministic binary matrix
  final_binary <- (avg_matrix >= threshold) * 1
  
  # Keep same row/col names
  dimnames(final_binary) <- dimnames(prob_matrix)
  
  return(final_binary)
}

valid_x_mc <- generate_best_mc_binary(all_sirius[,-c(1:2)], N = 10000)
valid_x_mc <- bind_cols(select(all_sirius, 1:2), valid_x_mc)


# Predict and plot validation data
valid_pred  <- predict(xgb_model, valid_x)

evaluate_model(xgb_model, valid_x, valid_y)

plot_valid <- data.frame(measured = valid_y, predicted = valid_pred) %>%
  bind_cols(select(valid_df, species))

# Create the new column
plot_valid <- plot_valid %>%
  mutate(top_species = if_else(species %in% top3, species, "Other"))

# annotate measured and predicted values
plot_valid <- plot_valid %>%
  mutate(measured_ann = case_when(measured < log10(2000) ~ "nB",
                                  measured >= log10(5000) ~ "vB",
                                  TRUE ~ "B"),
         predicted_ann = case_when(predicted < log10(2000) ~ "nB",
                                   predicted >= log10(5000) ~ "vB",
                                   TRUE ~ "B"))

# Convert to ordinal factors (numerical for ranking)
plot_valid$measured_ann <- factor(plot_valid$measured_ann,
                                 levels = c("nB", "B", "vB"),
                                 ordered = TRUE)

plot_valid$predicted_ann <- factor(plot_valid$predicted_ann,
                                  levels = c("nB", "B", "vB"),
                                  ordered = TRUE)

# Convert to numeric ranks
true  <- as.numeric(plot_valid$measured_ann)
pred  <- as.numeric(plot_valid$predicted_ann)

# Calculate confustion matrix
cm <- caret::confusionMatrix(plot_valid$predicted_ann,
                             plot_valid$measured_ann)

cm

# Calculate the quadratic weighted kappa (QWK)
qwk <- irr::kappa2(
  data.frame(true, pred),
  weight = "squared"
)

qwk

p3 <- ggplot(data = plot_valid, aes(x = measured, y = predicted, colour = factor(top_species, levels = c("Cyprinus_carpio", "Oncorhynchus_mykiss", "Poecilia_reticulata", "Other")))) +
  geom_point(size = 1, shape = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -2, linetype = "dotted", color = "grey40") +
  scale_x_continuous(name = "Measured logBCF",
                     limits = c(-2, 6)) +
  scale_y_continuous(name = "Predicted logBCF",
                     limits = c(-1, 6)) +
  scale_color_discrete(name = "",
                       labels = c(expression(italic("C. carpio")), expression(italic("O. mykiss")), expression(italic("P. reticulata")), "Other")) +
  theme_bw() +
  theme(
    legend.position = "none",
    # legend.background = element_rect(fill = "transparent")
    ) +
  labs(title = "Validation") +
  guides(colour = guide_legend(ncol=2)) +
  coord_equal(ratio = 1)

p3

ggsave("valid_plot.svg",
       plot = p3,
       dpi = 1200,
       width = 4,
       height = 4,
       unit = "in")

# USEPA TEST Suite

epaTest_df <- read.csv("TEST/MyToxicity/Bioconcentration_factor_Consensus.csv", na.strings = "N/A") %>%
  janitor::clean_names()

epaTest_y_pred <- epaTest_df %>%
  select(pred_value_log10)%>%
  bind_cols(select(valid_df, log_bcf)) %>%
  drop_na(pred_value_log10)

rmse(as.numeric(epaTest_y_pred$log_bcf), as.numeric(epaTest_y_pred$pred_value_log10))
mae(as.numeric(epaTest_y_pred$log_bcf), as.numeric(epaTest_y_pred$pred_value_log10))
cor(as.numeric(epaTest_y_pred$log_bcf), as.numeric(epaTest_y_pred$pred_value_log10))^2

# annotate measured and predicted values
epaTest_y_pred <- epaTest_y_pred %>%
  mutate(measured_ann = case_when(log_bcf < log10(2000) ~ "nB",
                                  log_bcf >= log10(5000) ~ "vB",
                                  TRUE ~ "B"),
         predicted_ann = case_when(pred_value_log10 < log10(2000) ~ "nB",
                                   pred_value_log10 >= log10(5000) ~ "vB",
                                   TRUE ~ "B"))

# Convert to ordinal factors (numerical for ranking)
epaTest_y_pred$measured_ann <- factor(epaTest_y_pred$measured_ann,
                                  levels = c("nB", "B", "vB"),
                                  ordered = TRUE)

epaTest_y_pred$predicted_ann <- factor(epaTest_y_pred$predicted_ann,
                                   levels = c("nB", "B", "vB"),
                                   ordered = TRUE)

# Convert to numeric ranks
true  <- as.numeric(epaTest_y_pred$measured_ann)
pred  <- as.numeric(epaTest_y_pred$predicted_ann)

# Calculate confustion matrix
cm <- caret::confusionMatrix(epaTest_y_pred$predicted_ann,
                             epaTest_y_pred$measured_ann)

cm

# Calculate the quadratic weighted kappa (QWK)
qwk <- irr::kappa2(
  data.frame(true, pred),
  weight = "squared"
)

qwk

p4 <- ggplot(data = epaTest_y_pred, aes(x = log_bcf, y = pred_value_log10)) +
  geom_point(size = 1, shape = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -2, linetype = "dotted", color = "grey40") +
  scale_x_continuous(name = "Measured logBCF",
                     limits = c(-2, 6)) +
  scale_y_continuous(name = "Predicted logBCF",
                     limits = c(-1, 6)) +
  theme_bw() +
  coord_equal(ratio = 1)

p4

ggsave("epatest_plot.svg",
       plot = p4,
       dpi = 1200,
       width = 3.3,
       unit = "in")

# LogKow
kow_df <- valid_df %>%
  select(-fp_cols)

rmse(as.numeric(kow_df$log_bcf), as.numeric(kow_df$XLogP))
mae(as.numeric(kow_df$log_bcf), as.numeric(kow_df$XLogP))
cor(as.numeric(kow_df$log_bcf), as.numeric(kow_df$XLogP))^2

# annotate measured and predicted values
kow_df <- kow_df %>%
  mutate(measured_ann = case_when(log_bcf < log10(2000) ~ "nB",
                                  log_bcf >= log10(2000) ~ "B",
                                  TRUE ~ "B"),
         predicted_ann = case_when(XLogP < 4.5 ~ "nB",
                                   XLogP >= 4.5 ~ "B",
                                   TRUE ~ "B"))

# Convert to ordinal factors (numerical for ranking)
kow_df$measured_ann <- factor(kow_df$measured_ann,
                                      levels = c("nB", "B", "vB"),
                                      ordered = TRUE)

kow_df$predicted_ann <- factor(kow_df$predicted_ann,
                                       levels = c("nB", "B", "vB"),
                                       ordered = TRUE)

# Convert to numeric ranks
true  <- as.numeric(kow_df$measured_ann)
pred  <- as.numeric(kow_df$predicted_ann)

# Calculate confustion matrix
cm <- caret::confusionMatrix(kow_df$predicted_ann,
                             kow_df$measured_ann)

cm

# Calculate the quadratic weighted kappa (QWK)
qwk <- irr::kappa2(
  data.frame(true, pred),
  weight = "squared"
)

qwk

p5 <- ggplot(data = kow_df, aes(x = log_bcf, y = XLogP)) +
  geom_point(size = 1, shape = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -2, linetype = "dotted", color = "grey40") +
  scale_x_continuous(name = "Measured logBCF",
                     #limits = c(-2, 6)
                     ) +
  scale_y_continuous(name = "Predicted logKow",
                     #limits = c(-1, 6)
                     ) +
  theme_bw() +
  coord_equal(ratio = 1)

p5

ggsave("kow_plot.svg",
       plot = p5,
       dpi = 1200,
       width = 3.3,
       unit = "in")

