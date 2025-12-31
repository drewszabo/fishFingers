library(tidyverse)
library(caret)
library(RColorBrewer)

# all_species_levels is a factor/vector of 112 species
chem_species_df <- expand.grid(
  SMILES = pull(distinct(df, SMILES)),
  species = all_species_levels
)

chem_info <- df %>%
  select(3:6) %>%
  distinct(SMILES, .keep_all = TRUE)

# Merge fingerprints onto the expanded grid
full_features <- df %>%
  distinct(SMILES, .keep_all = TRUE) %>%
  select(-c(1,2,3,5,6))

chem_species_df <- chem_species_df %>%
  left_join(full_features, by = "SMILES")

# Create dummy variables for species
species_matrix <- model.matrix(~ species - 1, data = chem_species_df)  # removes intercept
colnames(species_matrix) <- sub("^species", "", colnames(species_matrix))

# Combine fingerprints and species dummies
# Assuming fingerprint columns are from column 3 onward
full_x <- cbind(chem_species_df[, -(1:2)], species_matrix)

# Suppose your model is 'bcf_model'
predicted_log_bcf <- predict(xgb_model, newdata = full_x)

# Add predictions back to the data frame
pubchem_species_df <- chem_species_df %>%
  select(SMILES, species) %>%
  left_join(chem_info, by = "SMILES")
pubchem_species_df$log_bcf_pred <- as.numeric(predicted_log_bcf)

# BCF Scarcity Data

pubchem_species_df <- pubchem_species_df %>%
  mutate(bcf_level = case_when(log_bcf_pred < log10(2000) ~ "nB",
                               log_bcf_pred > log10(5000) ~ "vB",
                               TRUE ~ "B"))

pubchem_species_df <- pubchem_species_df %>%
  mutate(species = str_replace(species, "_", " "),
         species_short = str_replace(species, 
                                     "^([A-Za-z])[a-z]*\\s+",  # match genus
                                     "\\1. "))

# Compute mean predicted log BCF per species and per chemical
species_order <- pubchem_species_df %>%
  group_by(species_short) %>%
  summarise(mean_bcf = mean(log_bcf_pred, na.rm = TRUE)) %>%
  arrange(desc(mean_bcf))

chemical_order <- pubchem_species_df %>%
  group_by(PubChem_CID) %>%
  summarise(mean_bcf = mean(log_bcf_pred, na.rm = TRUE)) %>%
  arrange(desc(mean_bcf))

pubchem_species_df <- pubchem_species_df %>%
  mutate(
    species_short = factor(species_short, levels = species_order$species_short),
    PubChem_CID = factor(PubChem_CID, levels = chemical_order$PubChem_CID)
  )

ggplot(data = pubchem_species_df) +
  geom_tile(aes(x = species_short, y = as.factor(PubChem_CID), fill = predicted_log_bcf)) +
  scale_x_discrete(name = "Species") +
  scale_fill_gradientn(name = "Predicted \nlog10BCF",
                       colours = rev(brewer.pal(7, "RdBu"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6, face = "italic"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 6),
        legend.position = c(0.5,0.95),
        legend.direction = "horizontal",
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.box.spacing = unit(0, "cm"),
        legend.background = element_rect(fill = "transparent"),
        legend.margin = margin(0,0,0,0, "cm"),
        legend.key.height = unit(0.2, "cm"),
        plot.margin=grid::unit(c(1,0.05,0.05,0.05), "cm"),
        panel.spacing = unit(0, "mm")
  )

ggsave("data_predicted_bcf.png",
       width = 7,
       height = 5,
       units = "in",
       dpi = 320)
