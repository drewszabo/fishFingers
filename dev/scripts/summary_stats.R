# Summary of top species

top_species <- df %>%
  select(1:6) %>%
  group_by(species) %>%
  distinct(SMILES, .keep_all = TRUE) %>%
  summarise(n = dplyr::n(), .groups = "drop")  

smiles_multi_species <- df %>%
  distinct(SMILES, species) %>%      # remove repeated SMILES–species combos
  count(SMILES, name = "n_species") %>%  
  filter(n_species > 1)

smiles_multi_species %>%
  filter(str_detect(SMILES, "C")) %>%      # keep only SMILES with carbon
  summarise(
    n_Cl = sum(str_detect(SMILES, "Cl")),
    n_F  = sum(str_detect(SMILES, "F")),
    n_Br = sum(str_detect(SMILES, "Br"))
  )
