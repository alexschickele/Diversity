
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name
FOLDER_NAME = "DIVERSITY_PAPER"

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")
library(worrms)

# --- 2. Open PHYTOBASE
# --- 2.1. Build input data
# Filter the required columns
phyto <- vroom("/net/sea/work/public/shared/AtlantECO/BASE/v1_legacy/AtlantECO-BASE-v1_microbiome_traditional_phytoplankton_species_occurrences_PhytoBasev2_20220905.csv")
colnames(phyto) <- tolower(colnames(phyto))
phyto_in <- phyto %>% 
  .[!grepl("CPR|Recorder|Richardson|270", .$samplingprotocol),] %>% # remove CPR
  dplyr::filter(measurementvalue != 0) %>%  # remove zero
  dplyr::filter(depth <= 200) %>% 
  dplyr::filter(year > 1900) %>% 
  dplyr::select(c("scientificname", "worms_id", "decimallatitude", "decimallongitude", 
                  "depth", "year", "month", "measurementvalue", "measurementunit", "taxonrank")) %>% 
  dplyr::filter(measurementvalue == "Presence" & taxonrank == "Species" & depth < 50) %>% 
  mutate(worms_id = as.character(worms_id))

# --- 2.2. Get taxonomy as well
phyto_taxo_comp <- phyto_in %>% 
  group_by(worms_id, scientificname) %>% 
  summarize(n = n())

tmp <- wm_classification_(name = phyto_taxo_comp$scientificname) %>% 
  dplyr::filter(rank == "Class") %>% 
  dplyr::select(id, scientificname)
colnames(tmp) <- c("worms_id","class")

phyto_taxo_comp <- phyto_taxo_comp %>% 
  left_join(tmp) %>% 
  group_by(class) %>% 
  summarise(n = sum(n)) %>% 
  ungroup()

# --- 3. Open ZOOBASE
# --- 3.1. Build input data
# Filter the required columns
zoo <- vroom("/net/sea/work/public/shared/AtlantECO/BASE/v1_legacy/AtlantECO-BASE-v1_microbiome_traditional_zooplankton_species_occurrences_ZooBasev2_20220909.csv")
colnames(zoo) <- tolower(colnames(zoo))
zoo_in <- zoo %>%
  .[!grepl("CPR|Recorder|Richardson|270", .$samplingprotocol),] %>% # remove CPR
  dplyr::filter(measurementvalue != 0) %>%  # remove zero
  dplyr::filter(depth <= 50) %>% 
  dplyr::filter(year > 1900) %>% 
  dplyr::select(c("scientificname", "worms_id", "decimallatitude", "decimallongitude", 
                  "depth", "year", "month", "measurementvalue", "measurementunit", "taxonrank")) %>% 
  dplyr::filter(measurementvalue == "Presence" & taxonrank == "Species" & depth < 50) %>% 
  mutate(worms_id = as.character(worms_id))

# --- 3.2. Get taxonomy as well
zoo_taxo_comp <- zoo_in %>% 
  group_by(worms_id, scientificname) %>% 
  summarize(n = n())

tmp <- wm_classification_(name = zoo_taxo_comp$scientificname) %>% 
  dplyr::filter(rank == "Class") %>% 
  dplyr::select(id, scientificname)
colnames(tmp) <- c("worms_id","class")

zoo_taxo_comp <- zoo_taxo_comp %>% 
  left_join(tmp) %>% 
  group_by(class) %>% 
  summarise(n = sum(n)) %>% 
  ungroup()

# --- 4. Assemble & save
occ <- bind_rows(phyto_in, zoo_in) # all species occurrence above 50m
write.csv(occ, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/OCCURRENCE_raw_input.csv"), row.names = F)
write.csv(occ, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/OCCURRENCE_norm_input.csv"), row.names = F)

occ_taxo_comp <- bind_rows(phyto_taxo_comp, zoo_taxo_comp) %>% 
  dplyr::filter(!is.na(class))

write.csv(occ_taxo_comp, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/OCCURRENCE_taxo.csv"), row.names = F)

message("--- Normalized save done; now normalizing the data --- \n")



# END