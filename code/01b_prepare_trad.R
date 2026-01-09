
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME = "DIVERSITY_PAPER"
HILL_NB <- seq(0,5,0.25)

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

# --- 1.4. List files
ATLANTECO_filenames <- list.files("/nfs/kryo/work/public/shared/AtlantECO/BASE/v1_legacy/", full.names = TRUE) %>% 
  .[grepl("biomass_", .)] %>% 
  .[grepl("_traditional_", .)] %>% 
  .[grepl(".csv", .)]

ATLANTECO_filenames <- ATLANTECO_filenames[-5] # remove old copepoda version

ATLANTECO_shortnames <- lapply(1:length(ATLANTECO_filenames), function(x)(x = strsplit(gsub("^[^_]*_[^_]*_[^_]*_(.*)_.+", "\\1", ATLANTECO_filenames), "_")[[x]][1])) %>% 
  unlist()

# --- 1.5. List other parameters
DATA_SOURCE <- c("abundance","biomass")

# --- 1.6. Concatenate list
PARAMETER_TBL <- expand_grid(PFG = 1:length(ATLANTECO_shortnames),
                             DATA_SOURCE = 1:length(DATA_SOURCE))

message(paste(Sys.time(), "--- Initialization of the R environment - DONE"))
message(">>> Starting to process the datasets")

# --- 2. Extract files
# and concatenate from PFG-level to class-level

ATLANTECO_all_aggr <- mclapply(1:nrow(PARAMETER_TBL), function(x){
  
  # --- 2.1. Open file
  ATLANTECO_file <- vroom(ATLANTECO_filenames[PARAMETER_TBL$PFG[x]])
  colnames(ATLANTECO_file) <- tolower(colnames(ATLANTECO_file)) # harmonize to lower
  
  if(ATLANTECO_filenames[PARAMETER_TBL$PFG[x]] == "/nfs/kryo/work/public/shared/AtlantECO/BASE/v1_legagy/AtlantECO-BASE-v1_microbiome_traditional_Euphausiacea_abund+biomass_20221220.csv"){
    # Repair euphausid data set
    # EcoMon and SSRF cruises present systematic abundance outliers. We remove the cruises. 
    # cf the reference cruise papaer, the abundance were given in ind/X m3 where X is the total water filter by the net, not ind/m3
    id <- grep("EcoMon|SSRF", ATLANTECO_file$origcollectionid)
    ATLANTECO_file <- ATLANTECO_file[-id,]
  } # euphausid repair
  
  # --- 2.2. Depth repair
  # If depth is empty, replace by min depth for the filtering
  id <- which(is.na(ATLANTECO_file$depth))
  ATLANTECO_file$depth[id] <- ATLANTECO_file$mindepth[id]
  
  # --- 2.3. Do some filtering
  # Depth; year; zero's; CPR
  ATLANTECO_file <- ATLANTECO_file %>% 
    .[!grepl("CPR|Recorder|Richardson|270", .$samplingprotocol),] %>% # remove CPR
    dplyr::filter(measurementvalue > 0) %>%  # remove zero
    dplyr::filter(depth <= 50) %>% 
    dplyr::filter(year > 1950)
  
  # --- 2.4. Get column names to select + select + rename
  names_qc <- list(names_abundance <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","measurementvalue","measurementunit", "species","class"),
                   names_biomass <- c("scientificname","worms_id","decimallatitude","decimallongitude","depth","year","month","meanbiomass","biomass_mgcm3","biomassunit","species","class"))
  
  ATLANTECO_df <- ATLANTECO_file %>% 
    dplyr::select(any_of(names_qc[[PARAMETER_TBL$DATA_SOURCE[x]]])) 
  
  # --- 2.5. Name repair if no unit associated with biomass_mgcm3
  if(colnames(ATLANTECO_df[8]) == "biomass_mgcm3" | colnames(ATLANTECO_df[8]) == "meanbiomass" & colnames(ATLANTECO_df[9]) == "species"){
    message(paste("NAME REPAIR --- Added biomassunit column for", ATLANTECO_shortnames[PARAMETER_TBL$PFG[x]]))
    ATLANTECO_df <- ATLANTECO_df %>% 
      mutate(biomassunit = "mgC m-3") %>% 
      dplyr::select(any_of(names_qc[[PARAMETER_TBL$DATA_SOURCE[x]]])) 
  } # end if
  
  colnames(ATLANTECO_df) <- names_qc[[1]] # rename to the standard measurementvalue + measurementunit
  
  # --- 2.6. Unit repair for abundance
  if(PARAMETER_TBL$DATA_SOURCE[x] == 1){
    # --- 2.6.1. Repair with no consequences on the measurementvalue
    ATLANTECO_df$measurementunit <- gsub("#", "ind", ATLANTECO_df$measurementunit) # replace
    ATLANTECO_df$measurementunit <- gsub("cell", "ind", ATLANTECO_df$measurementunit) # replace
    ATLANTECO_df$measurementunit <- gsub("/", " ", ATLANTECO_df$measurementunit) # replace
    ATLANTECO_df$measurementunit <- gsub("m3", "m-3", ATLANTECO_df$measurementunit) # replace
    
    # --- 2.6.1. Repair with consequences on the measurementvalue
    # L to m-3
    id <- grep("L", ATLANTECO_df$measurementunit)
    ATLANTECO_df$measurementunit <- gsub("L", "m-3", ATLANTECO_df$measurementunit) # replace
    ATLANTECO_df$measurementvalue[id] <- ATLANTECO_df$measurementvalue[id]*1e-3 # convert
    
    # cm3 to m-3
    id <- grep("cm-3", ATLANTECO_df$measurementunit)
    ATLANTECO_df$measurementunit <- gsub("cm-3", "m-3", ATLANTECO_df$measurementunit) # replace
    ATLANTECO_df$measurementvalue[id] <- ATLANTECO_df$measurementvalue[id]*1e-6 # convert
    
  } # if abundance
  
  # --- 2.7. Reduce spatial resolution
  res <- 1 # resolution
  digit <- nchar(sub('^0+','',sub('\\.','',res)))-1
  
  ATLANTECO_df <- ATLANTECO_df %>% 
    mutate(decimallatitude = round(decimallatitude+0.5*res, digits = digit)-0.5*res) %>%
    mutate(decimallongitude = round(decimallongitude+0.5*res, digits = digit)-0.5*res)
  
  # --- 2.8. Median by worms ID in a grid x time
  # Median is more representative and less affected by outliers
  ATLANTECO_taxa_avg <- ATLANTECO_df %>% 
    group_by(worms_id, decimallatitude, decimallongitude, month) %>% 
    mutate(measurementvalue = median(measurementvalue, na.rm = TRUE)) %>% 
    ungroup() %>% 
    distinct()
  
  # --- 2.9. Aggregate at the species level
  ATLANTECO_aggr <- ATLANTECO_taxa_avg %>% 
    dplyr::group_by(decimallatitude, decimallongitude, month, measurementunit, species, class) %>% 
    summarize(measurementvalue = sum(measurementvalue, na.rm = TRUE), 
              depth = mean(depth),
              year = mean(year)) %>% 
    dplyr::filter(!is.na(species))
  
  # --- 2.10. Return
  return(ATLANTECO_aggr)
  
}, mc.cores = 10) %>% bind_rows()

# --- 2.10. Secure zero filtering
# For some reason, I need a second check as the parallel function skips some zero's
ATLANTECO_all_aggr <- ATLANTECO_all_aggr %>% 
  dplyr::filter(measurementvalue > 0)

message(paste(Sys.time(), "--- Aggregate at the class level - DONE"))
message(">>> Starting to process the diversity")

# --- 3. Get top taxa worldwide
# --- 3.1. Define sample_id and build contingency table
sample_taxa <- ATLANTECO_all_aggr %>% 
  dplyr::group_by(decimallongitude, decimallatitude, month, measurementunit) %>% 
  mutate(site_id = cur_group_id()) %>% 
  ungroup() %>% 
  group_by(site_id, measurementunit, class) %>% 
  mutate(n_species = n_distinct(species), .groups = "drop") %>% 
  dplyr::select(site_id, measurementunit, class, measurementvalue, n_species)

contingency_taxa <- expand_grid(sample_taxa$site_id %>% unique(),
                                ATLANTECO_all_aggr$class %>% unique(),
                                ATLANTECO_all_aggr$measurementunit %>% unique())
colnames(contingency_taxa) <- c("site_id","class","measurementunit")

all_taxa <- contingency_taxa %>% left_join(sample_taxa)
# We count abscence as a 0 here, as it would be done in the diversity estimate
# This is to have a more representative contribution of each taxa class to our diversity estimates
all_taxa$measurementvalue[is.na(all_taxa$measurementvalue)] <- 0

# --- 3.2. For abundance
abundance_taxa_comp <- all_taxa %>% 
  dplyr::filter(measurementunit == "ind m-3") %>% 
  group_by(site_id, class) %>% 
  summarise(measurementvalue0 = sum(measurementvalue), n_species = mean(n_species, na.rm = T)) %>% # sum all species by class
  group_by(class) %>% 
  summarise(measurementvalue = mean(measurementvalue0), n = sum(measurementvalue0 != 0), n_species = mean(n_species, na.rm = T)) %>% # mean across class samples
  ungroup()

write.csv(abundance_taxa_comp, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_ABUNDANCE_taxo.csv"), row.names = F)

# --- 3.2. For biomass
biomass_taxa_comp <- all_taxa %>% 
  dplyr::filter(measurementunit == "mgC m-3") %>% 
  group_by(site_id, class) %>% 
  summarise(measurementvalue0 = sum(measurementvalue), n_species = mean(n_species, na.rm = T)) %>% # sum all species by class
  group_by(class) %>% 
  summarise(measurementvalue = mean(measurementvalue0), n = sum(measurementvalue0 != 0), n_species = mean(n_species, na.rm = T)) %>% # mean across class samples
  ungroup()

write.csv(biomass_taxa_comp, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_BIOMASS_taxo.csv"), row.names = F)

# --- 3. Compute diversity indices
# Loop over hill numbers

ATLANTECO_all_diversities <- lapply(1:length(HILL_NB), function(x){
  
  # --- 3.1. ID by site x measurementunit
  ATLANTECO_diversity <- ATLANTECO_all_aggr %>% 
    dplyr::group_by(decimallongitude, decimallatitude, month, measurementunit) %>% 
    mutate(site_id = cur_group_id()) %>% 
    ungroup()
  
  ID <- ATLANTECO_diversity$site_id %>% unique()
  
  # --- 3.2. Compute diversity
  ATLANTECO_diversity <- mclapply(1:length(ID), function(z){
    tmp <- ATLANTECO_diversity %>% 
      dplyr::filter(site_id == z) %>% 
      dplyr::select(species, measurementvalue) %>% 
      pivot_wider(names_from = species, values_from = measurementvalue, values_fn = mean)
    
    hill_value <- hill_taxa(comm = tmp, q = HILL_NB[x])
    
    out <- ATLANTECO_diversity %>% 
      dplyr::filter(site_id == z) %>% 
      dplyr::select(-species) %>% 
      mutate(measurementvalue = hill_value,
             scientificname = paste("Hill", HILL_NB[x], "(", measurementunit, ")"),
             worms_id = paste("Hill", HILL_NB[x], "(", measurementunit, ")"),
             taxonrank = "Species") %>% 
      distinct()
    
    return(out)
  }, mc.cores = 10) %>% bind_rows() %>% 
    dplyr::select(-site_id)
  
}) %>% bind_rows() %>% distinct()

message(paste(Sys.time(), "--- Processing diversity - DONE"))

# --- 4. Intermediate save
if(!file.exists(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME))){
  dir.create(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME))  # Create the folder if it does not exist
}
write.csv(ATLANTECO_all_diversities, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_raw_input.csv"), row.names = F)
message("--- Intermediate save done; now normalizing the data --- \n")

# --- 5. Quantile rescale
# We care about modelling the pattern, between 0 and 100
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

target_all_diversities_norm <- ATLANTECO_all_diversities %>% 
  group_by(scientificname) %>% 
  mutate(measurementvalue = quantile_scale(measurementvalue)) %>% 
  ungroup()

# --- 6. Save normalized
write.csv(target_all_diversities_norm, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_norm_input.csv"), row.names = F)
message("--- Normalized save done; now normalizing the data --- \n")

