
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME = "DIVERSITY_PAPER"
HILL_NB <- seq(0,5,0.25)

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

# --- 1.4. Define prokaryotic size fractions to consider
size_to_consider <- c("0.22-3", "0.22-1.6", "0.2<-") # Size fractions relevant to prokaryotes from different cruises

# --- 2. Extract base data
# --- 2.1. Retrieve and filter sample metadata
samples <- vroom("/net/meso/work/aschickele/Nitromics/data/samples-metadata.tsv") %>% 
  dplyr::select(Sample, dataset, station, latitude, longitude, date, size_fraction, depth_num) %>% 
  dplyr::filter(depth_num < 50) %>%  # Select samples from 0 - 50 m layer
  dplyr::filter(size_fraction %in% size_to_consider) %>%  # Filter based on defined size fractions
  collect()

# Standardize column names to lowercase with underscores
names(samples) <- tolower(gsub(" ", "_", names(samples)))

# --- 2.2. Retrieve mOTU abundance data
# And extract the mOTU reference from taxonomy
profile <- vroom("/net/meso/work/aschickele/Diversity/data/OMD_v1/motus-profiles.tsv", skip = 2) %>% 
  dplyr::select(c(`#consensus_taxonomy`, any_of(unique(samples$sample))))
colnames(profile)[1] <- "consensus_taxonomy"
profile <- profile %>% 
  mutate(mOTU_taxonomy = str_trim(str_extract(consensus_taxonomy, "^[^\\[]+")),
         mOTUs_Species_Cluster = str_extract(consensus_taxonomy, "(?<=\\[)[^\\]]+(?=\\])")) %>% 
  dplyr::select(-consensus_taxonomy) %>% 
  dplyr::select(mOTU_taxonomy, mOTUs_Species_Cluster, everything())
colnames(profile)[2] <- "mOTUs Species Cluster"

# --- 2.3. Rarefy the profiles
# Get the number of positive samples per location - we chose the threshold accordingly
sample_size <- apply(profile[, -c(1:2)], 2, function(x)(x = length(x[x > 0]))) %>% quantile(seq(0,1,0.01))
sample_size <- 1000

# Sample of equal sample sizes
id <- sample(1:nrow(profile), size = sample_size)
profile <- profile[id, ]

# --- 2.5. Reshape to long format
# And remove zero's
data <- profile %>% 
  as.data.frame() %>% 
  pivot_longer(cols = 3:ncol(.),  # Transform wide table to long format
               names_to = "sample",
               values_to = "reads") %>% 
  dplyr::filter(reads > 0)

# --- 2.6. Merge sample metadata with KO abundance data
data <- inner_join(samples, data)

# --- 2.4. Extract taxonomy
# --- 2.4.1. Base table
taxonomy <- vroom("/net/meso/work/aschickele/Nitromics/data/genomes-summary.csv") %>% 
  dplyr::select(`mOTUs Species Cluster`, `GTDB Taxonomy`) %>%
  collect()

# --- 2.4.2. Function to split a taxonomic annotation into separate columns, handling NA values
split_taxonomy <- function(taxon) {
  ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")  # Expected ranks
  
  if (is.na(taxon)) return(setNames(rep(NA, length(ranks)), ranks))  # Handle NA input
  
  # Split taxonomic annotation by semicolon, then extract rank and name using "__"
  taxon_split <- strsplit(taxon, ";")[[1]]
  taxon_split <- sapply(taxon_split, function(x) strsplit(x, "__")[[1]], simplify = FALSE)
  taxonomy_vector <- setNames(rep(NA, length(ranks)), ranks)  # Initialize output with NA
  
  # Map rank codes to full names and extract corresponding taxonomic names
  rank_map <- c("d"="Domain", "p"="Phylum", "c"="Class", "o"="Order", 
                "f"="Family", "g"="Genus", "s"="Species")
  
  for (item in taxon_split) {
    if (length(item) == 2 && item[1] %in% names(rank_map)) {
      taxonomy_vector[rank_map[item[1]]] <- item[2] } # end if
  } # end for
  
  return(taxonomy_vector)
} # end function

# --- 2.4.3. Apply the function to the dataset
taxonomy_matrix <- t(apply(taxonomy[,2], 1, split_taxonomy))
taxonomy <- cbind(taxonomy, taxonomy_matrix) %>% dplyr::select(`mOTUs Species Cluster`, Species, Class) %>% 
  distinct()

# --- 2.4.4. Join
data <- inner_join(data, taxonomy)

# --- 2.4. Sum to 1 per sample
data <- data %>% 
  group_by(sample) %>% 
  mutate(reads = reads / sum(reads)) %>% 
  ungroup()

# --- 2.5. Remove 0
data <- data %>% 
  dplyr::filter(reads > 0)

rm(annotations, profile, samples)
gc()

# --- 2.6. Fix colnames
colnames(data) <- tolower(colnames(data))

# --- 3. Get top taxa worldwide
omic_taxa_comp <- data %>% 
  group_by(station, class) %>% 
  summarize(reads = sum(reads)) %>% 
  group_by(class) %>% 
  summarise(reads = median(reads), n = n()) %>% 
  ungroup()

write.csv(omic_taxa_comp, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/METAGENOMIC_taxo.csv"), row.names = F)

# --- 3. Match CEPHALOPOD input requirements
# --- 3.1. Extract the input table
cephalopod_input <- data %>% 
  dplyr::select(species, latitude, longitude, date, depth_num, reads, size_fraction) %>% 
  mutate(date = as.POSIXct(date, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),  # Convert to date-time
         month = as.integer(format(date, "%m")), # extract numeric month
         year = as.integer(format(date, "%Y"))) %>% # extract numeric year
  dplyr::select(-date) 

# --- 3.1.2. Sum reads per Genome x sample
# if multiple genomes of the sample
cephalopod_input <- cephalopod_input %>% 
  group_by(species, latitude, longitude, year, month, depth_num, size_fraction) %>% 
  summarise(sum_reads = sum(reads)) %>% 
  ungroup()

# --- 3.1.3. Match column names
names(cephalopod_input) <- c("scientificname","decimallatitude","decimallongitude","year","month","depth","size_fraction","measurementvalue")
cephalopod_input <- cephalopod_input %>% 
  mutate(worms_id = scientificname,
         measurementunit = "relative metagenomic reads",
         taxonrank = "Species") %>% 
  filter(!is.na(scientificname))
cephalopod_input$depth <- as.numeric(cephalopod_input$depth)

# --- 4. Assemble first hill
target_all_diversities <- lapply(1:length(HILL_NB), function(x){
  
  # --- 4.1. ID by site x measurementunit
  target_diversity <- cephalopod_input %>% 
    dplyr::group_by(decimallongitude, decimallatitude, month, measurementunit, size_fraction) %>% 
    mutate(site_id = cur_group_id()) %>% 
    ungroup()
  
  ID <- target_diversity$site_id %>% unique()
  
  # --- 4.2. Compute diversity
  target_diversity <- mclapply(1:length(ID), function(z){
    tmp <- target_diversity %>% 
      dplyr::filter(site_id == z) %>% 
      dplyr::select(scientificname, measurementvalue) %>% 
      pivot_wider(names_from = scientificname, values_from = measurementvalue, values_fn = mean)
    
    hill_value <- hill_taxa(comm = tmp, q = HILL_NB[x])

    out <-target_diversity %>% 
      dplyr::filter(site_id == z) %>% 
      mutate(measurementvalue = hill_value,
             scientificname = paste("Hill", HILL_NB[x], "(", measurementunit, ")"),
             worms_id = paste("Hill", HILL_NB[x], "(", measurementunit, ")"),
             taxonrank = "Species")
    
    return(out)
  }, mc.cores = 10) %>% bind_rows() %>% 
    dplyr::select(-site_id)
  
}) %>% bind_rows() %>% distinct()

# --- 6. Intermediate save
if(!file.exists(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME))){
  dir.create(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME))  # Create the folder if it does not exist
}
write.csv(target_all_diversities, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/METAGENOMIC_raw_input.csv"), row.names = F)
message("--- Intermediate save done; now normalizing the data --- \n")

# --- 6. Quantile rescale
# We care about modelling the pattern, between 0 and 100
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

target_all_diversities_norm <- target_all_diversities %>% 
  group_by(scientificname) %>% 
  mutate(measurementvalue = quantile_scale(measurementvalue)) %>% 
  ungroup()

# --- 7. Save normalized
write.csv(target_all_diversities_norm, file = paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/METAGENOMIC_norm_input.csv"), row.names = F)
message("--- Normalized save done; now normalizing the data --- \n")

