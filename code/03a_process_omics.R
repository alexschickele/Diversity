
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME <- "DIVERSITY_PAPER"
HILL_NB <- seq(0,5,0.25)
# CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_METAGENOMIC_2025-05-16 16:17:58.023044"
# CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_METAGENOMIC_RAREFIED_2025-06-16 16:33:14.774913"
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524"

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

# --- 2. Build ensembles
# --- 2.1. Extract file information
# Which subfolder list and which model in the ensemble
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]
ensemble_files <- mclapply(model_files, function(x){
  memory_cleanup() # low memory use
  
  load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/MODEL.RData"))
  if(length(MODEL$MODEL_LIST) >= 1){
    load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/QUERY.RData"))
    return(list(SUBFOLDER_NAME = x, MODEL_LIST = MODEL$MODEL_LIST, Y = QUERY$Y, MESS = QUERY$MESS))
  } else {
    return(NULL)
  } # if model list
}, mc.cores = MAX_CLUSTERS) %>% .[lengths(.) != 0]

 # --- 2.2. Loop over the files
 message(paste0(Sys.time(), "--- METAGENOMICS: build the ensembles - loop over files"))
 tmp <- mclapply(ensemble_files, function(x){
  memory_cleanup() # low memory use
  
  # --- 2.2.1. Load MODEL files
  load(paste0(CEPHALOPOD_OUTPUT,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 2.2.2. Extract projections in a matrix
  # If there is more than 1 algorithm, we extract and average across algorithm
  # Output matrix is cell x bootstrap x month
  if(length(x$MODEL_LIST) > 1){
    m <- lapply(x$MODEL_LIST, function(y){
      MODEL[[y]][["proj"]]$y_hat
    }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
  } else {
    m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
  } # end if
}, mc.cores = MAX_CLUSTERS, mc.cleanup = TRUE)

# --- 2.3. Stack in a cell x species x bootstrap x month matrix
# --- 2.3.1. Re-arrange the array
message(paste0(Sys.time(), "--- METAGENOMICS: build the ensembles - format to array"))
data <- tmp %>% 
  abind(along = 4) %>% 
  aperm(c(1,4,2,3))
# --- 2.3.2. Pretty dimensions
dimnames(data)[[2]] <- lapply(ensemble_files, function(x){out <- x$SUBFOLDER_NAME}) %>% unlist() %>% as.character()
dimnames(data)[[4]] <- 1:12 %>% as.character()

# --- 2.4. Memory cleanup
rm(tmp)
gc() # clean garbage and temporary files
message(paste0(Sys.time(), "--- METAGENOMICS: build the ensembles - DONE"))

# --- 3. Save
save(data, file = paste0("./output/", FOLDER_NAME, "/METAGENOMICS_raw_diversity.RData"))
message(paste0(Sys.time(), "--- METAGENOMICS: build the ensembles - DONE"))






# END