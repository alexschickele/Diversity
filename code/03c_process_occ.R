
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME <- "DIVERSITY_PAPER"
HILL_NB <- seq(0,5,0.25)
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_OCCURRENCE_2025-06-12 16:07:51.160496"

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

# --- 1.4. Fake raster for later
r0 <- terra::rast(nrows = 180, ncols = 360, xmin = -180, xmax = 180, ymin = -90, ymax = 90)

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
message(paste0(Sys.time(), "--- OCCURRENCE: build the ensembles - loop over files"))
tmp <- mclapply(ensemble_files, function(x){
  memory_cleanup() # low memory use
  
  # --- 2.2.1. Load MODEL files
  load(paste0(CEPHALOPOD_OUTPUT,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
  
  # --- 2.2.2. Extract projections in a matrix
  # If there is more than 1 algorithm, we extract and average across algorithm
  # Output matrix is cell x bootstrap x month
  if(!is.null(MODEL[[1]][["proj"]]$y_hat)){
    if(length(x$MODEL_LIST) >= 1){
      m <- lapply(x$MODEL_LIST, function(y){
        MODEL[[y]][["proj"]]$y_hat
      }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
    } else {
      m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
    } # if size of model list
  } else {
    m <- NULL
  } # security if there is any projections
  
    # --- 2.2.3. Extract VIP
  # If there is more than 1 algorithm, we extract from the computed ensemble
  # Else we extract from the algorithm element
  if(length(x$MODEL_LIST) > 1){
    vip <- MODEL$ENSEMBLE$vip
  } else {
    vip <- MODEL[[x$MODEL_LIST]][["vip"]]
  } # end if
  
  # --- 2.2.4. Return
  return(list(m = m, vip = vip))
}, mc.cores = 30, mc.cleanup = TRUE) %>% .[lengths(.) != 0]

# --- 2.3. Stack in a cell x species x bootstrap x month matrix
# --- 2.3.1. Re-arrange the array
message(paste0(Sys.time(), "--- OCCURRENCE: build the ensembles - format to array"))
data <- lapply(tmp, function(x)(x = x[[1]])) %>% 
  abind(along = 4) %>% 
  aperm(c(1,4,2,3))

# --- 2.3.2. Pretty dimensions
# Does not matter much as we are going to compute diversity across dimensions
dimnames(data)[[4]] <- 1:12 %>% as.character()
message(paste0(Sys.time(), "--- OCCURRENCE: build the ensembles - done"))

# --- 2.4. Stack the VIP
vip <- lapply(tmp, function(x)(x = x[[2]])) %>% 
  bind_rows() %>% 
  group_by(variable) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  mutate(value = value / sum(value))

# --- 3. Save
save(vip, file = paste0("./output/", FOLDER_NAME, "/OCCURRENCE_vip.RData"))
save(data, file = paste0("./output/", FOLDER_NAME, "/OCCURRENCE_projections.RData"))
message(paste0(Sys.time(), "--- OCCURRENCE: build the ensembles - DONE"))

# --- 3. Initialize diversity computing
full_cell_id <- which(!is.na(data[,1,1,1]))
cell_vector <- terra::values(r0) %>% as.numeric()
loop_over <- expand.grid(HILL_NB, dimnames(data)[[3]], dimnames(data)[[4]])
names(loop_over) <- c("Hill_value","Bootstrap","Month")

# --- 4. Alpha diversities
# --- 4.1. Compute the diversities
message(paste0(Sys.time(), "--- OCCURRENCE: compute alpha diversity - START"))
alpha_div_list <- mclapply(1:nrow(loop_over), function(x){
  memory_cleanup() # low memory use
  
  # --- 4.1.1. Subset full cell x bootstrap x month
  df0 <- data[full_cell_id,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
  
  # --- 4.1.2. Compute the diversity
  div0 <- hill_taxa(comm = df0, q = loop_over$Hill_value[x]) # compute diversity
  
  # --- 6.1.3. Assign cells back
  div <- cell_vector # base
  div[full_cell_id] <- div0 # fill non empty cells
  return(div)
  
}, mc.cores = 30) # limited due to memory use
message(paste0(Sys.time(), "--- OCCURRENCE: compute alpha diversity - DONE"))

# --- 4.2. Back transformation to array
alpha_div_data <- array(data = NA,
                        dim = c(length(cell_vector), length(HILL_NB), dim(data)[[3]], dim(data)[[4]]),
                        dimnames = list(NULL,
                                        as.character(HILL_NB),
                                        dimnames(data)[[3]], 
                                        dimnames(data)[[4]]))

for(x in 1:nrow(loop_over)){
  alpha_div_data[, as.character(loop_over$Hill_value[x]), loop_over$Bootstrap[x], loop_over$Month[x]] <- alpha_div_list[[x]]
} # end for

alpha_div_data[is.infinite(alpha_div_data)] <- NA # security
message(paste0(Sys.time(), "--- OCCURRENCE: format alpha diversity - DONE"))

# --- 5. Save
data <- alpha_div_data # same name as the others
save(data, file = paste0("./output/", FOLDER_NAME, "/OCCURRENCE_raw_diversity.RData"))
message(paste0(Sys.time(), "--- OCCURRENCE: save - DONE"))

