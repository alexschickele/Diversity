#' ============================= MASTER SCRIPT =================================
#' For diversity assessment based on CEPHALOPOD outputs
#' A. Schickele 2023
#' =============================================================================

# --- 0. Start up and load functions
# --- Cleanup
rm(list=ls())
closeAllConnections()
setwd("/nfs/meso/work/aschickele/Diversity")
source(file = "./code/00_config.R")

# --- Directories
# run_name <- "PHYTOBASE_17092024"
run_name <- "ZOOBASE_17092024"
# run_name <- c("ZOOBASE_v2_pres_only_w_cpr","PHYTOBASE_v2_pres_only_w_cpr")
input_folder <- paste0("/nfs/meso/work/aschickele/CEPHALOPOD/output/", run_name) # dynamic
output_folder <- paste0("/nfs/meso/work/aschickele/Diversity/output/", run_name) # dynamic

# --- 1. Building the diversity matrices
# Load function
source(file = "./code/01_diversity_assessment.R")

# Perform function
diversity_assessment(INPUT_FOLDER = input_folder,
                     OUTPUT_FOLDER = output_folder,
                     HILL_NB = seq(0,5,0.25),
                     BETA_DIV = FALSE,
                     FOCAL_SIZE = 2,
                     MAX_CLUSTERS = 40)

# --- 2. Basic diversity plots
# Load function
source(file = "./code/02_diversity_plot.R")

# Perform function
diversity_plot(INPUT_FOLDER= input_folder,
               OUTPUT_FOLDER = output_folder,
               PROFILE = "climatology_P_SST_regridded")

# --- 3. Co-occurrence networks
# Load function






# --- END --- 