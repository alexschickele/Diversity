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
run_name <- "PHYTOBIAS_PA_testing_wmess13"
input_folder <- paste0("/nfs/meso/work/aschickele/Bluecloud_WB_local/output/", run_name) # dynamic
output_folder <- paste0("/nfs/meso/work/aschickele/Diversity/output/", run_name) # dynamic

# --- 1. Building the diversity matrices
# Load function
source(file = "./code/01_diversity_assessment.R")

# Perform function
diversity_assessment(INPUT_FOLDER = input_folder,
                     OUTPUT_FOLDER = output_folder,
                     HILL_NB = 0:3,
                     BETA_DIV = FALSE,
                     FOCAL_SIZE = 2,
                     MAX_CLUSTERS = 10)








# --- END --- 