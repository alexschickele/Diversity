# --- 1. System arguments

# --- 2. R Packages
# --- 2.1. General use
if(!require("devtools")){install.packages("devtools")}
if(!require("abind")){install.packages("abind")}

# --- 2.2. Tidy environment-related
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("parallel")){install.packages("parallel")}

# --- 2.4. Spatial data and object
if(!require("raster")){install.packages("raster")}
if(!require("virtualspecies")){install.packages("virtualspecies")}

# --- 2.6. Others
if(!require("RColorBrewer")){install.packages("RColorBrewer")}
if(!require("scales")){install.packages("scales")}

# --- Seed
set.seed(123)

# --- Necessary code steps

# --- Custom functions
source("./function/viridis.R")
source("./function/get_cell_neighbors.R")
source("./function/memory_cleanup.R")

# --- Model specific parameters
MAX_CLUSTERS <- 10
