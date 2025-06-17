# --- 1. System arguments

# --- 2. R Packages
# --- 2.1. General use
if(!require("devtools")){install.packages("devtools")}
if(!require("abind")){install.packages("abind")}
if(!require("vroom")){install.packages("vroom")}

# --- 2.2. Tidy environment-related
if(!require("tidyverse")){install.packages("tidyverse")}
if(!require("parallel")){install.packages("parallel")}

# --- 2.3. Spatial data and object
if(!require("terra")){install.packages("raster")}
if(!require("virtualspecies")){install.packages("virtualspecies")}
if(!require("sf")){install.packages("sf")}
if(!require("rnaturalearth")){install.packages("rnaturalearth")}

# --- 2.4. Diversity
if(!require("hillR")){install.packages("hillR")}
if(!require("vegan")){install.packages("vegan")}

# --- 2.5. Others
if(!require("RColorBrewer")){install.packages("RColorBrewer")}
if(!require("scales")){install.packages("scales")}
if(!require("RSQLite")){install.packages("RSQLite")}
if(!require("ellipse")){install.packages("ellipse")}

# --- Seed
set.seed(123)

# --- Necessary code steps

# --- Custom functions
source("./function/viridis.R")
source("./function/get_cell_neighbors.R")
source("./function/memory_cleanup.R")
source("./function/r_to_hatched_polygon.R")

# --- Model specific parameters
MAX_CLUSTERS <- 10
