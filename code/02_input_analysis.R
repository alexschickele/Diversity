
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME = "DIVERSITY_PAPER"
HILL_NB <- seq(0,5,0.25)

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

# --- 2. List files and open
# --- 2.1. List
norm_files <- list.files(paste0("./output/", FOLDER_NAME), full.names = T) %>% 
  .[grep("norm_input.csv", .)]

# --- 2.2. Open and stack
data <- lapply(norm_files, function(x){
  tmp <- vroom(x) %>% 
    mutate(measurementunit = x,
           worms_id = as.character(worms_id)) %>% 
    as.data.frame()
  if(is.character(tmp$measurementvalue)){
    tmp <- tmp %>% mutate(measurementvalue = 1)
  } # set to 1 if presence
  return(tmp)
}) %>% bind_rows()

# --- 3. Process plotting data
# --- 3.1. Robinson projection
robinson_proj <- "+proj=robin +lon_0=150 +datum=WGS84"

# --- 3.2. Project land mask
land <- ne_countries(scale = 50, returnclass = "sf") %>% .[,1]
land_raster <- terra::rasterize(land, terra::rast(res = 0.1))  # 0.1-degree resolution
land_rob <- terra::project(land_raster, robinson_proj)

# --- 3.3. Project samples
# Remove year and month as the projection is 2D here
tmp <- data %>% dplyr::select(-month, -year, -depth, -scientificname, -worms_id, - taxonrank, -measurementvalue) %>% distinct()
samples_sf <- sf::st_as_sf(tmp, coords = c("decimallongitude", "decimallatitude"), crs = 4326)
samples_proj <- sf::st_transform(samples_sf, crs = robinson_proj)

# --- 4. Point map
pal <- scales::alpha(c("chocolate2", "antiquewhite4", "darkolivegreen3"), 0.5)[as.factor(tmp$measurementunit)] # wes anderson

plot(land_rob, col = "gray20", box = FALSE, axes = FALSE)
points(sf::st_coordinates(samples_proj), col = pal, pch = 20, cex = 0.5) # add points
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, add = TRUE, lty = "dotted") # add grid

# --- 5. Profiles
data_type <- data$measurementunit %>% unique() # Get unique measurement units (e.g., types of measurements like Hill numbers)
par(mfrow = c(1,length(data_type))) # Set up plotting area: 1 row,N columns

# --- 5.1. Loop over each measurement unit
for(d in seq_along(data_type)){
  
  tmp <- data %>% filter(measurementunit == data_type[d]) # Filter data for the current measurement unit
  plot(tmp$measurementvalue, tmp$decimallatitude, pch = 15, col = scales::alpha("antiquewhite3", 0.01), xlim = c(0,100), ylim = c(-90,90), xlab = "", ylab = "") # Empty plot
  sc_names <- tmp$scientificname %>% unique() # Get unique scientific names (e.g., Hill numbers)
  pal <- viridis_pal(length(sc_names)) %>% rev() # Color palette
  
  # --- 5.2. Loop over each scientific name to draw smoothed lines
  lapply(seq_along(sc_names), function(x){
    
    # Filter data for current scientific name and sort by latitude
    df <- tmp %>% 
      filter(scientificname == sc_names[x]) %>% 
      arrange(decimallatitude)
    
    fit <- loess(measurementvalue ~ decimallatitude, data = df, span = 0.5) # Fit loess smoother
    lat_seq <- seq(min(df$decimallatitude), max(df$decimallatitude), length.out = 200) # Get latitude to predict on
    pred <- predict(fit, newdata = data.frame(decimallatitude = lat_seq)) # Predict
    lines(pred, lat_seq, col = pal[x], lwd = 3) # Plot
    return(NULL)
  }) # End lapply
  
  abline(h = seq(-90,90,30), lty = "dashed")
  
} # end for loop

# END
