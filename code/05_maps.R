
# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME <- "DIVERSITY_PAPER"

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")
 
# --- 1.4. Fake raster for later
r0 <- terra::rast(nrows = 180, ncols = 360, xmin = -180, xmax = 180, ymin = -90, ymax = 90)

# --- 1.5. Robinson projection
robinson_proj <- "+proj=robin +lon_0=150 +datum=WGS84"

# --- 1.6. Project land mask
land <- ne_countries(scale = 50, returnclass = "sf") %>% .[,1]
land_raster <- terra::rasterize(land, terra::rast(res = 0.1))  # 0.1-degree resolution
land_rob <- terra::project(land_raster, robinson_proj)

# --- 1.7. Quantile scale
# To make sure its correct across all data types
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

# --- 1.8. Hill clean names
hill_ref <- paste("Hill", seq(0, 5, 0.25))

# --- 2. Load data
# --- 2.1. Metagenomics
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/METAGENOMICS_raw_diversity.RData")
tmp <- apply(data, -1, quantile_scale) # rescale
m <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(m)[[2]] <- sub(" \\(.*", "", dimnames(m)[[2]])

rm(data)
gc()

# --- 2.2. Traditional abundance
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_ABUNDANCE_raw_diversity.RData")
tmp <- apply(data_abundance, -1, quantile_scale) # rescale
a <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(a)[[2]] <- sub(" \\(.*", "", dimnames(a)[[2]])

rm(data_abundance)
gc()

# --- 2.3. Traditional biomass
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_BIOMASS_raw_diversity.RData")
tmp <- apply(data_biomass, -1, quantile_scale) # rescale
b <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(b)[[2]] <- sub(" \\(.*", "", dimnames(b)[[2]])

rm(data_biomass)
gc()

# --- 2.4. Occurrence
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/OCCURRENCE_raw_diversity.RData")
tmp <- apply(data, -1, quantile_scale) # rescale
o <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(o)[[2]] <- paste("Hill", dimnames(o)[[2]])
rm(data)
gc()

# --- 3. Plot
# --- 3.1. Average Maps
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
par(mfrow = c(2,2), mar = c(1,1,4,1))

lapply(seq_along(to_map), function(x){
  # Prepare the rasters
  r <- r0
  terra::values(r) <- apply(to_map[[x]], 1, mean, na.rm = TRUE)
  r_rob <- terra::project(r, robinson_proj)
  
  # Plot it
  plot(r_rob, col = viridis_pal(100), axes = FALSE, main = names(to_map)[x], legend = T)
  plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
  grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
    vect() %>%  project(robinson_proj) 
  plot(grat, lty = "dotted", add = TRUE) # add grid
  contour(r_rob, add = TRUE, nlevels = 5)
  box("figure", col="black", lwd = 1) # box
  
}) # end lapply

# --- 3.2. Uncertainty Maps
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
par(mfrow = c(2,2), mar = c(1,1,4,1))

lapply(seq_along(to_map), function(x){
  # Prepare the rasters
  r <- r0
  terra::values(r) <- apply(to_map[[x]], c(1,2,4), sd, na.rm = TRUE) %>% apply(1, mean, na.rm = TRUE)
  r_rob <- terra::project(r, robinson_proj)
  
  # Plot it
  plot(r_rob, col = rocket_pal(100) %>% rev(), axes = FALSE, main = names(to_map)[x], range = c(0, 25))
  plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
  grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
    vect() %>%  project(robinson_proj) 
  plot(grat, lty = "dotted", add = TRUE) # add grid
  contour(r_rob, add = TRUE, nlevels = 5)
  box("figure", col="black", lwd = 1) # box
  
}) # end lapply

# --- 3.3. Hill profiles
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
par(mfrow = c(1,4), mar = c(3,3,4,1))

pal <- parula_pal(length(hill_ref)) # one color per hill profile
names(pal) <- hill_ref

lapply(seq_along(to_map), function(x){
  tmp <- to_map[[x]] %>% apply(c(1,2), mean, na.rm = T) # get hill
  lat <- xyFromCell(r0, 1:64800)[,2] # get latitude

  plot(x = 1, y = 1, xlim = c(45,85), ylim = c(-70, 70), type = "n", main = names(to_map)[x])
  
  for(i in 1:ncol(tmp)){
    
    fit <- loess(tmp[,i] ~ lat, span = 0.5) # Fit loess smoother
    lat_seq <- seq(min(lat), max(lat), length.out = 200) # Get latitude to predict on
    pred <- predict(fit, newdata = lat_seq) # Predict
    lines(pred, lat_seq, col = pal[colnames(tmp)[i]], lwd = 2) # Plot
    
  } # end for loop
  
  abline(h = c(-60, -30, 0, 30, 60), lwd = "dotted")
  box("figure", col="black", lwd = 1) # box
  
}) # end lapply

# --- 3.4. Global average
par(mfrow = c(1,1))

# Prepare the rasters
r <- r0
terra::values(r) <- abind(o,a,b,m, along = 2) %>% apply(1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# Plot it
plot(r_rob, col = viridis_pal(100), axes = FALSE, main = "global average")
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 3.5. Global Hill plots
tmp <- abind(o,a,b,m, along = 2) %>% apply(c(1,2), median, na.rm = TRUE) # global average
lat <- xyFromCell(r0, 1:64800)[,2] # get latitude

# Average columns with the same name
tmp <- sapply(unique(colnames(tmp)), function(n) {
  rows <- which(colnames(tmp) == n)
  if (length(rows) == 1) {
    tmp[, rows]
  } else {
    rowMeans(tmp[, rows, drop = FALSE], na.rm = TRUE)
  }
}) %>% as.matrix()

plot(x = 1, y = 1, xlim = c(45,85), ylim = c(-70, 70), type = "n", main = "Average Hill profile")
pal <- inferno_pal(21)

for(i in 1:ncol(tmp)){
  fit <- loess(tmp[,i] ~ lat, span = 0.5) # Fit loess smoother
  lat_seq <- seq(min(lat), max(lat), length.out = 200) # Get latitude to predict on
  pred <- predict(fit, newdata = lat_seq) # Predict
  lines(pred, lat_seq, col = pal[i], lwd = 2) # Plot
} # for hill loop






# END