# This script processes CEPHALOPOD outputs to extract global diversity patterns
# and their overlap with existing Marine Protected Area networks
# Corresponds to the analyses related to Figure 1 in the manuscript.

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
robinson_proj <- "+proj=robin +lon_0=210 +datum=WGS84"

# --- 1.6. Project land mask
land <- ne_countries(scale = 50, returnclass = "sf") %>% .[,1]
land_raster <- terra::rasterize(land, terra::rast(res = 0.1))  # 0.1-degree resolution
land_rob <- terra::project(land_raster, robinson_proj)
land_coarse <- terra::aggregate(land_raster, 10) # match res for MPA coverage

# --- 1.7. Quantile scale
# To make sure its correct across all data types
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

# --- 1.8. Helper Function
# Moving average function
moving_average <- function(x, n = 10) {
  if (length(x) > n) stats::filter(x, rep(1 / n, n), sides = 2) else NA
}

# --- 1.9. Hill clean names
hill_ref <- paste("Hill", seq(0, 5, 0.25))

# --- 1.10. MPA shape file
shp0 <- st_read("/net/meso/work/aschickele/Diversity/data/MPA/shp0/WDPA_WDOECM_Jun2026_Public_marine_shp-polygons.shp")
shp1 <- st_read("/net/meso/work/aschickele/Diversity/data/MPA/shp1/WDPA_WDOECM_Jun2026_Public_marine_shp-polygons.shp")
shp2 <- st_read("/net/meso/work/aschickele/Diversity/data/MPA/shp2/WDPA_WDOECM_Jun2026_Public_marine_shp-polygons.shp")
shp_all <- rbind(shp0, shp1, shp2)
rm(shp0, shp1, shp2)
gc()

# --- 2. Load data
# --- 2.1. Metagenomics
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/METAGENOMICS_raw_diversity.RData")
m <- apply(data, -1, quantile_scale) # rescale
dimnames(m)[[2]] <- sub(" \\(.*", "", dimnames(m)[[2]])

rm(data)
gc()

# --- 2.2. Traditional abundance
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_ABUNDANCE_raw_diversity.RData")
a <- apply(data_abundance, -1, quantile_scale) # rescale
dimnames(a)[[2]] <- sub(" \\(.*", "", dimnames(a)[[2]])

# --- 2.2.1. Interpolate the Hill 3 to be fully factorial (didnt pass QC)
# It just avoids weird averages later to disantangle all factors - no consequence on results
a_hill3 <- apply(a[,12:13,,], c(1,3,4), mean, na.rm = T)
a <- abind(a[,1:12,,], a_hill3, a[,13:20,,], along = 2)
dimnames(a)[[2]][13] <- "Hill 3"

rm(data_abundance, a_hill3)
gc()

# --- 2.3. Traditional biomass
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_BIOMASS_raw_diversity.RData")
b <- apply(data_biomass, -1, quantile_scale) # rescale
dimnames(b)[[2]] <- sub(" \\(.*", "", dimnames(b)[[2]])

rm(data_biomass)
gc()

# --- 2.4. Occurrence
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/OCCURRENCE_raw_diversity.RData")
o <- apply(data, -1, quantile_scale) # rescale
dimnames(o)[[2]] <- paste("Hill", dimnames(o)[[2]])
rm(data)
gc()

# --- 2.5. All together - now we can apply across all factors
data <- abind(m, a, b, o, along = 5)
dimnames(data)[[5]] <- c("m","a","b","o")
rm(m,a,b,o)
gc()

# --- 3 Global average
par(mfrow = c(1,1))

# --- 3.1. Prepare the rasters
r <- r0
terra::values(r) <- apply(data, 1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# --- 3.2. Plot it
plot(r_rob, col = parula_pal(100), axes = FALSE, main = "global average")
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 3.3. Latitudinal gradient test
r_df <- as.data.frame(r, xy = TRUE)
cor.test(r_df$y, r_df$lyr.1, method = "spearman")

# --- 3.4. All member profiles
# --- 3.4.1. Setup and consensus profile
par(mfrow = c(1,3))
plot(1, 1, type = 'n', xlim = c(0, 100), ylim = c(-90, 90), xlab = "Diversity index (%)", ylab = "Latitude", axes = FALSE)
profile_year <- apply(data, 1, mean, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T) # Consensus profile

# --- 3.4.2. Add all members
data_collapsed <- apply(data, c(1,2,4,5), mean, na.rm = T) # remove bootstrap for plot clarity
data_collapsed <- array(data_collapsed, dim = c(64800, 21*12*4))  # 2D
for(x in 1:dim(data_collapsed)[2]){
  profile_member <- data_collapsed[,x] %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T) # Member profile
  lines(x = moving_average(x = profile_member, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = scales::alpha("black", 0.01), lwd = 1) # Member profile
} # end for loop

lines(x = moving_average(x = profile_year, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = "black", lwd = 3) # Consensus profile

abline(h = c(60, 30, 0, -30, -60), lty = "dotted")
axis(side = 1)
axis(side = 2, at = c(60, 30, 0, -30, -60), las = 2)
box()

# --- 4. Global bootstrap uncertainty
par(mfrow = c(1,1))

# --- 4.1. Prepare the rasters
r <- r0
terra::values(r) <- apply(data, c(1,2,4,5), sd, na.rm = T) %>% apply(1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# --- 4.2. Plot it
plot(r_rob, col = rocket_pal(100) %>% rev(), axes = FALSE, main = "global SD", range = c(0,20), fill_range = T)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 5. MPA coverage of the diversity
# --- 5.1. Rasterize the MPAs
r_mpa_all <- shp_all %>% 
  vect(.) %>% rasterize(., r0, field = 1, background = 0)
r_mpa_all[land_coarse == 1] <- NA

r_mpa_no_take <- shp_all %>% 
  dplyr::filter(NO_TAKE == "All" | NO_TAKE == "Part") %>% 
  vect(.) %>% rasterize(., r0, field = 1, background = 0)
r_mpa_no_take[land_coarse == 1] <- NA

# --- 5.2. Initialize the matrices
cutoff <- c(1, seq(60,85,1))
cell_sizes <- as.vector(cellSize(r))
data_mean <- apply(data, c(1,3) , mean, na.rm = TRUE) 

# --- 5.3. Define function
get_cover <- function(PROJ, MASK, CUTOFF){
  mask_area <- sum(MASK*cell_sizes, na.rm = T)
  
  proj <- PROJ
  proj[proj < CUTOFF] <- NA
  proj[!is.na(proj)] <- 1
  
  proj_within_mask <- proj
  proj_within_mask[MASK != 1] <- NA

  proj_area <- sum(proj*cell_sizes, na.rm = T)
  proj_within_mask_area <- sum(proj_within_mask*cell_sizes, na.rm = T)

  out <- proj_within_mask_area/proj_area # how much hotspot are covered by MPA
  return(out)
} # end function

# --- 5.4. Loop over the output matrix
mpa_cover_all <- mclapply(cutoff, function(x){
  tmp <- apply(data_mean, 2, function(y)(y = get_cover(PROJ = y, MASK = as.vector(r_mpa_all), CUTOFF = x)))
}, mc.cores = 10) %>% abind(along = 2)
mpa_cover_no_take <- mclapply(cutoff, function(x){
  tmp <- apply(data_mean, 2, function(y)(y = get_cover(PROJ = y, MASK = as.vector(r_mpa_no_take), CUTOFF = x)))
}, mc.cores = 10) %>% abind(along = 2)

# --- 5.5. Plot
# --- 5.5.1. Boxplot per cutoff
# Cutoff indices and maximum cutoff coverage
cutoff_indices <- array(rep(1:dim(mpa_cover_all)[2], each = dim(mpa_cover_all)[1]),
  dim = dim(mpa_cover_all))
id <- which.max(apply(mpa_cover_no_take, 2, median, na.rm = T))

par(mfrow = c(1,3))
a <- boxplot(mpa_cover_all*100 ~ cutoff_indices, axes = F, col = scales::alpha("#FFB400", 0.5),
             xlab = "Cutoff value for diversity hotspot (%)", ylab = "Area coverage (%)")
b <- boxplot(mpa_cover_no_take*100 ~ cutoff_indices, axes = F, col = scales::alpha("#7AD151FF", 0.5), add = T)
axis(side = 1, at = c(1, seq(2, length(cutoff), 5)), labels = c("Global", seq(60,85,5)), las = 2)
axis(side = 2, las = 2)
abline(v = 1, h = a$stats[3,1], col = "#FFB400", lty = "longdash")
abline(v = id, h = a$stats[3,id], col = "#FFB400")
abline(v = 1, h = b$stats[3,1], col = "#7AD151FF", lty = "longdash")
abline(v = id, h = b$stats[3,id], col = "#7AD151FF")
abline(v = 1)
box()

# --- 5.5.2. Spatial location
# Prepare the rasters
r <- r0
terra::values(r) <- apply(data, 1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)
r_rob[r_rob < cutoff[id]] <- 0
r_rob[r_rob > 0] <- 1

r_mpa_all_rob <- terra::project(r_mpa_all, robinson_proj)
r_mpa_no_take_rob <- terra::project(r_mpa_no_take, robinson_proj)

# Plot it
par(mfrow = c(1,1))
plot(r_rob, col = scales::alpha("red", c(0, 0.5)), axes = FALSE, main = "MPA coverage")
contour(r_rob, add = TRUE, nlevels = 2, col = "red", label = "") 

plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
plot(r_mpa_all_rob, add = TRUE, col = scales::alpha("#FFB400", c(0, 0.3)))
contour(r_mpa_all_rob, add = TRUE, nlevels = 1, col = "#FFB400", label = "")

plot(r_mpa_no_take_rob, add = TRUE, col = scales::alpha("#7AD151FF", c(0, 0.3)))
contour(r_mpa_no_take_rob, add = TRUE, nlevels = 1, col = "#7AD151FF", label = "")

grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# END