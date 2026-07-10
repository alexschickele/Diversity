
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

# --- 3.4. Latitudinal profile plot
# --- 3.4.1. Prepare data
profile_year <- apply(data, 1, mean, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T)
profile_summer <- apply(data[,,,c(6,7,8),], 1, mean, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T)
profile_winter <- apply(data[,,,c(12,1,2),], 1, mean, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T)
profile_sd <- apply(data, c(1,3), mean, na.rm = T) %>% 
  apply(., 1, sd, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T)

# --- 3.4.2. Plot
par(mfrow = c(1,3))
plot(1, 1, type = 'n', xlim = c(0, 80), ylim = c(-90, 90), xlab = "Diversity index (%)", ylab = "Latitude", axes = FALSE)
# id <- complete.cases(c(profile_year-profile_sd, rev(profile_year+profile_sd)))
# polygon(x = c(profile_year-profile_sd*2, rev(profile_year+profile_sd*2))[id],
#         y = c(89.5:-89.5, -89.5:89.5)[id], col = "gray80", border = NA)
lines(x = moving_average(x = profile_year, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = "black", lwd = 3)
lines(x = moving_average(x = c(profile_summer[1:90], profile_winter[91:180]), n = 10) %>% as.numeric(), y = 89.5:-89.5, col = "red4", lwd = 2)
lines(x = moving_average(x = c(profile_winter[1:90], profile_summer[91:180]), n = 10) %>% as.numeric(), y = 89.5:-89.5, col = "skyblue3", lwd = 2)

abline(h = c(60, 30, 0, -30, -60), lty = "dotted")
axis(side = 1)
axis(side = 2, at = c(60, 30, 0, -30, -60), las = 2)
box()

# --- 3.5. All member profiles
# --- 3.5.1. Setup and consensus profile
par(mfrow = c(1,3))
plot(1, 1, type = 'n', xlim = c(0, 100), ylim = c(-90, 90), xlab = "Diversity index (%)", ylab = "Latitude", axes = FALSE)
profile_year <- apply(data, 1, mean, na.rm = T) %>% array(., dim = c(360,180)) %>% apply(., 2, mean, na.rm = T) # Consensus profile

# --- 3.5.2. Add all members
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
plot(r_rob, col = rocket_pal(100), axes = FALSE, main = "global SD", range = c(0,10), fill_range = T)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box







# END