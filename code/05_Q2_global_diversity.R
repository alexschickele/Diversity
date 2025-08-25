
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

# --- 1.8. Hill clean names
hill_ref <- paste("Hill", seq(0, 5, 0.25))

# --- 1.9. Open Biomes
# Based on Fay and McKinley at a coarser resolution
biomes <- nc_open("./data/Time_Varying_Biomes.nc")
biomes <- ncvar_get(biomes, "MeanBiomes") %>% t()

biomes <- terra::setValues(r, as.vector(biomes)) %>% flip()
biomes[biomes %in% c(16, 17)] <- 101 # so
biomes[biomes %in% c(1,8)] <- 102 # arctic
biomes[biomes %in% c(2,3, 9,10,15)] <- 103 # transition
biomes[biomes %in% c(6, 12)] <- 104 # westerly
biomes[biomes %in% c(4,5,7,8)] <- 105 # pacific
biomes[biomes %in% c(11,13)] <- 106 # atlantic
biomes[biomes %in% c(14)] <- 107 # indian

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

# --- 3. Global variance explained
# --- 3.1. Occurrence
load(paste0("./output/", FOLDER_NAME, "/OCCURRENCE_vip.RData"))
VIP <- vip %>% mutate(source = "occurrence")

# --- 3.2. Abundance
load(paste0("./output/", FOLDER_NAME, "/TRADITIONAL_ABUNDANCE_vip.RData"))
VIP <- bind_rows(VIP, vip_abundance %>% mutate(source = "abundance"))

# --- 3.3. Biomass
load(paste0("./output/", FOLDER_NAME, "/TRADITIONAL_BIOMASS_vip.RData"))
VIP <- bind_rows(VIP, vip_biomass %>% mutate(source = "biomass"))

# --- 3.4. Metagenomics
load(paste0("./output/", FOLDER_NAME, "/METAGENOMICS_vip.RData"))
VIP <- bind_rows(VIP, vip %>% mutate(source = "metagenomics"))

# --- 3.5. Cluster environmental variables
# At the global scale to ease interpretation on the maps
load("/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524/CALL.RData") # load omic' CALL
features <- CALL$ENV_DATA

# --- 3.5.1. Reshape as array cell * layer * month
features_array <- lapply(1:12, function(x)(x = features[[x]] %>% unwrap() %>% as.matrix())) %>% abind(along = 3)
dimnames(features_array) <- list(NULL, names(features[[1]] %>% unwrap()), as.character(1:12))

# --- 3.5.2. Yearly average
features_array_year <- apply(features_array, c(1,2), mean, na.rm = T)

# --- 3.5.3. Correlation clustering
library(dendextend)
m <- cor(features_array_year, method = "spearman", use = "pairwise.complete.obs") # Compute correlation matrix
d <- as.dist(1 - m) # Convert correlation to distance matrix (1 - correlation)
hc <- hclust(d, method = "complete") %>% cutree(4) # Hierarchical clustering
group <- data.frame(variable = names(hc), group = hc)

# --- 3.6. Prepare boxplot data
# We sum the variance explained per group and data source
df <- VIP %>% 
  left_join(group) %>% 
  group_by(source, group) %>% 
  summarise(value = sum(value)) %>% 
  ungroup()

# --- 3.7. Plot
boxplot(df$value ~ df$group, horizontal = TRUE, col = "antiquewhite3")
abline(v = seq(0, 1, 0.2), lty = "dotted", col = "gray20")


# Plot the dendrogram
par(mar = c(10,10,10,10), mfrow = c(1,3))
plot(as.dendrogram(hc), main = "Dendrogram of features", horiz = TRUE)













# --- 3. Plot
# --- 3.1. Global average
par(mfrow = c(1,1))

# Prepare the rasters
r <- r0
terra::values(r) <- abind(o,a,b,m, along = 2) %>% apply(1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# Plot it
plot(r_rob, col = parula_pal(100), axes = FALSE, main = "global average")
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 4. Extract diversity shape and composition
# --- 4.1. Extract coordinates of < 20 and > 80 quantiles
data <- abind(o,a,b,m, along = 2) %>% apply(1, mean, na.rm = TRUE) # global average
top20 <- which(data > quantile(data, 0.8, na.rm = TRUE))
bottom20 <- which(data < quantile(data, 0.2, na.rm = TRUE))


# Load omic dummy data
df <- vroom(paste0("./output/", FOLDER_NAME, "/METAGENOMIC_raw_input.csv"))











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