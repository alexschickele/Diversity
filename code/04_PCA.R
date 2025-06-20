
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

rm(land, land_raster)
gc()

# --- 1.7. Quantile scale
# To make sure its correct across all data types
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

# --- 2. Load data
# --- 2.1. Metagenomics
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/METAGENOMICS_raw_diversity.RData")
tmp <- apply(data, -1, quantile_scale) # rescale
data <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(data)[[2]] <- sub(" \\(.*", "", dimnames(data)[[2]])

data_list <- list(m = data) # start saving in a list for maps later

m <- apply(data, 1, c)  # Collapse
metadata <- expand.grid(method = sub(" \\(.*", "", dimnames(data)[[2]]),
                        bootstrap = dimnames(data)[[3]],
                        month = dimnames(data)[[4]],
                        data_type = "METAGENOMICS")
m <- cbind(m, metadata) # rows = diversity estimates, columns = cells + factors
rm(data)
gc()

# --- 2.2. Traditional abundance
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_ABUNDANCE_raw_diversity.RData")
tmp <- apply(data_abundance, -1, quantile_scale) # rescale
data_abundance <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(data_abundance)[[2]] <- sub(" \\(.*", "", dimnames(data_abundance)[[2]])

data_list[["a"]] <- data_abundance # start saving in a list for maps later

a <- apply(data_abundance, 1, c)  # Collapse
metadata <- expand.grid(method = sub(" \\(.*", "", dimnames(data_abundance)[[2]]),
                        bootstrap = dimnames(data_abundance)[[3]],
                        month = dimnames(data_abundance)[[4]],
                        data_type = "TRADITIONAL_ABUNDANCE")
a <- cbind(a, metadata) # rows = diversity estimates, columns = cells + factors
rm(data_abundance)
gc()

# --- 2.3. Traditional biomass
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/TRADITIONAL_BIOMASS_raw_diversity.RData")
tmp <- apply(data_biomass, -1, quantile_scale) # rescale
data_biomass <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(data_biomass)[[2]] <- sub(" \\(.*", "", dimnames(data_biomass)[[2]])

data_list[["b"]] <- data_biomass # start saving in a list for maps later

b <- apply(data_biomass, 1, c)  # Collapse
metadata <- expand.grid(method = sub(" \\(.*", "", dimnames(data_biomass)[[2]]),
                        bootstrap = dimnames(data_biomass)[[3]],
                        month = dimnames(data_biomass)[[4]],
                        data_type = "TRADITIONAL_BIOMASS")
b <- cbind(b, metadata) # rows = diversity estimates, columns = cells + factors
rm(data_biomass)
gc()

# --- 2.4. Occurrence
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/OCCURRENCE_raw_diversity.RData")
tmp <- apply(data, -1, quantile_scale) # rescale
data <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap
dimnames(data)[[2]] <- paste("Hill", dimnames(data)[[2]])

data_list[["o"]] <- data # start saving in a list for maps later

o <- apply(data, 1, c)  # Collapse
metadata <- expand.grid(method = paste("Hill", dimnames(data)[[2]]),
                        bootstrap = dimnames(data)[[3]],
                        month = dimnames(data)[[4]],
                        data_type = "OCCURRENCE")
o <- cbind(o, metadata) # rows = diversity estimates, columns = cells + factors
rm(data)
gc()

# --- 2.5. Concatenate
# Bind all individual diversity estimates as rows for PCA
data <- rbind(m, a, b, o)
rm(m, a, b, o, metadata)
gc()

# Add column names for the cells
colnames(data) <- c(paste0("C", 1:64800), "METHOD","BOOTSTRAP","MONTH","DATA_TYPE")

# Remove NA columns
id <- which(is.na(apply(data[, 1:64800], 2, sum)))
data <- data[,-id]

# --- 2.6. Split PCA and RDA data
metadata <- data[, (ncol(data)-3):ncol(data)]
data <- data[, 1:(ncol(data)-4)] %>% t()

# --- 2.7. Add cell metadata
r1 <- terra::rast("/net/meso/work/clercc/Predictors/PIPELINE_SET/climatology_S_PP_regridded.nc") %>% mean()
r2 <- terra::rast("/net/meso/work/clercc/Predictors/PIPELINE_SET/climatology_n_0_50.nc") %>% mean()
r3 <- terra::rast("/net/meso/work/clercc/Predictors/PIPELINE_SET/climatology_t_0_50.nc") %>% mean()
r <- terra::rast(list(r1, r2, r3))

r_point <- extract(x = r, y = 1:64800, xy = TRUE) %>% as.data.frame()

colnames(r_point) <- c("lon","lat","PP", "Nitrate", "Temperature")
site_metadata <- r_point[-id, -1]

rm(r1, r2, r3, r, r_point)
gc()

# --- 3. PCA on diversity estimates
# --- 3.1. Perform PCA // Species = diversity ; sites = geographical cells
# Site metadata such as environmental variables or latitude are supplementary
PCA <- vegan::pca(X = data, scale = T)
SUP <- vegan::envfit(PCA, site_metadata)

# --- 3.2. Extract contribution of each Principal Component (PC)
PC_percent <- (eigenvals(PCA)/sum(eigenvals(PCA))) # variance per PC

# --- 3.3. Extract scores
variable_scores <- vegan::scores(PCA, choices = c(1,2), display = "species") %>% 
  cbind(metadata) %>% # association
  group_by(METHOD, MONTH, DATA_TYPE) %>% 
  summarize(PC1 = mean(PC1),
            PC2 = mean(PC2)) %>% 
  ungroup()

variable_centroids <- lapply(c("DATA_TYPE"), function(x){
  out <- variable_scores %>% 
    group_by_at(x) %>% 
    summarise(X = mean(PC1), Y = mean(PC2)) %>% 
    ungroup() %>% as.data.frame()
  colnames(out) <- c("FACTOR","X","Y")
  return(out)
}) %>% bind_rows() # end lapply

sites_scores <- vegan::scores(PCA, choices = c(1,2), display = "sites")

# --- 3.4. Plotting the PCA
# --- 3.4.1. Extract the axis correlation
env_cor <- apply(site_metadata, 2, function(var) {cor(var, sites_scores, method = "pearson")})

# --- 3.4.2. Plot cloud of points
par(mfrow = c(1,1), mar = c(1,5,6,8))
plot(sites_scores[,1], sites_scores[,2], col = scales::alpha("coral", 0.1), 
     pch = 16, cex = 1, axes = F, xlab = "", ylab = "")
abline(h = 0, v = 0, lwd = 1, lty = "dotted")
# arrows(0,0,1.25,0, lwd = 2, length = 0.1, lty = "dotted")
# arrows(0,0,0, 2.8, lwd = 2, length = 0.1, lty = "dotted")

mtext(paste("PC1 (", round(PC_percent[1], 2)*100,"% )"), side = 4, line = 1, at = 0, las = 2)
mtext(paste("PC2 (", round(PC_percent[2], 2)*100,"% )"), side = 3, line = 1, at = 0)
mtext("Principal Component Analysis \n (Geographical cell f. of Diversity)", side = 3, line = 3, at = 0)

# --- 3.4.3. Factor points
points(variable_centroids$X, variable_centroids$Y, pch = 15)
text(x = variable_centroids$X, y = variable_centroids$Y, labels = variable_centroids$FACTOR, pos = 2, cex = 0.7)

# --- 3.4.4. Environmental arrows
lapply(c("PP","Nitrate","Temperature"), function(x){
  arrows(0,0, SUP$vectors$arrows[x, "PC1"], SUP$vectors$arrows[x, "PC2"], length = 0.1, col = "blue")
  text(x = SUP$vectors$arrows[x, "PC1"], y = SUP$vectors$arrows[x, "PC2"], labels = x, pos = 2, cex = 0.7)
}) # end lapply

# --- 4. Map the differences
par(mfrow = c(2,2), mar = c(0, 1, 2, 4))

# --- 4.1. Standard deviation across data types
# Prepare the data
sd_dt <- lapply(seq_along(data_list), function(x){
  out <- apply(data_list[[x]], 1, mean, na.rm = TRUE) # mean across month, bootstrap and hill number
}) %>% abind(along = 2) %>% apply(1, sd, na.rm = TRUE) # sd across data type

# Plot
setValues(r0, sd_dt) %>% terra::project(robinson_proj) %>% plot(col = rocket_pal(100) %>% rev(), range = c(0,50), axes = F, main = "Hotspot probability disagreement \n (SD across data type)")
setValues(r0, sd_dt) %>% terra::project(robinson_proj) %>% contour( add = TRUE, nlevels = 5)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# --- 4.2. Standard deviation across hill numbers
# Prepare the data
tmp <- abind(data_list, along = 2)
sd_h <- lapply(dimnames(tmp)[[2]] %>% unique(), function(x){
  id <- which(dimnames(tmp)[[2]] == x)
  out <- tmp[,id,,] %>% apply(1, mean, na.rm = TRUE) # mean across month, bootstrap and data type
}) %>% abind(along = 2) %>% apply(1, function(x)(x = mean(x[1:18] - x[4:21], na.rm = TRUE))) # sd across hill numbers

# Plot
setValues(r0, sd_h) %>% terra::project(robinson_proj) %>% plot(col = curl_pal(100) %>% rev(), range = c(-5,5), axes = F, main = "Hotspot probability change \n per unit Hill scaling factor")
setValues(r0, sd_h) %>% terra::project(robinson_proj) %>% contour( add = TRUE, nlevels = 5)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# --- 4.3. Standard deviation across month
# Prepare the data
sd_m <- abind(data_list, along = 2) %>% 
  apply(c(1,4), mean, na.rm = TRUE) %>% # mean across bootstrap, hill number and data type
  apply(1, function(x)(x = mean(x[1:3]) - mean(x[7:9]))) # sd across month

# Plot
setValues(r0, sd_m) %>% terra::project(robinson_proj) %>% plot(col = curl_pal(100) %>% rev(), range = c(-50,50), axes = F, main = "Hotspot probability \n (summer - winter)")
setValues(r0, sd_m) %>% terra::project(robinson_proj) %>% contour(add = TRUE, nlevels = 5)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box





# Significancy test






#rda# END
