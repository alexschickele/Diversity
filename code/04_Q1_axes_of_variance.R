
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

# --- 1.8. Helper Function
# Moving average function
moving_average <- function(x, n = 10) {
  if (length(x) > n) stats::filter(x, rep(1 / n, n), sides = 2) else NA
}

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
metadata <- expand.grid(method = dimnames(data)[[2]],
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
plot(sites_scores[,1], sites_scores[,2], 
     col = scales::alpha(rocket_pal(100) %>% rev(), 0.1)[cut(site_metadata$lat %>% abs(), breaks = 100)], 
     pch = 16, cex = 3, axes = F, xlab = "", ylab = "", xlim = c(-1.5, 1.5), ylim = c(-3,3))
abline(h = 0, v = 0, lwd = 1, lty = "dotted")

mtext(paste("PC1 (", round(PC_percent[1], 2)*100,"% )"), side = 4, line = 1, at = 0, las = 2)
mtext(paste("PC2 (", round(PC_percent[2], 2)*100,"% )"), side = 3, line = 1, at = 0)
mtext("Principal Component Analysis \n (Geographical cell f. of Diversity)", side = 3, line = 3, at = 0)

# --- 3.4.3. Factor points
factor_pal <- c("chocolate2", "#CCC591", "darkolivegreen4", "antiquewhite4") # wes anderson
points(variable_centroids$X, variable_centroids$Y, pch = 22, bg = factor_pal, cex = 1.5)
text(x = variable_centroids$X, y = variable_centroids$Y, labels = gsub("TRADITIONAL_", "", variable_centroids$FACTOR), pos = c(1,3,1,1), cex = 0.7)

# --- 3.4.4. Environmental arrows
arrow_names <- c("PRIMARY \n PRODUCTION","                NITRATE","TEMPERATURE              ")
lapply(seq_along(arrow_names), function(x){
  arrows(0,0, SUP$vectors$arrows[x+1, "PC1"], SUP$vectors$arrows[x+1, "PC2"], length = 0.1, col = "black")
  text(x = SUP$vectors$arrows[x+1, "PC1"], y = SUP$vectors$arrows[x+1, "PC2"], labels = arrow_names[x], pos = c(1,1,2)[x], cex = 0.7)
}) # end lapply

# --- 3.5. Plot density on PC2 - by data type
# --- 3.5.1. Prepare the data
sites_profile <- data
colnames(sites_profile) <- metadata$DATA_TYPE # name column after data type
sites_profile <- sapply(unique(colnames(sites_profile)), function(x) {
  rowMeans(sites_profile[, colnames(sites_profile) == x, drop = FALSE]) # average across identical data type
}) %>% cbind(sites_scores[, 2]) %>% as.data.frame()
colnames(sites_profile)[5] <- "X"

# --- 3.5.2 Iterative plot
profile_names <- variable_centroids$FACTOR
factor_pal <- c("chocolate2", "#CCC591", "darkolivegreen4", "antiquewhite4") # wes anderson

par(mfrow = c(1,4), mar = c(3,3,1,1))
plot(x = 1, y = 1, type = 'n', xlim = c(0, 100), ylim = c(-3,3), xlab = "Diversity hotspot probability (%)", ylab = "PC2", main = "", axes = FALSE)
lapply(seq_along(profile_names), function(x){
  fit <- glm(formula = as.formula(paste(profile_names[x], "~ poly(X, 2)")), data = sites_profile)
  pred <- predict(fit, newdata = data.frame(X = seq(-3, 3, 0.1)), type = "response")  # or type = "link" if you want linear predictor
  lines(x = pred, y = seq(-3, 3, 0.1), col = factor_pal[x], lwd = 2)
  print(fit$coefficients) # plot information
  anova(fit) # plot information
})
abline(h = 0, lty = "dotted", lwd = 1)
axis(side = 1)

# --- 3.6. Factor importance plot
# --- 3.6.1. Extract the information that DATA_TYPE is the most important one
# As the squared scores across PC1 and PC2
contribution <- variable_scores %>% 
  mutate(contribution = (PC1**2) + (PC2**2)) %>% 
  dplyr::select(-PC1, -PC2)

# Now we compute the average score per factor
contribution <- lapply(c("METHOD","MONTH","DATA_TYPE"), function(x){
  tmp <- contribution %>% 
    group_by_at(x) %>% 
    summarise(FACTOR = x,
              contribution_total = sum(contribution)) %>% 
    ungroup()
  colnames(tmp) <- c("MODALITY","FACTOR","contribution_total")
  return(tmp)
}) %>% bind_rows() %>% # end lapply
  group_by(FACTOR) %>% 
  summarise(contribution_mean = mean(contribution_total)) %>% # total contrib per factor
  ungroup() %>% 
  mutate(contribution_perc = contribution_mean / sum(contribution_mean)) # divide by total across factor

# --- 3.6.2. Plot
par(mfrow = c(1,2), mar = c(10,5,5,5))
plot(x = 1:3, y = contribution$contribution_perc, pch = 16, cex = 4, ylim = c(0,1), xlim = c(0,4), axes = FALSE, xlab = "", ylab = "Contribution (%)")
segments(x0 = 1:3, x1 = 1:3, y0 = rep(0,3), y1 = contribution$contribution_perc, lwd = 2)
axis(side = 1, at = 1:3, labels = contribution$FACTOR, las = 2)
axis(side = 2, at = seq(0, 1, 0.2), las = 1)
box()

# --- 4. Map the differences
par(mfrow = c(2,2), mar = c(0, 1, 2, 4))

# --- 4.1. Standard deviation across data types
# Prepare the data
val_dt <- lapply(seq_along(data_list), function(x){
  out <- apply(data_list[[x]], 1, mean, na.rm = TRUE) # mean across month, bootstrap and hill number
}) %>% abind(along = 2) 

arr_dt <- array(as.matrix(val_dt), dim = c(360, 180, 4)) # array for ANOVA
val_dt <- val_dt %>% apply(1, sd, na.rm = TRUE) # sd across data type

# Significancy test
signif_dt <- moving_window_anova(A = arr_dt, R = 3, THRESHOLD = 0.05) %>% 
  st_as_sf(., coords = c(1,2), crs = 4326) %>% st_transform(., crs = robinson_proj)

# Plot
setValues(r0, val_dt) %>% terra::project(robinson_proj) %>% plot(col = rocket_pal(100) %>% rev(), range = c(0,50), axes = F, main = "Hotspot probability disagreement \n (SD across data type)")
setValues(r0, val_dt) %>% terra::project(robinson_proj) %>% contour( add = TRUE, nlevels = 5)
points(signif_dt, pch = 16, col = scales::alpha("black", 0.5), cex = 0.5) # significancy test
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20") # land
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# --- 4.2. Standard deviation across hill numbers
# Prepare the data
tmp <- abind(data_list, along = 2)
val_h <- lapply(dimnames(tmp)[[2]] %>% unique(), function(x){
  id <- which(dimnames(tmp)[[2]] == x)
  out <- tmp[,id,,] %>% apply(1, mean, na.rm = TRUE) # mean across month, bootstrap and data type
}) %>% abind(along = 2) 

arr_h <- val_h %>% apply(1, function(x)(x = c(mean(x[1:18], na.rm = TRUE), mean(x[4:21], na.rm = TRUE)))) %>% 
  t() %>% array(., dim = c(360, 180, 2)) # array for ANOVA
val_h <- val_h %>% apply(1, function(x)(x = mean(x[1:18] - x[4:21], na.rm = TRUE))) # sd across hill numbers

# Significancy test
signif_h <- moving_window_anova(A = arr_h, R = 3, THRESHOLD = 0.05) %>% 
  st_as_sf(., coords = c(1,2), crs = 4326) %>% st_transform(., crs = robinson_proj)

# Plot
setValues(r0, val_h) %>% terra::project(robinson_proj) %>% plot(col = curl_pal(100) %>% rev(), range = c(-5,5), axes = F, main = "Hotspot probability change \n per unit Hill scaling factor")
setValues(r0, val_h) %>% terra::project(robinson_proj) %>% contour( add = TRUE, nlevels = 5)
points(signif_h, pch = 16, col = scales::alpha("black", 0.5), cex = 0.5) # significancy test
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20") # land
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# --- 4.3. Standard deviation across month
# Prepare the data
val_m <- abind(data_list, along = 2) %>% 
  apply(c(1,4), mean, na.rm = TRUE) # mean across bootstrap, hill number and data type

arr_m <- val_m %>% apply(1, function(x)(x = c(mean(x[7:9]), mean(x[1:3])))) %>% 
  t() %>% array(., dim = c(360, 180, 2))
val_m <- val_m %>% apply(1, function(x)(x = mean(x[7:9]) - mean(x[1:3]))) # summer - winter

# Significancy test
signif_m <- moving_window_anova(A = arr_m, R = 3, THRESHOLD = 0.05) %>% 
  st_as_sf(., coords = c(1,2), crs = 4326) %>% st_transform(., crs = robinson_proj)

# Plot
setValues(r0, val_m) %>% terra::project(robinson_proj) %>% plot(col = curl_pal(100) %>% rev(), range = c(-50,50), axes = F, main = "Hotspot probability \n (summer - winter)")
setValues(r0, val_m) %>% terra::project(robinson_proj) %>% contour(add = TRUE, nlevels = 5)
points(signif_m, pch = 16, col = scales::alpha("black", 0.5), cex = 0.5) # significancy test
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20") # land
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

# --- 5. Map the profiles
# Taking the data from the previous generated arrays
par(mfrow = c(1,3))

# --- 5.1. Data type
# --- Prepare data
arr_dt <- lapply(seq_along(data_list), function(x){
  out <- apply(data_list[[x]], 1, mean, na.rm = TRUE) # mean across month, bootstrap and hill number
}) %>% abind(along = 2) %>%  array(., dim = c(360, 180, 4)) # array for ANOVA

# --- Plot
pal <- brewer.pal(4, "BrBG")
plot(1, 1, type = 'n', xlim = c(20, 80), ylim = c(-90, 90), xlab = "Hotspot probability", ylab = "Latitude", axes = FALSE)
lapply(1:dim(arr_dt)[[3]], function(x){
  val <- arr_dt[,,x] %>% apply(c(2), mean, na.rm = TRUE)
  lines(x = moving_average(x = val, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = pal[x], lwd = 2)
}) # end lapply
abline(h = c(60, 30, 0, -30, -60), lty = "dotted")
axis(side = 1)
axis(side = 2, at = c(60, 30, 0, -30, -60), las = 2)
box()

# --- 5.2. Hill
# --- Prepare the data
tmp <- abind(data_list, along = 2)
arr_h <- lapply(dimnames(tmp)[[2]] %>% unique(), function(x){
  id <- which(dimnames(tmp)[[2]] == x)
  out <- tmp[,id,,] %>% apply(1, mean, na.rm = TRUE) # mean across month, bootstrap and data type
}) %>% abind(along = 2) %>%  array(., dim = c(360, 180, 21)) # array

# --- Plot
pal <- rocket_pal(21)
plot(1, 1, type = 'n', xlim = c(20, 80), ylim = c(-90, 90), xlab = "Hotspot probability", ylab = "Latitude", axes = FALSE)
lapply(1:dim(arr_h)[[3]], function(x){
  val <- arr_h[,,x] %>% apply(c(2), mean, na.rm = TRUE)
  lines(x = moving_average(x = val, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = pal[x], lwd = 1)
}) # end lapply
abline(h = c(60, 30, 0, -30, -60), lty = "dotted")
axis(side = 1)
axis(side = 2, at = c(60, 30, 0, -30, -60), las = 2)
box()

# --- 5.3. Month
# --- Prepare data
arr_m <- abind(data_list, along = 2) %>% 
  apply(c(1,4), mean, na.rm = TRUE) %>% array(., dim = c(360, 180, 12))

# --- Plot
pal <- circular_pal(12)
plot(1, 1, type = 'n', xlim = c(20, 80), ylim = c(-90, 90), xlab = "Hotspot probability", ylab = "Latitude", axes = FALSE)
lapply(1:dim(arr_m)[[3]], function(x){
  val <- arr_m[,,x] %>% apply(c(2), mean, na.rm = TRUE)
  lines(x = moving_average(x = val, n = 10) %>% as.numeric(), y = 89.5:-89.5, col = pal[x], lwd = 1)
}) # end lapply
abline(h = c(60, 30, 0, -30, -60), lty = "dotted")
axis(side = 1)
axis(side = 2, at = c(60, 30, 0, -30, -60), las = 2)
box()

# --- 6. Global uncertainty
# Standard deviation on the bootstrap across all factors
par(mfrow = c(1,1))
val_sd <- abind(data_list, along = 2) %>% apply(c(1,2,4), sd, na.rm = TRUE) %>% apply(1, mean, na.rm = TRUE)

# Plot
setValues(r0, val_sd) %>% terra::project(robinson_proj) %>% plot(col = rocket_pal(100) %>% rev(), range = c(0,20), axes = F, main = "Global uncertainty \n (SD across all bootstraps)")
setValues(r0, val_sd) %>% terra::project(robinson_proj) %>% contour(add = TRUE, nlevels = 5)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20") # land
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
box("figure", col="black", lwd = 1) # box

#rda# END
