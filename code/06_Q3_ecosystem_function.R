
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

# --- 3. Load ecosystem properties
data_ecosystem <- list()

# --- 3.1. Marine biodiversity Tittensor
s <- read_sf("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Tittensor&al._2010/01_Data/WCMC-019-PatternsBiodiversity2010-AcrossTaxa.shp")
r_interp <- akima::interp(x = s$X_COORD, y = s$Y_COORD, z = s$AllNorm,
                          xo = -180:179, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["Biodiversity"]] <- terra::rast(r_interp$z %>% t(), extent = c(-180,180,-90,90)) %>% flip() # build raster

# --- 3.2. Mean annual catch of small pelagic fish
# Pelagic data from Watson 2017
load("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Fisheries_Watson2017/catch_data/clim_pelagics<30cm_catches_1990-2019.Rdata")
data_ecosystem[["Catch"]] <- terra::rasterize(x = clim[,2:3] %>% as.matrix(), y = r0, values = clim[,7])

# --- 3.3. Carbon dynamics // NPP
nc <- nc_open("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/DeVries&Weber_2017/Cexp_deVries_2017.nc")
r_interp <- akima::interp(x = ncvar_get(nc, "LON")[,,1], y = ncvar_get(nc, "LAT")[,,1], z = ncvar_get(nc, varid = "NPP") %>% apply(c(1,2), mean),
                          xo = 0:359, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["NPP"]] <- terra::rast(r_interp$z %>% t(), extent = c(0,360,-90,90)) %>% flip() %>% terra::rotate() # build raster

# --- 3.4. Carbon dynamics // FPOC
nc <- nc_open("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/DeVries&Weber_2017/Cexp_deVries_2017.nc")
z_ok <- ncvar_get(nc, varid = "FPOCex") %>% apply(c(1,2), mean)
z_ok[is.na(z_ok)] <- 0
r_interp <- akima::interp(x = ncvar_get(nc, "LON")[,,1], y = ncvar_get(nc, "LAT")[,,1], z_ok,
                          xo = 0:359, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["FPOC"]] <- terra::rast(r_interp$z %>% t(), extent = c(0,360,-90,90)) %>% flip() %>% terra::rotate() # build raster

# --- 3.5. Carbon dynamics // efficiency
data_ecosystem[["Efficiency"]] <- data_ecosystem[["FPOC"]]/data_ecosystem[["NPP"]]

# --- 3.6. Particle size distribution
nc <- nc_open("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Kostadinov&al._2009/PSD_PFT_C_mission_climatology_SeaWiFS_Kostadinov.nc")
r <- terra::rast(nrows = 2160, ncols = 4320, xmin = -180, xmax = 180, ymin = -90, ymax = 90)
terra::values(r) <- ncvar_get(nc, "PSD_slope")
data_ecosystem[["PSD"]] <- terra::aggregate(r, fact = c(12,12))

# --- 3.7. Plot ecosystem properties
par(mfrow = c(3,2))
lapply(seq_along(data_ecosystem), function(x){
  
  # --- 3.7.1. Prepare the rasters
  crs(data_ecosystem[[x]]) <- "+proj=longlat +datum=WGS84"
  r_rob <- terra::project(data_ecosystem[[x]], robinson_proj)
  
  # --- 3.7.2. Plot it
  plot(r_rob, col = parula_pal(100), axes = FALSE, main = names(data_ecosystem)[x], 
       range = quantile(values(r_rob), c(0.05, 0.95), na.rm = TRUE), fill_range = TRUE)
  plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
  grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
    vect() %>%  project(robinson_proj) 
  plot(grat, lty = "dotted", add = TRUE) # add grid
  # contour(r_rob, add = TRUE, nlevels = 5)
  box("figure", col="black", lwd = 1) # box
})

# --- 4. Extract Hill maps
# Geographical cell x hill x data type - loess plot for smoothness
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
data_div <- lapply(to_map, function(x) apply(x, c(1,2), mean, na.rm = TRUE))

BEF_cor <- lapply(seq_along(data_ecosystem), function(x){
  val <- data_ecosystem[[x]] %>% as.vector() # extract ecosystem function values
  lapply(seq_along(data_div), function(i){
      corval <- cor(val, hill_maps[[i]], use = "pairwise.complete.obs", method = "spearman") # Do correlation
      fit <- loess(corval[1,] ~ str_sub(colnames(corval), 6, -1), span = 0.9) # Fit loess smoother
      pred <- predict(fit, newdata = seq(0, 5, 0.25)) # Predict
      return(pred)
  }) %>% bind_cols() %>% t() # data type loop
}) # ecosystem function loop

# --- 5. Correlation plots
# Loop over Biodiversity Ecosystem Functions
par(mfrow = c(4,2), mar = c(2,7,2,1))
lapply(seq_along(BEF_cor), function(x){
  terra::rast((BEF_cor[[x]]**2)*sign(BEF_cor[[x]])) %>% # sqr scale for the colorbar
    terra::image(., col = curl_pal(100), zlim = c(-1,1), 
                 main = names(data_ecosystem[x]), axes = FALSE, legend = FALSE) # terra to have nice colorbars and shaping
  axis(side = 1, at = c(0,4,8,12,16,20)+0.5, labels = 0.5:5.5, line = -1, lwd = 0) # Hill numbers
  axis(side = 2, at = 0.5:3.5, labels = names(data_div) %>% rev(), las = 2, line = -0.5, lwd = 0) # Data type
  abline(v = 0:21, col = "white", lwd = 2) # White borders / gaps
  abline(h = 0:5, col = "white", lwd = 6) # White borders / gaps
  
}) # end lapply

# Fake plot for colorbars
z <- seq(-1,1,0.01) %>% as.matrix()
par(mar = c(4,8,4,4))
terra::image((z**2)*sign(z), col = curl_pal(100), axes = FALSE)
axis(side = 1, at = seq(0,1,0.1), labels = seq(-1,1,0.2), las = 1, lwd = 0, line = -0.5)

# --- 6. Biodiversity Ecosystem Function XY plots
# --- 6.1. Setup (INTERACTIVE)
id <- 5 # choose which ecosystem property
logged <- FALSE # choose if log or not
pal <- c(occurrence = "skyblue",abundance = "antiquewhite3", biomass = "gray20",metagenomics = "chocolate")

# --- 6.2. Ecosystem property
x <- data_ecosystem[[id]] %>% as.vector()
if(logged == TRUE){x <- log(x+1e-10)} # e-10 to avoid infinite values // more a plotting issue
q_ext <- quantile(x, c(0.025, 0.975), na.rm = TRUE)
x[x > q_ext[2]] <- q_ext[2]
x[x < q_ext[1]] <- q_ext[1]

# --- 6.3. Start plot
par(mfrow = c(1,1), mar = c(5,5,5,5))
plot(x = 1, y = 1, type = 'n', xlim = range(x, na.rm = TRUE), ylim = c(1,100), 
     axes = FALSE, xlab = names(data_ecosystem)[id], ylab = "Diversity hotspot probability (lines)",
     xaxs = "i", yaxs = "i")

# --- 6.4. Histogram of X
x_hist <- hist(x, breaks = 50, plot = FALSE)
x_hist$counts <- x_hist$counts / max(x_hist$counts)*100
plot(x_hist, add = TRUE, col = "gray95", border = "gray")

# --- 6.5. Fix axis
axis(side = 1, at = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 5),
     labels = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 5) %>% signif(2))
axis(side = 2, seq(0,100,20), las = 2)
axis(side = 4, seq(0, 100, 20), las = 2)
mtext("Normalized counts (histogram)", side = 4, line = 2)
abline(h = seq(0,100,20), lty = "dotted", col = "black")
box()

# --- 6.6. Plot by data type
# Loess by data type, superposed
lapply(seq_along(data_div), function(i){
  tmp <- data_div[[i]]
  lapply(1:dim(tmp)[[2]], function(j){
    good <- is.finite(x) & is.finite(tmp[, j]) # only complete case
    fit <- smooth.spline(x[good], tmp[good, j], spar = 1, cv = TRUE) # way faster than loess and no extrapolation
    pred <- predict(fit, seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100))$y
    
    lines(x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100), y = pred, col = scales::alpha(pal[i], 0.5), lwd = 1)
    if(j == 1){
      lines(x = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 100), y = pred, col = scales::alpha(pal[i], 0.5), lwd = 3) # thick line for richness
      text(x = max(x, na.rm = TRUE), y = tail(pred, 1),"q0", xpd = TRUE, pos = 4, col = pal[i])
      } # end if
  }) # end lapply hill
}) # end lapply data type



# END
