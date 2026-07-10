
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

# --- 1.9. Moving average function
# Used only in this script, hence not in ./function
ma <- function(x, n = 10){if(length(x) > n) stats::filter(x, rep(1 / n, n), sides = 2) else NA}  # end function

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

# --- 3. Load ecosystem properties
data_ecosystem <- list()

# --- 3.1. Marine biodiversity Tittensor
s <- read_sf("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Tittensor&al._2010/01_Data/WCMC-019-PatternsBiodiversity2010-AcrossTaxa.shp")
r_interp <- akima::interp(x = s$X_COORD, y = s$Y_COORD, z = s$AllNorm,
                          xo = -180:179, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["Biodiversity"]] <- terra::rast(r_interp$z %>% t(), extent = c(-180,180,-90,90)) %>% flip() # build raster

# --- 3.2. Carbon dynamics // NPP
nc <- nc_open("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/DeVries&Weber_2017/Cexp_deVries_2017.nc")
r_interp <- akima::interp(x = ncvar_get(nc, "LON")[,,1], y = ncvar_get(nc, "LAT")[,,1], z = ncvar_get(nc, varid = "NPP") %>% apply(c(1,2), mean),
                          xo = 0:359, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["NPP"]] <- terra::rast(r_interp$z %>% t(), extent = c(0,360,-90,90)) %>% flip() %>% terra::rotate() # build raster

# --- 3.3. Carbon dynamics // FPOC
nc <- nc_open("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/DeVries&Weber_2017/Cexp_deVries_2017.nc")
z_ok <- ncvar_get(nc, varid = "FPOCex") %>% apply(c(1,2), mean)
z_ok[is.na(z_ok)] <- 0
r_interp <- akima::interp(x = ncvar_get(nc, "LON")[,,1], y = ncvar_get(nc, "LAT")[,,1], z_ok,
                          xo = 0:359, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
data_ecosystem[["FPOC"]] <- terra::rast(r_interp$z %>% t(), extent = c(0,360,-90,90)) %>% flip() %>% terra::rotate() # build raster

# --- 3.4. Carbon dynamics // efficiency
data_ecosystem[["Efficiency"]] <- data_ecosystem[["FPOC"]]/data_ecosystem[["NPP"]]

# --- 3.5. Plot ecosystem properties
par(mfrow = c(2,2))
lapply(seq_along(data_ecosystem), function(x){
  
  # --- 3.5.1. Prepare the rasters
  crs(data_ecosystem[[x]]) <- "+proj=longlat +datum=WGS84"
  r_rob <- terra::project(data_ecosystem[[x]], robinson_proj)
  
  # --- 3.5.2. Plot it
  plot(r_rob, col = viridis_pal(100), axes = FALSE, main = names(data_ecosystem)[x], 
       range = quantile(values(r_rob), c(0.025, 0.975), na.rm = TRUE), fill_range = TRUE)
  plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
  grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
    vect() %>%  project(robinson_proj) 
  plot(grat, lty = "dotted", add = TRUE) # add grid
  box("figure", col="black", lwd = 1) # box
})

# --- 4. Plot the emergent biodiversity ecosystem function
# --- 4.1. Preparation
pal <- c( m = "chocolate1", a = "#287E8C", b = "#9FD744", o = "#413078")
par(mfrow = c(2,2), mar = c(2,2,3,1))

# --- 4.2. Start the loop for plots
for (dt in 1:dim(data)[[5]]) {
  for (e in 1:length(data_ecosystem)) {
    # --- 4.2.1. Data preparation
    bio <- apply(data[,,,,dt], 1, mean, na.rm = TRUE) # biodiversity
    eco <- as.vector(data_ecosystem[[e]]) # ecosystem function
    ok <- is.finite(eco) & is.finite(bio)
    eco_ok <- eco[ok]
    bio_ok <- bio[ok]
    
    # --- 4.2.2. Fit polynomial
    df <- data.frame(eco_ok = eco_ok, bio_ok = bio_ok)
    fit <- lm(eco_ok ~ poly(bio_ok, 2, raw = T), data = df)
    
    # --- 4.2.3. Prediction grid (used for line and binning)
    x_pred <- seq(quantile(bio_ok, 0.05), quantile(bio_ok, 0.95), length.out = 100)
    y_pred <- predict(fit, newdata = data.frame(bio_ok = x_pred))
    
    # --- 4.2.4. Compute residuals for SD
    resid_df <- data.frame(resid = resid(fit), x = bio_ok)
    
    # --- 4.2.5. Bin residuals along the same x_pred grid
    x_bins <- cut(resid_df$x, breaks = x_pred, include.lowest = TRUE)
    bin_stats <- resid_df %>%
      group_by(x_bin = x_bins) %>%
      summarize(x_center = mean(x, na.rm = TRUE),
        sd_resid = sd(resid, na.rm = TRUE)) %>%
      arrange(x_center)
    
    # --- 4.2.6. Plot smoothScatter
    smoothScatter(x = bio_ok, y = eco_ok, nbin = 256,
      xlim = c(10,90), ylim = quantile(eco_ok, c(0.025, 0.975)),
      main = paste(dimnames(data)[[5]][dt], "*", names(data_ecosystem)[e]), cex.main = 0.6,
      colramp = colorRampPalette(c("white", pal[dt])), transformation = function(z) z^1, nrpoints = 0)
    
    # --- 4.2.7. Draw polynomial line
    lines(x_pred, y_pred, col = "white", lwd = 4)
    lines(x_pred, y_pred, col = "black", lwd = 2)
    
    # --- 4.2.8. Draw 1 SD dashed lines using bin_stats
    lines(bin_stats$x_center, y_pred[1:length(bin_stats$x_center)] + ma(bin_stats$sd_resid, 10),
          col = scales::alpha("black", 0.4), lwd = 1, lty = "longdash")
    lines(bin_stats$x_center, y_pred[1:length(bin_stats$x_center)] - ma(bin_stats$sd_resid, 10),
          col = scales::alpha("black", 0.4), lwd = 1, lty = "longdash")
    
    # --- 4.2.9. Print fit
    r2 <- summary(fit)$r.squared
    spearman <- cor(eco_ok, bio_ok, method = "spearman")
    rmse <- sqrt(mean(residuals(fit)^2))
    text(x = 90, y = quantile(eco_ok, 0.975)*0.95, paste("r-squarred:", signif(r2, 2)), pos = 2)
    text(x = 90, y = quantile(eco_ok, 0.975)*0.85, paste("spearman:", signif(spearman, 2)), pos = 2)
 
  } # data type loop
} # ecosystem function loop

# END
