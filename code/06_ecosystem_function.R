
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

# --- 3. Load ecosystem function
# --- 3.1. Henson 2012 carbon
henson_rasters <- lapply(c(export = "/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Henson&al._2012/export_Henson2012.nc",
                           PEeff = "/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Henson&al._2012/PEeff_Henson2012.nc",
                           poc2000 = "/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Henson&al._2012/poc2000_Henson2012.nc",
                           Teff = "/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Henson&al._2012/Teff_Henson2012.nc"),
       function(INPUT_FILE){
         
         library(ncdf4)
         library(akima)
         
         nc <- nc_open(INPUT_FILE)
         ncvar <- ncvar_get(nc) # assuming the first one
         
         lon <- nc[["dim"]][["lon"]][["vals"]] %% 360 # modulo 360 trick
         lat <- nc[["dim"]][["lat"]][["vals"]]
         id <- which(is.infinite(ncvar) | is.na(ncvar))
         
         r_interp <- akima::interp(x = rep(lon, 200)[-id], y = rep(lat, each = 360)[-id], z = as.vector(ncvar)[-id],
                                   xo = 0:359, yo = -89.5:89.5, linear = T, extrap = F) # interpolate
         
         r <- terra::rast(r_interp$z %>% t(), extent = c(0,360,-90,90)) %>% flip() %>% rotate() # build raster
         message(paste("--- DONE:", INPUT_FILE))
         return(r)
         
       }) # end lapply

# --- 4. Extract Hill maps
# Geographical cell x hill x data type - loess plot for smoothness
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
hill_maps <- lapply(seq_along(to_map), function(x){
  tmp <- apply(to_map[[x]], c(1,2), mean, na.rm = TRUE)})

# Compute the global average
tmp <- abind(o,a,b,m, along = 2) %>% apply(c(1,2), mean, na.rm = TRUE) # global average

# Average dimnames with the same name
all <- sapply(unique(dimnames(tmp)[[2]]), function(n) {
  id <- which(dimnames(tmp)[[2]] == n)
  if (length(id) == 1) {tmp[, id]} else {apply(tmp[, id], 1, mean, na.rm = TRUE)}
}) %>% as.matrix()

hill_maps[[5]] <- all # add 5th dimension as an average

par(mfrow = c(1,4), mar = c(1,1,3,1))

lapply(seq_along(henson_rasters), function(x){
  # Extract values
  val <- henson_rasters[[x]] %>% as.vector() # extract ecosystem function
  
  pal <- c("antiquewhite4", "#CCC591", "darkolivegreen3", "chocolate2","black") # wes anderson
  plot(x = 1, y = 1, xlim = c(0,5), ylim = c(-1, 1), type = "n")
  
  lapply(seq_along(hill_maps), function(i){
    
    # Do correlation
    corval <- cor(val, hill_maps[[i]], use = "pairwise.complete.obs", method = "spearman")
    
    # Do loess
    fit <- loess(corval[1,] ~ str_sub(colnames(corval), 6, -1), span = 1) # Fit loess smoother
    pred <- predict(fit, newdata = seq(0, 5, 0.25)) # Predict
    lines(seq(0, 5, 0.25), pred, col = pal[i], lwd = 3) # Plot
    
    }) # data type loop
  
  abline(h = seq(-1, 1, 0.2), lty = "dotted")
  
  
}) # ecossytem function loop



# END
