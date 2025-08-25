
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
raster_considered <- lapply(c(export = "/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Henson&al._2012/export_Henson2012.nc",
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

# --- 3.2. Fisheries
# Pelagic data from Watson 2017
load("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Fisheries_Watson2017/catch_data/clim_all_pelagics_catches_1990-2019.Rdata")
raster_considered[["pelagic_catch"]] <- terra::rasterize(x = clim_all_pel[,2:3] %>% as.matrix(),
                                                         y = r0,
                                                         values = clim_all_pel[,7])

# Demersal data from Watson 2017
load("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Fisheries_Watson2017/catch_data/clim_demersal_catches_1990-2019.Rdata")
raster_considered[["demersal_catch"]] <- terra::rasterize(x = clim_dem[,3:4] %>% as.matrix(),
                                                          y = r0,
                                                          values = clim_dem[,8])

# --- 3.3. Particle size spectrum
# Slope
load("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Clements&al.2023/table_clim_mon_part_BV_slope_carbon_Clements&al.2023_07_08_23.RData")
tab <- tab %>% dplyr::filter(variable == "Slope")
raster_considered[["slope"]] <- terra::rasterize(x = tab[,3:4] %>% as.matrix(),
                                                          y = r0,
                                                          values = tab[,5])

# Biovolume
load("/net/sea/work/fabioben/OVERSEE/data/env_predictors/Global_ecosystem_properties/Carbon/Clements&al.2023/table_clim_mon_part_BV_slope_carbon_Clements&al.2023_07_08_23.RData")
tab <- tab %>% dplyr::filter(variable == "Biovolume")
raster_considered[["biovolume"]] <- terra::rasterize(x = tab[,3:4] %>% as.matrix(),
                                                 y = r0,
                                                 values = tab[,5])

# --- 4. Extract Hill maps
# Geographical cell x hill x data type - loess plot for smoothness
to_map <- list(occurrence = o, abundance = a, biomass = b, metagenomic = m)
hill_maps <- lapply(seq_along(to_map), function(x){
  tmp <- apply(to_map[[x]], c(1,2), mean, na.rm = TRUE)
  })

# Compute the global average
tmp <- abind(o,a,b,m, along = 2) %>% apply(c(1,2), mean, na.rm = TRUE) # global average (summer / winter)

# Average dimnames with the same name
all <- sapply(unique(dimnames(tmp)[[2]]), function(n) {
  id <- which(dimnames(tmp)[[2]] == n)
  if (length(id) == 1) {tmp[, id]} else {apply(tmp[, id], 1, mean, na.rm = TRUE)}
}) %>% as.matrix()

hill_maps[[5]] <- all # add 5th dimension as an average

plot_data <- lapply(seq_along(raster_considered), function(x){
  # Extract values
  val <- raster_considered[[x]] %>% as.vector() # extract ecosystem function
  
  lapply(seq_along(hill_maps), function(i){
      # Do correlation
      corval <- cor(val, hill_maps[[i]], use = "pairwise.complete.obs", method = "spearman") %>% abs()
      
      # Do loess
      fit <- loess(corval[1,] ~ str_sub(colnames(corval), 6, -1), span = 1) # Fit loess smoother
      pred <- predict(fit, newdata = seq(0, 5, 0.25)) # Predict
      return(pred)

  }) %>% bind_cols() %>% t() # data type loop
}) # ecosystem function loop

# Plot
par(mfrow = c(4,2), mar = c(1,1,3,1))
lapply(seq_along(plot_data), function(x){
  plot(plot_data[[x]] %>% rast(), col = viridis_pal(100), range = c(0.5,1), main = names(raster_considered)[x])
  abline(h = c(0:5), v = c(0:21), col = "white", lwd = 3)
})


# END
