
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

biomes <- terra::setValues(r0, as.vector(biomes)) %>% flip()
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
c <- cor(features_array_year, method = "spearman", use = "pairwise.complete.obs") # Compute correlation matrix
d <- as.dist(1 - c) # Convert correlation to distance matrix (1 - correlation)
hc <- hclust(d, method = "complete") %>% cutree(4) # Hierarchical clustering
group <- data.frame(variable = names(hc), group = hc)

# --- 3.6. Prepare boxplot data
# We sum the variance explained per group and data source
df <- VIP %>% 
  left_join(group) %>% 
  group_by(source, group) %>% 
  summarize(value = sum(value)) %>% 
  ungroup()

# --- 3.7. Plot
par(mfrow = c(2,2), mar = c(3, 15, 2, 2))
pal <- c(occurrence = "skyblue",abundance = "antiquewhite1",biomass = "antiquewhite3",metagenomics = "chocolate")
lapply(c("occurrence","abundance","biomass","metagenomics"), function(x){
  tmp <- df %>% dplyr::filter(source == x)
  plot(x = tmp$group, xlim = c(0.5,4.5), xlab = "", 
       y = tmp$value, ylim = c(0,1), ylab = "", 
       axes = FALSE, col = pal[x], pch = 20, cex = 3)
  abline(h = seq(0,1,0.2), lty = "dotted", col = "gray20")
  segments(x0 = tmp$group, y0 = 0, x1 = tmp$group, y1 = tmp$value)
  points(x = tmp$group, y = tmp$value, bg = pal[x], pch = 21, cex = 3)
  axis(side = 1, at = 1:4, labels = 1:4)
  axis(side = 2, las = 2)
  box()
})

# --- 4 Global average
par(mfrow = c(1,1))

# --- 4.1. Prepare the rasters
r <- r0
terra::values(r) <- abind(o,a,b,m, along = 2) %>% apply(1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# --- 4.2. Plot it
plot(r_rob, col = parula_pal(100), axes = FALSE, main = "global average")
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 4.3. Latitudinal gradiet test
r_df <- as.data.frame(r, xy = TRUE)
cor.test(r_df$y, r_df$lyr.1, method = "spearman")

# --- 5. Global bootstrap uncertainty
par(mfrow = c(1,1))

# --- 5.1. Prepare the rasters
r <- r0
terra::values(r) <- abind(o,a,b,m, along = 2) %>% apply(c(1,2,4), sd, na.rm = TRUE) %>% apply(1, mean, na.rm = TRUE) # global average
r_rob <- terra::project(r, robinson_proj)

# --- 5.2. Plot it
plot(r_rob, col = viridis_pal(100), axes = FALSE, main = "global SD", range = c(0,10), fill_range = T)
plot(land_rob, add = TRUE, legend = FALSE, axes = FALSE, col = "gray20")
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(r_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box


# --- 5. MESS
# Load all MESS analysis ever made...
# --- 5.1. Occurrence data
# Which subfolder list and which model in the ensemble
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_OCCURRENCE_2025-06-12 16:07:51.160496"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# Extract the vector of MESS values across all
mess_occ <- mclapply(model_files, function(x){
  memory_cleanup() # low memory use
  
  load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/MODEL.RData"))
  if(length(MODEL$MODEL_LIST) >= 1){
    load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/QUERY.RData"))
    out <- QUERY$MESS %>% unwrap() %>% median() %>% as.vector()
    return(out)
  } else {
    return(NULL)
  } # if model list
}, mc.cores = MAX_CLUSTERS) %>% 
  .[lengths(.) != 0] %>% 
  bind_cols() %>% 
  apply(1, median, na.rm = TRUE)

# --- 5.2. Abundance and Biomass data
# Observations are the same, we count double as there is two diversity from it
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_TRADITIONAL_2025-05-19 16:33:03.338843"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# Extract the vector of MESS values across all
mess_trad <- mclapply(model_files, function(x){
  memory_cleanup() # low memory use
  
  load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/MODEL.RData"))
  if(length(MODEL$MODEL_LIST) >= 1){
    load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/QUERY.RData"))
    out <- QUERY$MESS %>% unwrap() %>% median() %>% as.vector()
    return(out)
  } else {
    return(NULL)
  } # if model list
}, mc.cores = MAX_CLUSTERS) %>% 
  .[lengths(.) != 0] %>% 
  bind_cols() %>% 
  apply(1, median, na.rm = TRUE)

# --- 5.3. Metagenomic data
# Which subfolder list and which model in the ensemble
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# Extract the vector of MESS values across all
mess_omics <- mclapply(model_files, function(x){
  memory_cleanup() # low memory use
  
  load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/MODEL.RData"))
  if(length(MODEL$MODEL_LIST) >= 1){
    load(paste0(CEPHALOPOD_OUTPUT,"/", x, "/QUERY.RData"))
    out <- QUERY$MESS %>% unwrap() %>% median() %>% as.vector()
    return(out)
  } else {
    return(NULL)
  } # if model list
}, mc.cores = MAX_CLUSTERS) %>% 
  .[lengths(.) != 0] %>% 
  bind_cols() %>% 
  apply(1, median, na.rm = TRUE)

# --- 5.4. Compute the average MESS
# Prepare raster (set boundaries at 0 and -10 for colorbar)
mess_all <- abind(mess_occ, mess_trad, mess_trad, mess_omics, along = 2) %>% apply(1, median, na.rm = TRUE)
mess_all[mess_all > 0] <- 0
mess_all[mess_all < -10] <- -10
mess_rob <- r0 %>% setValues(mess_all) %>% terra::project(., robinson_proj)

# Do the actual plot
plot(mess_rob, col = rocket_pal(11)[-11], axes = FALSE)
plot(land_rob, col = "gray20", add = TRUE, legend = FALSE)
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(mess_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 6. PERMANANOVA on axis of variance
data <- abind(o,a,b,m, along = 2) %>% apply(c(1,2,4), mean, na.rm = TRUE)
id <- which(is.na(data[,1,1]))

library(vegan)
library(coop)

# Collapse into 2D: rows = (hill × month), cols = cells
mat <- apply(data[-id,,], c(1), function(x){x = c(x)}) 

# Build metadata
meta <- expand.grid(
  hill = dimnames(data)[[2]],
  month = dimnames(data)[[3]]
)
meta$hill  <- factor(meta$hill)
meta$month <- factor(meta$month)
meta$type <- c(rep("occurrence", dim(o)[[2]]*dim(o)[[4]]),
               rep("abundance", dim(a)[[2]]*dim(a)[[4]]),
               rep("biomass", dim(b)[[2]]*dim(b)[[4]]),
               rep("metagenomic", dim(m)[[2]]*dim(m)[[4]]))
meta$type <- as.factor(meta$type)

dist_mat <- dist(mat)

# PERMANOVA
adonis_res <- adonis2(dist_mat ~ hill + month + type, data = meta, method = "euclidian", na.rm = TRUE, parallel = 20, by = "onedf")
print(adonis_res)

# Extract R² and turn into percentages
var_exp <- 100 * adonis_res$R2
var_exp





# END