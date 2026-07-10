
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

# --- 2. Quality check success
# --- 2.1. Load quality checks success
load(file = paste0("./output/", FOLDER_NAME, "/TRADITIONAL_ABUNDANCE_success.RData"))
load(file = paste0("./output/", FOLDER_NAME, "/TRADITIONAL_BIOMASS_success.RData"))
load(file = paste0("./output/", FOLDER_NAME, "/METAGENOMICS_success.RData"))
load(file = paste0("./output/", FOLDER_NAME, "/OCCURRENCE_success.RData"))

# --- 2.2. Plot
# --- 2.2.1. Prepare
par(mfrow = c(1,1), mar = c(5,5,5,5))
pal <- scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6)
shift <- c(-0.3, -0.1, 0.1, 0.3)
success_list <- list(model_success_occurrence, model_success_abundance, model_success_biomass, model_success_metagenomic)

plot(1,1, xlim = c(3.74,100), ylim = c(0.5,6.5), type = 'n', axes = F, xlab = "Quality check success rate (%)", ylab = "Algorithm")

# --- 2.2.2. Add bars
lapply(seq_along(success_list), function(x){
  rect(xleft = 0, 
       ybottom = (1:6)-0.1+shift[x],
       xright = success_list[[x]][c(3,2,5,1,4,6)],
       ytop = (1:6)+0.1+shift[x], 
       col = pal[x])
  
}) # end lapply

# --- 2.2.3. Axes and esthetic
axis(side = 1, at = seq(0,100,20), labels = seq(0,100,20))
axis(side = 2, at = 1:6, labels = c("SVM","MLP","BRT","RF","GAM","GLM"), las = 2)
abline(v = seq(0,100,20), lty = "dotted", col = "gray20")
box()

# --- 2.3. Global average
mean(unlist(success_list))
sd(unlist(success_list))

# --- 3. MESS Analysis
# --- 3.1. Occurrence data
# --- 3.1.1. Which subfolder list and which model in the ensemble
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_OCCURRENCE_2025-06-12 16:07:51.160496"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# --- 3.1.2. Extract the vector of MESS values across all
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

# --- 3.2. Abundance and Biomass data
# --- 3.2.1. Observations are the same, we count double as there is two diversity from it
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_TRADITIONAL_2025-05-19 16:33:03.338843"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# --- 3.2.2. Extract the vector of MESS values across all
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

# --- 3.3. Metagenomic data
# --- 3.3.1. Which subfolder list and which model in the ensemble
CEPHALOPOD_OUTPUT <- "/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524"
all_files <- list.files(CEPHALOPOD_OUTPUT, recursive = TRUE)
model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]

# --- 3.3.2. Extract the vector of MESS values across all
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

# --- 3.4. Compute the average MESS
# --- 3.4.1. Prepare raster (set boundaries at 0 and -10 for colorbar)
mess_all <- abind(mess_occ, mess_trad, mess_trad, mess_omics, along = 2) %>% apply(1, median, na.rm = TRUE)
mess_all[mess_all > 0] <- 0
mess_all[mess_all < -10] <- -10
mess_rob <- r0 %>% setValues(mess_all) %>% terra::project(., robinson_proj)

# --- 3.4.2. Do the actual plot
plot(mess_rob, col = rocket_pal(11)[-11], axes = FALSE)
plot(land_rob, col = "gray20", add = TRUE, legend = FALSE)
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid
contour(mess_rob, add = TRUE, nlevels = 5)
box("figure", col="black", lwd = 1) # box

# --- 4. Global variance explained
# --- 4.1. Occurrence
load(paste0("./output/", FOLDER_NAME, "/OCCURRENCE_vip.RData"))
VIP <- vip %>% mutate(source = "occurrence")

# --- 4.2. Abundance
load(paste0("./output/", FOLDER_NAME, "/TRADITIONAL_ABUNDANCE_vip.RData"))
VIP <- bind_rows(VIP, vip_abundance %>% mutate(source = "abundance"))

# --- 4.3. Biomass
load(paste0("./output/", FOLDER_NAME, "/TRADITIONAL_BIOMASS_vip.RData"))
VIP <- bind_rows(VIP, vip_biomass %>% mutate(source = "biomass"))

# --- 4.4. Metagenomics
load(paste0("./output/", FOLDER_NAME, "/METAGENOMICS_vip.RData"))
VIP <- bind_rows(VIP, vip %>% mutate(source = "metagenomics"))

# --- 4.5. Cluster environmental variables
# At the global scale to ease interpretation on the maps
load("/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524/CALL.RData") # load omic' CALL
features <- CALL$ENV_DATA

# --- 4.5.1. Reshape as array cell * layer * month
features_array <- lapply(1:12, function(x)(x = features[[x]] %>% unwrap() %>% as.matrix())) %>% abind(along = 3)
dimnames(features_array) <- list(NULL, names(features[[1]] %>% unwrap()), as.character(1:12))

# --- 4.5.2. Yearly average
features_array_year <- apply(features_array, c(1,2), mean, na.rm = T)

# --- 4.5.3. Correlation clustering
library(dendextend)
c <- cor(features_array_year, method = "spearman", use = "pairwise.complete.obs") # Compute correlation matrix
d <- as.dist(1 - c) # Convert correlation to distance matrix (1 - correlation)
hc <- hclust(d, method = "complete") %>% cutree(4) # Hierarchical clustering
group <- data.frame(variable = names(hc), group = hc)

# --- 4.6. Prepare boxplot data
# We sum the variance explained per group and data source
df <- VIP %>% 
  left_join(group) %>% 
  group_by(source, group) %>% 
  summarize(value = sum(value)) %>% 
  ungroup()

# --- 4.7. Plot
# --- 4.7.1. Prepare
par(mfrow = c(1,1), mar = c(5,5,5,5))
pal <- scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6)
shift <- c(-0.3, -0.1, 0.1, 0.3)
plot_vect <- c("occurrence","abundance","biomass","metagenomics")

plot(1,1, xlim = c(3.74,60), ylim = c(0.5,4.5), type = 'n', axes = F, xlab = "Sum of variance explained (%)", ylab = "Environmental features cluster")

# --- 4.7.2. Add bars
lapply(seq_along(plot_vect), function(x){
  tmp <- df %>% dplyr::filter(source == plot_vect[x])
  rect(xleft = 0, 
       ybottom = (1:4)-0.1+shift[x],
       xright = tmp$value*100,
       ytop = (1:4)+0.1+shift[x], 
       col = pal[x])
})

# --- 4.7.3. Axes and esthetic
axis(side = 1, at = seq(0,100,20), labels = seq(0,100,20))
axis(side = 2, at = 1:4, labels = 1:4, las = 2)
abline(v = seq(0,100,20), lty = "dotted", col = "gray20")
box()


# --- 4.7. Plot normalized by group size
# --- 4.7.1. Prepare
par(mfrow = c(1,1), mar = c(5,5,5,5))
pal <- scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6)
shift <- c(-0.3, -0.1, 0.1, 0.3)
plot_vect <- c("occurrence","abundance","biomass","metagenomics")

plot(1,1, xlim = c(0.37,10), ylim = c(0.5,4.5), type = 'n', axes = F, xlab = "Mean of variance explained (%)", ylab = "Environmental features cluster")

# --- 4.7.2. Add bars
lapply(seq_along(plot_vect), function(x){
  tmp <- df %>% dplyr::filter(source == plot_vect[x])
  rect(xleft = 0, 
       ybottom = (1:4)-0.1+shift[x],
       xright = tmp$value*100 / table(group$group),
       ytop = (1:4)+0.1+shift[x], 
       col = pal[x])
})

# --- 4.7.3. Axes and esthetic
axis(side = 1, at = seq(0,10,2), labels = seq(0,10,2))
axis(side = 2, at = 1:4, labels = 1:4, las = 2)
abline(v = seq(0,10,2), lty = "dotted", col = "gray20")
box()
