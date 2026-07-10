# Plot the status of the data

# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name
FOLDER_NAME = "DIVERSITY_PAPER"

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

rm(land, land_raster)
gc()

# --- 1.7. Quantile scale
# To make sure its correct across all data types
quantile_scale <- function(x){
  quantiles <- quantile(x, probs = seq(0, 1, by = 0.01), na.rm = TRUE) %>% unique()  # Change the 'probs' argument as needed
  quantile_values <- cut(x, breaks = quantiles, include.lowest = TRUE, labels = FALSE)
  return(quantile_values)
} # function

# --- 2. Load all datasets
# --- 2.1. Raw datasets
# just a name repair on the column type
occ_data <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/OCCURRENCE_norm_input.csv")) %>% 
  mutate(worms_id = as.character(worms_id),
         measurementvalue = 1,
         measurementunit = "Presence")
trad_data <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_norm_input.csv")) %>% 
  mutate(worms_id = as.character(worms_id))
omic_data <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/METAGENOMIC_norm_input.csv")) %>% 
  mutate(worms_id = as.character(worms_id))

# --- 2.2. Concatenate
# Put all together in the same data frame ; regrid at 1 x 1 degree for FAST
res <- 1 # geographical resolution
digit <- 0 # digits to keep
all_data <- bind_rows(occ_data, trad_data, omic_data) %>% 
  dplyr::filter(measurementunit %in% c("Presence", "ind m-3", "mgC m-3", "relative metagenomic reads")) %>% 
  mutate(decimallatitude = round(decimallatitude + 0.5 * res, digits = digit) - 0.5 * res,
         decimallongitude = round(decimallongitude + 0.5 * res, digits = digit) - 0.5 * res) %>% 
  dplyr::select(-year, -depth, -month) %>% 
  distinct()

# --- 3. Plot map data
# --- 3.1. Process samples as spatial points
samples_sf <- sf::st_as_sf(all_data, coords = c("decimallongitude", "decimallatitude"), crs = 4326)
samples_proj <- sf::st_transform(samples_sf, crs = robinson_proj)
type_colors <- scales::alpha(c("ind m-3" = "#287E8C", "mgC m-3" = "#9FD744", "Presence" = "#413078", "relative metagenomic reads" = "chocolate1"), c(0.05, 0.05, 0.2, 1))
type_cex <- c("ind m-3" = 0.3, "mgC m-3" = 0.3, "Presence" = 0.6, "relative metagenomic reads" = 1)

# --- 3.2. Do the plot
plot(land_rob, col = "gray20", box = FALSE, axes = FALSE)
points(sf::st_coordinates(samples_proj), col = type_colors[all_data$measurementunit], pch = 15, cex = type_cex[all_data$measurementunit]) # add points
plot(land_rob, col = "gray20", box = FALSE, axes = FALSE, add = TRUE)
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), 29), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, lty = "dotted", add = TRUE) # add grid

# --- 4. Number of species per latitude bin ; hypothesis of the LDG
# --- 4.0. Graphical parameters
par(mfrow = c(1, 4), mar = c(3, 3, 2, 1))

# --- 4.1. Occurrence data
tmp <- all_data %>% 
  dplyr::filter(measurementunit == "Presence") %>% 
  dplyr::select(decimallatitude, decimallongitude, scientificname) %>% 
  group_by(decimallatitude, decimallongitude) %>% 
  summarise(measurementvalue = n()) %>% 
  ungroup() %>% 
  mutate(measurementvalue = quantile_scale(measurementvalue),
         measurementunit = "Presence") %>% 
  mutate(latbin = cut(decimallatitude, breaks = seq(-90, 90, 10)),
         latid = as.numeric(latbin)) 

boxplot(tmp$measurementvalue ~ tmp$latid, horizontal = T, axes = FALSE, 
        col = scales::alpha("#413078", 0.5), outline = F)
axis(side = 1, at = seq(0, 100, 20), labels = seq(0, 100, 20), cex.axis = 0.7)
axis(side = 2, at = 0.5:18.5, labels = seq(-90, 90, 10), las = 2, cex.axis = 0.7)
box()

# Spearman correlation
cor.test(tmp$measurementvalue, abs(tmp$decimallatitude), method = "spearman")

# --- 4.2. Abundance data
tmp <- all_data %>% 
  dplyr::filter(measurementunit == "ind m-3",
                scientificname == "Hill 0 ( ind m-3 )") %>% 
  dplyr::select(decimallatitude, decimallongitude, measurementvalue) %>% 
  mutate(latbin = cut(decimallatitude, breaks = seq(-90, 90, 10)),
         latid = as.numeric(latbin)) 

boxplot(tmp$measurementvalue ~ tmp$latid, horizontal = T, axes = FALSE, 
        col = scales::alpha("#287E8C", 0.5), outline = F)
axis(side = 1, at = seq(0, 100, 20), labels = seq(0, 100, 20), cex.axis = 0.7)
axis(side = 2, at = 0.5:18.5, labels = seq(-90, 90, 10), las = 2, cex.axis = 0.7)
box()

# Spearman correlation
cor.test(tmp$measurementvalue, abs(tmp$decimallatitude), method = "spearman")

# --- 4.3. Biomass data
tmp <- all_data %>% 
  dplyr::filter(measurementunit == "mgC m-3",
                scientificname == "Hill 0 ( mgC m-3 )") %>% 
  dplyr::select(decimallatitude, decimallongitude, measurementvalue) %>% 
  mutate(latbin = cut(decimallatitude, breaks = seq(-90, 90, 10)),
         latid = as.numeric(latbin)) 

boxplot(tmp$measurementvalue ~ tmp$latid, horizontal = T, axes = FALSE, 
        col = scales::alpha("#9FD744", 0.5), outline = F)
axis(side = 1, at = seq(0, 100, 20), labels = seq(0, 100, 20), cex.axis = 0.7)
axis(side = 2, at = 0.5:18.5, labels = seq(-90, 90, 10), las = 2, cex.axis = 0.7)
box()

# Spearman correlation
cor.test(tmp$measurementvalue, abs(tmp$decimallatitude), method = "spearman")

# --- 4.4. Genomic data
tmp <- all_data %>% 
  dplyr::filter(measurementunit == "relative metagenomic reads", 
                scientificname == "Hill 0 ( relative metagenomic reads )") %>% 
  dplyr::select(decimallatitude, decimallongitude, measurementvalue) %>% 
  mutate(latbin = cut(decimallatitude, breaks = seq(-90, 90, 10)),
         latid = as.numeric(latbin)) 

boxplot(tmp$measurementvalue ~ tmp$latid, horizontal = T, axes = FALSE, 
        col = scales::alpha("chocolate1", 0.5), outline = F)
axis(side = 1, at = seq(0, 100, 20), labels = seq(0, 100, 20), cex.axis = 0.7)
axis(side = 2, at = 0.5:18.5, labels = seq(-90, 90, 10), las = 2, cex.axis = 0.7)
box()

# Spearman correlation
cor.test(tmp$measurementvalue, abs(tmp$decimallatitude), method = "spearman")

# --- 5. Taxonomy
# --- 5.1. Open datasets
# --- 5.1.1. Metagenomic
omic_taxo <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/METAGENOMIC_taxo.csv")) %>% 
  mutate(value = reads/sum(reads), n = n/sum(n)) %>% arrange(n)
tmp <- omic_taxo %>% slice_head(n = -5) %>% 
  summarise(class = "Others", value = sum(value), n = sum(n))
omic_taxo <- omic_taxo %>% slice_tail(n = 5) %>% bind_rows(tmp)
cor.test(omic_taxo$value, omic_taxo$n, method = "spearman")

# --- 5.1.2. Abundance
abundance_taxo <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_ABUNDANCE_taxo.csv")) %>% 
  mutate(value_a = measurementvalue/sum(measurementvalue), n_a = n/sum(n)) %>% arrange(n_a)
tmp <- abundance_taxo %>% slice_head(n = -5) %>% 
  summarise(class = "Others", value_a = sum(value_a), n_a = sum(n_a))
abundance_taxo <- abundance_taxo %>% slice_tail(n = 5) %>% bind_rows(tmp) %>% dplyr::select(class, n_a, value_a)
cor.test(abundance_taxo$value_a, abundance_taxo$n_a, method = "spearman")

# --- 5.1.3. Biomass
biomass_taxo <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/TRADITIONAL_BIOMASS_taxo.csv")) %>% 
  mutate(value_b = measurementvalue/sum(measurementvalue), n_b = n/sum(n)) %>% arrange(n_b)
tmp <- biomass_taxo %>% slice_head(n = -5) %>% 
  summarise(class = "Others", value_b = sum(value_b), n_b = sum(n_b))
biomass_taxo <- biomass_taxo %>% slice_tail(n = 5) %>% bind_rows(tmp) %>% dplyr::select(class, n_b, value_b)
cor.test(biomass_taxo$value_b, biomass_taxo$n_b, method = "spearman")

# --- 5.1.4. Occurrence
occ_taxo <- vroom(paste0("/nfs/meso/work/aschickele/Diversity/output/",FOLDER_NAME,"/OCCURRENCE_taxo.csv"), delim = ",") %>% 
  mutate(n_o = n/sum(n)) %>% arrange(n_o)
tmp <- occ_taxo %>% slice_head(n = -5) %>% 
  summarise(class = "Others", n_o = sum(n_o))
occ_taxo <- occ_taxo %>% slice_tail(n = 5) %>% bind_rows(tmp) %>% dplyr::select(class, n_o)

# --- 5.2. Setup
library(plotrix)
par(mfrow = c(1,2))

# --- 5.3. Omic plot separately
# Does not share taxonomy with other groups
pal <- c(inferno_pal(5), "white")
plot(0, 0, xlim = c(-1,1), ylim = c(-1,1), type = "n", main = "Metagenomics", axes = F, xlab = "", ylab = "")
out_pie <- floating.pie(0, 0, x = omic_taxo$value, col = pal, radius = 0.9, border = "white") # outer pie
pie.labels(0, 0, out_pie, labels = omic_taxo$value %>% round(1), radius = 0.7, col = "black", cex = 1)
in_pie <- floating.pie(0, 0, x = omic_taxo$n, radius = 0.5, col = pal, border = "white") # inner pie
pie.labels(0, 0, in_pie, labels = omic_taxo$n %>% round(1), radius = 0.25, col = "black", cex = 1)
points(0, 0, pch = 1, col = "white", cex = 22, lwd = 20) # inner border

# Print colors for paper figure
print(pal)
print(omic_taxo$class)

# --- 5.4. Others
# Common colors as they share taxonomy
joined_taxo <- occ_taxo %>%
  dplyr::select(class, n_o) %>% 
  full_join(abundance_taxo, by = "class") %>%
  full_join(biomass_taxo, by = "class")
joined_taxo <- rbind(joined_taxo[-6,], joined_taxo[6,])
pal <- c(viridis_pal(6) %>% rev(), "white")

# Abundance
plot(0, 0, xlim = c(-1,1), ylim = c(-1,1), type = "n", main = "Abundance", axes = F, xlab = "", ylab = "")
out_pie <- floating.pie(0, 0, x = joined_taxo$value_a, col = pal, radius = 0.9, border = "white") # outer pie
pie.labels(0, 0, out_pie, labels = joined_taxo$value_a %>% round(1) %>% .[!is.na(.)], radius = 0.7, col = "black", cex = 1)
in_pie <- floating.pie(0, 0, x = joined_taxo$n_a, radius = 0.5, col = pal, border = "white") # inner pie
pie.labels(0, 0, in_pie, labels = joined_taxo$n_a %>% round(1) %>% .[!is.na(.)], radius = 0.25, col = "black", cex = 1)
points(0, 0, pch = 1, col = "white", cex = 22, lwd = 20) # inner border

# Biomass
plot(0, 0, xlim = c(-1,1), ylim = c(-1,1), type = "n", main = "Biomass", axes = F, xlab = "", ylab = "")
out_pie <- floating.pie(0, 0, x = joined_taxo$value_b, col = pal, radius = 0.9, border = "white") # outer pie
pie.labels(0, 0, out_pie, labels = joined_taxo$value_b %>% round(1) %>% .[!is.na(.)], radius = 0.7, col = "black", cex = 1)
in_pie <- floating.pie(0, 0, x = joined_taxo$n_b, radius = 0.5, col = pal, border = "white") # inner pie
pie.labels(0, 0, in_pie, labels = joined_taxo$n_b %>% round(1) %>% .[!is.na(.)], radius = 0.25, col = "black", cex = 1)
points(0, 0, pch = 1, col = "white", cex = 22, lwd = 20) # inner border

# Occurrence
plot(0, 0, xlim = c(-1,1), ylim = c(-1,1), type = "n", main = "Occurrence", axes = F, xlab = "", ylab = "")
in_pie <- floating.pie(0, 0, x = joined_taxo$n_o, radius = 0.5, col = pal, border = "white") # inner pie
pie.labels(0, 0, in_pie, labels = joined_taxo$n_o %>% round(1) %>% .[!is.na(.)], radius = 0.25, col = "black", cex = 1)
points(0, 0, pch = 1, col = "white", cex = 22, lwd = 20) # inner border

print(joined_taxo)
print(pal)

# END

