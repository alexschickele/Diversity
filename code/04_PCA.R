
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

# --- 2. Load data
# --- 2.1. Metagenomics
load("/net/meso/work/aschickele/Diversity/output/DIVERSITY_PAPER/METAGENOMICS_raw_diversity.RData")
tmp <- apply(data, -1, quantile_scale) # rescale
data <- abind(tmp[1:32400,,,], tmp[32401:64800,,,c(7:12,1:6)], along = 1) # southern hemisphere swap

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

o <- apply(data, 1, c)  # Collapse
metadata <- expand.grid(method = paste("Hill", dimnames(data)[[2]]),
                        bootstrap = dimnames(data)[[3]],
                        month = dimnames(data)[[4]],
                        data_type = "OCCURRENCE")
o <- cbind(o, metadata) # rows = diversity estimates, columns = cells + factors
rm(data)
gc()

# --- 2.5. Concatenate
# Bind all individual diversity estimates as rows
data <- rbind(m, a, b, o)
rm(m, a, b, o, metadata)
gc()

# Add column names for the cells
colnames(data) <- c(paste0("C", 1:64800), "METHOD","BOOTSTRAP","MONTH","DATA_TYPE")

# Remove NA columns
id <- which(is.na(apply(data[, 1:64800], 2, sum)))
data <- data[,-id]

# --- 2.5. Split PCA and RDA data
metadata <- data[, (ncol(data)-3):ncol(data)]
data <- data[, 1:(ncol(data)-4)]

# --- 3. PCA on diversity estimates
# --- 3.1. Perform PCA // Species = diversity ; sites = geographical cells
PCA <- vegan::pca(X = t(data), scale = T)

# --- 3.2. Extract contribution of each Principal Component (PC)
PC_percent <- (eigenvals(PCA)/sum(eigenvals(PCA))) # variance per PC

# --- 3.3. Extract coordinates
variable_scores <- vegan::scores(PCA, choices = c(1,2), display = "species") %>% 
  cbind(metadata) %>% # association
  group_by(METHOD, MONTH, DATA_TYPE) %>% 
  summarize(PC1 = mean(PC1),
            PC2 = mean(PC2)) %>% 
  ungroup()
  
centroids <- lapply(c("METHOD","MONTH","DATA_TYPE"), function(x){
  variable_scores %>% 
    group_by_at(x) %>% 
    summarise(X = mean(PC1), Y = mean(PC2)) %>% 
    ungroup()
}) # end lapply

# --- 3.4. Estimate contribution
# As the squared scores across PC1 and PC2
contribution <- variable_scores %>% 
  mutate(contribution = (PC1**2) + (PC2**2)) %>% 
  dplyr::select(-PC1, -PC2) %>% 
  mutate(contribution_perc = contribution / sum(contribution))

# Now we compute the average score per factor
contribution <- lapply(c("METHOD","MONTH","DATA_TYPE"), function(x){
  tmp <- contribution %>% 
    group_by_at(x) %>% 
    summarise(FACTOR = x,
              contribution_perc = mean(contribution_perc)) %>% 
    ungroup()
  colnames(tmp) <- c("MODALITY","FACTOR","contribution_perc")
  return(tmp)
}) %>% bind_rows() # end lapply

# Relative to the mean to be able to compare
contribution <- contribution %>% 
  mutate(contribution_perc = (contribution_perc / mean(contribution_perc) -1) %>% abs()) %>% 
  arrange(contribution_perc)

# --- 4. Graphical output
# --- 4.1. Set up
# pal <- mako_pal(6)[2:5]
# pal <- brewer.pal(4, "Set2")
pal <- c("red3", "orange2", "antiquewhite4") # wes anderson

# --- 4.1. Contribution plot
par(mar = c(8,3,2,2))
plot(x = 1:nrow(contribution), y = contribution$contribution_perc %>% rev(),
     xlab = "", ylab = "PCA Scores rel. to the mean", axes = FALSE, type = "n")
abline(v = 1:nrow(contribution), lwd = 25, col = c("white","gray95"))
segments(x0 = 1:nrow(contribution), y0 = 0, y1 = contribution$contribution_perc %>% rev(), lwd = 2)
axis(side = 2, las = 2)
axis(side = 1, las = 2, at = 1:nrow(contribution), labels = contribution$MODALITY %>% rev(), cex.axis = 0.6)
abline(h = seq(-0.4, 0.2, 0.1), lty = c("dashed","dashed","dashed","dashed","solid","dashed","dashed"),
       col = c("gray","gray","gray","gray","black","gray","gray"))
points(x = 1:nrow(contribution), y = contribution$contribution_perc %>% rev(),
       pch = 21, lwd = 2, cex = 3, bg = pal[as.factor(contribution$FACTOR %>% rev())])
box()

# --- 4.2. PCA
# --- 4.2.1. By data type
pal <- c("chocolate2", "#CCC591", "darkolivegreen4", "antiquewhite4") # wes anderson
# pal <- inferno_pal(5)[1:4]

# We first perform the point plot
par(mar = c(3,3,3,7))
plot(variable_scores$PC1, variable_scores$PC2, type = "n", axes = FALSE, xlim = c(-0.8, 1.6), ylim = c(-0.8, 1.2), xlab = "", ylab = "")
points(variable_scores$PC1, variable_scores$PC2,
     pch = 16, lwd = 1, col = pal[as.factor(variable_scores$DATA_TYPE)] %>% scales::alpha(0.3))
abline(h = 0, v = 0, lty = "dashed")

mtext("0", side = 1, line = 1, at = 0)
mtext("0", side = 2, line = 1, at = 0, las = 2)
mtext(paste("PC1 (", round(PC_percent[1], 2)*100,"% )"), side = 4, line = 1, at = 0, las = 2)
mtext(paste("PC2 (", round(PC_percent[2], 2)*100,"% )"), side = 3, line = 1, at = 0)

# Add the ellipses
lapply(1:length(unique(variable_scores$DATA_TYPE)), function(x){
  df <- variable_scores %>% 
    dplyr::select(PC1, PC2, DATA_TYPE) %>% 
    dplyr::filter(DATA_TYPE == unique(variable_scores$DATA_TYPE)[x]) %>% 
    dplyr::select(-DATA_TYPE) 
  lines(ellipse(cov(df), centre = colMeans(df)), col = pal[x], lwd = 2)
}) # end ellipse loop

# --- 4.2.2. By Hill number
# We skip ellipses here as they are too many
pal <- rocket_pal(21)

par(mar = c(3,3,3,7))
plot(variable_scores$PC1, variable_scores$PC2, type = "n", axes = FALSE, xlim = c(-0.8, 1.6), ylim = c(-0.8, 1.2), xlab = "", ylab = "")
points(variable_scores$PC1, variable_scores$PC2,
       pch = 16, lwd = 1, col = pal[as.factor(variable_scores$METHOD)] %>% scales::alpha(0.3))
abline(h = 0, v = 0, lty = "dashed")

mtext("0", side = 1, line = 1, at = 0)
mtext("0", side = 2, line = 1, at = 0, las = 2)
mtext(paste("PC1 (", round(PC_percent[1], 2)*100,"% )"), side = 4, line = 1, at = 0, las = 2)
mtext(paste("PC2 (", round(PC_percent[2], 2)*100,"% )"), side = 3, line = 1, at = 0)

# --- 4.2.3. By Month
# We skip ellipses here as they are too many
pal <- circular_pal(12)

par(mar = c(3,3,3,7))
plot(variable_scores$PC1, variable_scores$PC2, type = "n", axes = FALSE, xlim = c(-0.8, 1.6), ylim = c(-0.8, 1.2), xlab = "", ylab = "")
points(variable_scores$PC1, variable_scores$PC2,
       pch = 16, lwd = 1, col = pal[as.factor(variable_scores$MONTH)] %>% scales::alpha(0.3))
abline(h = 0, v = 0, lty = "dashed")

mtext("0", side = 1, line = 1, at = 0)
mtext("0", side = 2, line = 1, at = 0, las = 2)
mtext(paste("PC1 (", round(PC_percent[1], 2)*100,"% )"), side = 4, line = 1, at = 0, las = 2)
mtext(paste("PC2 (", round(PC_percent[2], 2)*100,"% )"), side = 3, line = 1, at = 0)

# --- 4.3. Principal component spatial pattern
# --- 4.3.1. Compute site scores and default cell vector
site_scores <- vegan::scores(PCA, choices = c(1,2), display = "site") %>% as.data.frame()
colnames(site_scores) <- c("PC1", "PC2")
val <- rep(NA, 64800)

# --- 4.3.1. Map of PC1
val[-id] <- site_scores$PC1 # assign values
r <- setValues(r0, val) %>% terra::project(robinson_proj) # project

plot(r, col = c(inferno_pal(100), "white") %>% rev(), range = c(0,1.5), 
     axes = FALSE, fill_range = TRUE, main = "PC1") # plot
plot(land_rob, col = "gray20", add = TRUE, legend = FALSE)
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, add = TRUE, lty = "dotted") # add grid

# --- 4.3.1. Map of PC2
val[-id] <- site_scores$PC2 # assign values
r <- setValues(r0, val) %>% terra::project(robinson_proj) # project

plot(r, col = c(inferno_pal(100), "white") %>% rev(), range = c(0,1.5), 
     axes = FALSE, fill_range = TRUE, main = "PC2") # plot
plot(land_rob, col = "gray20", add = TRUE, legend = FALSE)
grat <- sf::st_graticule(lon = c(seq(-180,180, 30), -31), lat = c(seq(-90,90, 30), 89)) %>%
  vect() %>%  project(robinson_proj) 
plot(grat, add = TRUE) # add grid




#rda# END
