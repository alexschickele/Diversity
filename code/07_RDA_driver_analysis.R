# This script contains the driver analysis associated to the diversity modelled
# by CEPHALOPOD. It corresponds to the Figure 3 and associated supplementary
# materials.

# --- 1. Initialize
# --- 1.1. Set correct working directory
setwd("/net/meso/work/aschickele/Diversity")

# --- 1.2. Set folder name and Hill numbers to test
FOLDER_NAME <- "DIVERSITY_PAPER"

# --- 1.3. Source all libraries & functions
source(file = "./code/00_config.R")

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

# --- 2.6. Collapse month into first dimension
data_long <- aperm(data, c(1,4,2,3,5))
data_long <- array(data_long, dim = c(64800*12, 21*10*4))
data_names_long <- expand_grid(dimnames(data)[[5]],dimnames(data)[[3]],dimnames(data)[[2]])
dimnames(data_names_long)[[2]] <- c("Type","Bootstrap","Hill")

# --- 3. Environmental features
# --- 3.1. All features
load("/net/meso/work/aschickele/CEPHALOPOD/output/DIVERSITY_MOTU_RAREFIED_2025-06-17 14:21:26.693524/CALL.RData") # load omic' CALL
features <- CALL$ENV_DATA

# --- 3.2. Reshape as array cell * layer * month
features_array <- lapply(1:12, function(x)(x = features[[x]] %>% unwrap() %>% as.matrix())) %>% abind(along = 3)
dimnames(features_array) <- list(NULL, names(features[[1]] %>% unwrap()), as.character(1:12))
features_array[, "climatology_fsle_aviso_2001_2020", ] <- features_array[, "climatology_fsle_aviso_2001_2020", ] * (-1)

# --- 3.3. Collapse months
features_long <- aperm(features_array, c(1,3,2))
features_long <- array(features_long, dim = c(64800*12, 39))
dimnames(features_long)[[2]] <- dimnames(features_array)[[2]]

rm(features, CALL, features_array)
gc()

# --- 4. Re-compute the axis of variance
tmp <- aperm(data, c(1,4,2,3,5)) %>% array(., dim = c(64800*12, 21, 10, 4))
val_dt <- tmp %>% apply(., c(1,4), mean, na.rm = T) %>% apply(., 1, sd, na.rm = T) # sd across data types
val_h <- tmp %>% apply(., c(1,2), mean, na.rm = T) %>% apply(., 1, function(x)(x = mean(x[1:18] - x[4:21], na.rm = TRUE))) # hill difference

sup_long <- data.frame(dt = val_dt, h = val_h) # to be used as supplementary in the RDA

# --- 5. Compute multivariate analysis
# --- 5.1. Prepare data
ok <- complete.cases(data_long, features_long, sup_long)
X <- data_long[ok, ]

E <- features_long[ok, ]
nzv <- apply(E, 2, function(x) var(x, na.rm = TRUE) > 0)
E <- E[, nzv]

Z <- sup_long[ok,]

# --- 5.2. Perform PCA
# We use scaling = 1 because we care about correlation, so the angle between species x env
RDA <- vegan::rda(X = X, Y = as.data.frame(E), scale = T)
bio_scores <- scores(RDA, display = "species", scaling = 2, choices = c(1,2,3,4)) # bio scores
env_scores <- scores(RDA, display = "bp", scaling = 2, choices = c(1,2,3,4)) # env scores
sup_scores <- vegan::envfit(RDA, Z, permutations = 100, choices = c(1,2,3,4))$vectors$arrows # supplementary variables

# reverse RDA2 to have the association with primary production positive to ease interpretation
# (avoids double negative associations)
bio_scores[,2] <- bio_scores[,2]*(-1)
env_scores[,2] <- env_scores[,2]*(-1)
sup_scores[,2] <- sup_scores[,2]*(-1)

# --- 5.3. Extract contribution(s)
PC_percent <- (eigenvals(RDA)/sum(eigenvals(RDA))) # variance per PC
env_magnitude <- sqrt(rowSums(env_scores[, 1:2]^2)) # representativity of env. on PC
env_alpha <- (env_magnitude/max(env_magnitude))**2 # translate to alpha value on plot
sup_scores # print for sup fit

# --- 5.4. Plot
# --- 5.4.1. Pretty names
env_names <- row.names(env_scores) %>% 
  gsub("climatology_", "", .) %>% 
  gsub("_regridded", "", .) %>% 
  gsub("_SODA", "", .) %>% 
  gsub("_1998_2020_CMEMS", "", .) %>% 
  gsub("_aviso", "", .) %>% 
  gsub("_2001_2020", "", .) %>% 
  gsub("_CMEMS", "", .) %>% 
  gsub("woa18_", "", .) %>% 
  gsub("all_|decav_|decav81B0_", "", .) %>% 
  gsub("_m", "", .) %>% 
  gsub("_0_50", "", .) %>% 
  gsub("_factor", "", .) %>% 
  gsub("_100_1000", "", .)
row.names(env_scores) <- env_names

# --- 5.4.2. Variable colors
pal <- c("#413078","#287E8C","#9FD744","chocolate") %>% rev() %>% rep(., each = 21*10) %>% scales::alpha(., 0.1)
shape <- 16

# --- 5.4.3. Plot window 
par(mfrow = c(1,1), mar = c(6,6,6,1))
plot(bio_scores[,1]/max(abs(bio_scores)), bio_scores[,2]/max(abs(bio_scores)), main = "RDA Analysis of diversity index and environmental conditions",
     cex = 1, col = pal, pch = shape, axes = FALSE, xlim = c(-1.3,1.3), ylim = c(-0.8, 0.8))
arrows(x0 = 0, x1 = env_scores[,1], 
       y0 = 0, y1 = env_scores[,2],
       code = 2, length = 0.1, lwd = 1, col = scales::alpha("black", env_alpha))
text(x = env_scores[,1],
     y = env_scores[,2], 
     labels = env_names, col = scales::alpha("black", env_alpha), cex = 0.7,
     pos = ifelse(env_scores[,1] > 0, 4, 2))
abline(h = 0, v = 0, lty = "dotted", lwd = 2, col = "gray20")

points(sup_scores[,1:2], col = c("#DD2C45FF", "#368DAE"), cex = 1.3, pch = 15)
text(sup_scores[,1:2], labels = c("dt", "h"), pos = c(3,1), col = c("#DD2C45FF", "#368DAE"))

# --- 5.5. Environmental loadings
# --- 5.5.1. RDA1
par(mfrow = c(2,1))
pal <- rep("gray20", length(env_names))
pal[34] <- "red4"

barplot(sort(env_scores[,1]),
        ylim = c(-1,1), las = 2, col = pal, border = NA)
abline(v = seq(0.7,45.1, length.out = 38), col = pal, lty = "dotted")

# --- 5.5.2. RDA2
pal <- rep("gray20", length(env_names))
pal[38] <- "red4"

barplot(sort(env_scores[,2]),
        ylim = c(-0.4,0.6), las = 2, col = pal, border = NA)
abline(v = seq(0.7,45.1, length.out = 38), col = pal, lty = "dotted")

# 5.5.3. RDA3
pal <- rep("gray20", length(env_names))
pal[1] <- "red4"

barplot(sort(env_scores[,3]),
        ylim = c(-0.6,0.3), las = 2, col = pal, border = NA)
abline(v = seq(0.7,45.1, length.out = 38), col = pal, lty = "dotted")

# --- 5.6. Supplementary informations
print(PC_percent[1:4]) # RDA axis
print(cor(E, method = "pearson")) # env correlation
print(cor(X)) # pathway correlation

# --- 6. Response plot to environment
# --- 6.1. Define the function

response_2D <- function(x, y, z, col, range, xlab, ylab, main){
  # --- 6.1.1. Prepare the data
  df <- data.frame(x = x, y = y, z = z) # base data frame
  df$x_bin <- cut(df$x, breaks = 100) # bin environment
  df$y_bin <- cut(df$y, breaks = 100) # bin environment
  df <- df %>% 
    dplyr::group_by(x_bin, y_bin) %>% 
    dplyr::summarise(x = mean(x),
                     y = mean(y),
                     z = mean(z)) %>% 
    ungroup() # bin diversity
  
  r <- terra::rast(nrows = 100, ncols = 100,
    xmin = min(df$x), xmax = max(df$x), ymin = min(df$y), ymax = max(df$y)) # create raster
  cells <- terra::cellFromXY(r, as.matrix(df[, c("x", "y")])) # get cells id
  r[cells] <- df$z # assign values to cells id
  
  # --- 6.1.2. Do the plot
  plot(r, asp = NA, col = col, xlab = xlab, ylab = ylab, main = main, 
       range = range, fill_range = T)
  contour(r, add = T)
  grid(col = "black")
  
  return(NULL)
} # end function

# --- 6.2. Apply the plot
par(mfrow = c(3,2)) # zoom out to get them all

# --- 6.2.1. Global diversity
response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_t_0_50"],
            z = apply(X, 1, mean, na.rm = T), col = parula_pal(100), range = c(0,90),
            xlab = "Primary production (log)", ylab = "Temperature", main = "Diversity index response")

response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_fsle_aviso_2001_2020"],
            z = apply(X, 1, mean, na.rm = T), col = parula_pal(100), range = c(0,90),
            xlab = "Primary production (log)", ylab = "FSLE", main = "Diversity index response")

# --- 6.2.1. Data type SD
response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_t_0_50"],
            z = Z$dt, col = rocket_pal(100) %>% rev(), range = c(0,50),
            xlab = "Primary production (log)", ylab = "Temperature", main = "Data type SD response")

response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_fsle_aviso_2001_2020"],
            z = Z$dt, col = rocket_pal(100) %>% rev(), range = c(0,50),
            xlab = "Primary production (log)", ylab = "FSLE", main = "Data type SD response")

# --- 6.2.1. Hill scaling factor change
response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_t_0_50"],
            z = Z$h, col = curl_pal(100) %>% rev() %>% .[1:78], range = c(-5, 5),
            xlab = "Primary production (log)", ylab = "Temperature", main = "Hill change response")

response_2D(x = log(E[, "climatology_S_PP_regridded"]), y = E[, "climatology_fsle_aviso_2001_2020"],
            z = Z$h, col = curl_pal(100) %>% rev(), range = c(-5, 5),
            xlab = "Primary production (log)", ylab = "FSLE", main = "Hill change response")







# END




