
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
# data_long <- aperm(data_diff, c(1,4,2,3,5))
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

# --- 4. Compute multivariate analysis
# --- 4.1. Prepare data
ok <- complete.cases(data_long, features_long)
X <- data_long[ok, ]

E <- features_long[ok, ]
nzv <- apply(E, 2, function(x) var(x, na.rm = TRUE) > 0)
E <- E[, nzv]

# --- 4.2. Perform PCA
# We use scaling = 1 because we care about correlation, so the angle between species x env
RDA <- vegan::rda(X = X, Y = as.data.frame(E), scale = T)
bio_scores <- scores(RDA, display = "species", scaling = 2, choices = c(1,2,3,4)) # bio scores
env_scores <- scores(RDA, display = "bp", scaling = 2, choices = c(1,2,3,4)) # env scores

# reverse RDA2 to have the association with primary production positive to ease interpretation
# (avoids double negative associations)
bio_scores[,2] <- bio_scores[,2]*(-1)
env_scores[,2] <- env_scores[,2]*(-1)

# --- 4.3. Extract contribution(s)
PC_percent <- (eigenvals(RDA)/sum(eigenvals(RDA))) # variance per PC
env_magnitude <- sqrt(rowSums(env_scores[, 1:2]^2)) # representativity of env. on PC
env_alpha <- (env_magnitude/max(env_magnitude))**2 # translate to alpha value on plot

# --- 4.4. Plot
# --- 4.4.1. Pretty names
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
  gsub("_100_1000", "", .)
row.names(env_scores) <- env_names

# --- 4.4.2. Variable colors
pal <- rocket_pal(21) %>% rep(., 10*4) %>% scales::alpha(., 0.3) # hill palette
shape <- 15:18 %>% rep(., each = 21*10) # type point

# --- 4.4.3. Plot window 
par(mfrow = c(1,1), mar = c(6,6,6,1))
plot(bio_scores[,1]/max(abs(bio_scores)), bio_scores[,2]/max(abs(bio_scores)), main = "RDA Analysis of diversity index and environmental conditions",
     cex = 2, col = pal, pch = shape, axes = FALSE, xlim = c(-1.3,1.3), ylim = c(-0.45, 0.6))
arrows(x0 = 0, x1 = env_scores[,1], 
       y0 = 0, y1 = env_scores[,2],
       code = 2, length = 0.1, lwd = 1, col = scales::alpha("black", env_alpha))
text(x = env_scores[,1],
     y = env_scores[,2], 
     labels = env_names, col = scales::alpha("black", env_alpha), cex = 1.2,
     pos = ifelse(env_scores[,1] > 0, 4, 2))
abline(h = 0, v = 0, lty = "dotted", lwd = 2, col = "gray20")

# --- 4.5. Environmental loadings
# --- 4.5.1. RDA1
par(mfrow = c(2,1))
pal <- rep("gray20", length(env_names))
pal[34] <- "red4"

barplot(sort(env_scores[,1]),
        ylim = c(-1,1), las = 2, col = pal, border = NA)
abline(v = seq(0.7,45.1, length.out = 38), col = pal, lty = "dotted")

# --- 4.5.2. RDA2
pal <- rep("gray20", length(env_names))
pal[1] <- "red4"

barplot(sort(env_scores[,2]) %>% rev(),
        ylim = c(-0.4,0.6), las = 2, col = pal, border = NA)
abline(v = seq(0.7,45.1, length.out = 38), col = pal, lty = "dotted")

# --- 4.6. Supplementary informations
print(PC_percent[1:4]) # RDA axis
print(cor(E, method = "pearson")) # env correlation
print(cor(X)) # pathway correlation

# --- 5. Biological loadings
# --- 5.1. Prepare dataset
bio_loadings <- bio_scores %>%  cbind(data_names_long) 
bio_loadings$Type <- factor(bio_loadings$Type, levels = c("m","a","b","o"))
bio_loadings$Hill_num <- as.numeric(gsub("Hill ", "", bio_loadings$Hill))

pal <- rocket_pal(21) %>% rep(., 10*4) %>% scales::alpha(., 0.3) # hill palette
pal <- scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6) %>% rep(., each = 21*10) # type point

# --- 5.2. Global diversity
par(mar = c(10,5,10,5))
boxplot(bio_loadings$RDA2, bio_loadings$RDA1, outline = FALSE,
        horizontal = TRUE, col = scales::alpha("gray20", 0.5) , border = "gray20")
stripchart(list(bio_loadings$RDA2, bio_loadings$RDA1),
           method = "jitter",
           jitter = 0.2,
           pch = 16,
           col = scales::alpha("black", 0.2),
           add = TRUE)
abline(v = 0)

# --- 5.3. By data type
# --- 5.3.1. RDA1
par(mar = c(5,5,5,5))
boxplot(bio_loadings$RDA1 ~ bio_loadings$Type, outline = FALSE, 
        horizontal = TRUE, col = scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6) %>% rev(.), border = "gray20")

stripchart(bio_loadings$RDA1 ~ bio_loadings$Type,
           method = "jitter",
           jitter = 0.2,
           pch = 16,
           col = scales::alpha("black", 0.2),
           add = TRUE)

# --- 5.3.2. RDA2
par(mar = c(5,5,5,5))
boxplot(bio_loadings$RDA2 ~ bio_loadings$Type, outline = FALSE,
        horizontal = TRUE, col = scales::alpha(c("#413078", "#287E8C", "#9FD744","chocolate1"), 0.6) %>% rev(.), border = "gray20")

stripchart(bio_loadings$RDA2 ~ bio_loadings$Type,
           method = "jitter",
           jitter = 0.2,
           pch = 16,
           col = scales::alpha("black", 0.2),
           add = TRUE)
abline(v = 0)

# --- 5.4. By Hill
# --- 5.3.1. RDA1
pal <- c(o = "#413078", a = "#287E8C", b = "#9FD744", m = "chocolate1")
par(mfrow = c(1,1), mar = c(3,5,3,5))
boxplot(bio_loadings$RDA1 ~ bio_loadings$Hill, outline = FALSE, ylim = c(1,4.5),
        horizontal = TRUE, col = scales::alpha("white", 0.6), border = "gray20")

lapply(c("o","a","b","m"), function(x){
  tmp <- bio_loadings %>% 
    dplyr::filter(Type == x)
  
  stripchart(tmp$RDA1 ~ tmp$Hill,
             method = "jitter",
             jitter = 0.2,
             pch = 16,
             col = pal[x] %>% scales::alpha(., 0.3),
             add = TRUE)
  abline(v = 0)
  
  fit <- lm(RDA1 ~ Hill_num, data = tmp)
  pred <- stats::predict(fit, 
                         newdata = data.frame(Hill_num = seq(0, 5, length.out = 100)))
  lines(pred, seq(1, length(hill_ref), length.out = 100), col = pal[x], lwd = 3)
})

# --- 5.3.1. RDA2
pal <- c(o = "#413078", a = "#287E8C", b = "#9FD744", m = "chocolate1")
par(mfrow = c(1,1), mar = c(3,5,3,5))
boxplot(bio_loadings$RDA2 ~ bio_loadings$Hill, outline = FALSE,
        horizontal = TRUE, col = scales::alpha("white", 0.6), border = "gray20")

lapply(c("o","a","b","m"), function(x){
  tmp <- bio_loadings %>% 
    dplyr::filter(Type == x)
  
  stripchart(tmp$RDA2 ~ tmp$Hill,
             method = "jitter",
             jitter = 0.2,
             pch = 16,
             col = pal[x] %>% scales::alpha(., 0.3),
             add = TRUE)
  abline(v = 0)
  
  fit <- lm(RDA2 ~ Hill_num, data = tmp)
  pred <- stats::predict(fit, 
                         newdata = data.frame(Hill_num = seq(0, 5, length.out = 100)))
  lines(pred, seq(1, length(hill_ref), length.out = 100), col = pal[x], lwd = 3)
})

