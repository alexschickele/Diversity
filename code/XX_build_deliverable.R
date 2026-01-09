
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

rm(data_abundance)
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

# --- 3. Build the 5D array
# Hill 3 didnt work for abundance, I will just skip it for all
m <- m[,-13,,]
o <- o[,-13,,]
b <- b[,-13,,]

all_div <- abind(m,a,b,o, along = 5) # all together
dim(all_div) <- c(360, 180, dim(all_div)[2:5]) # reshape

all_div[is.na(all_div)] <- -9999.9
all_div <- all_div %>% aperm(c(1,2,5,3,4,6))

# --- 4. Create the netcdf
# --- 4.1. Set up
library(RNetCDF)
nc <- create.nc("./output/B5D_M15_global_diversity", format = "netcdf4")

# --- 4.2. Global attributes
att.put.nc(nc, variable = "NC_GLOBAL", name = "Conventions", type = "NC_CHAR", value = "CF-1.12")
att.put.nc(nc, variable = "NC_GLOBAL", name = "title", type = "NC_CHAR", value = "Merged diversity dataset - Biocean5D")
att.put.nc(nc, variable = "NC_GLOBAL", name = "institution", type = "NC_CHAR", value = "ETH Zurich")
att.put.nc(nc, variable = "NC_GLOBAL", name = "source", type = "NC_CHAR", value = "ETH Zurich")
att.put.nc(nc, variable = "NC_GLOBAL", name = "history", type = "NC_CHAR", value = paste(Sys.time(), "File created"))
att.put.nc(nc, variable = "NC_GLOBAL", name = "comment", type = "NC_CHAR", value = "Uses attributes recommended by http://cfconventions.org")
att.put.nc(nc, variable = "NC_GLOBAL", name = "references", type = "NC_CHAR", value = "https://bio-oracle.org")

# --- 4.3. Dimensions
# Longitude dimension
dim.def.nc(nc, dimname = "lon", dimlength = 360)
var.def.nc(nc, varname = "lon", vartype = "NC_FLOAT", dimensions = "lon")
att.put.nc(nc, variable = "lon", name = "standard_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "long_name", type = "NC_CHAR", value = "longitude")
att.put.nc(nc, variable = "lon", name = "units", type = "NC_CHAR", value = "degrees_east")

# Latitude dimension
dim.def.nc(nc, dimname = "lat", dimlength = 180)
var.def.nc(nc, varname = "lat", vartype = "NC_FLOAT", dimensions = "lat")
att.put.nc(nc, variable = "lat", name = "standard_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "long_name", type = "NC_CHAR", value = "latitude")
att.put.nc(nc, variable = "lat", name = "units", type = "NC_CHAR", value = "degrees_north")

# Time dimension
dim.def.nc(nc, dimname = "time", dimlength = 12)
var.def.nc(nc, varname = "time", vartype = "NC_FLOAT", dimensions = "time")
att.put.nc(nc, variable = "time", name = "standard_name", type = "NC_CHAR", value = "month")
att.put.nc(nc, variable = "time", name = "long_name", type = "NC_CHAR", value = "month")
att.put.nc(nc, variable = "time", name = "units", type = "NC_CHAR", value = "month from 1 to 12")

# Bootstrap dimension
dim.def.nc(nc, dimname = "bootstrap_uncertainty", dimlength = 10)
var.def.nc(nc, varname = "bootstrap_uncertainty", vartype = "NC_FLOAT", dimensions = "bootstrap_uncertainty")
att.put.nc(nc, variable = "bootstrap_uncertainty", name = "standard_name", type = "NC_CHAR", value = "bootstrap_uncertainty")
att.put.nc(nc, variable = "bootstrap_uncertainty", name = "long_name", type = "NC_CHAR", value = "bootstrap replicate for uncertainty")
att.put.nc(nc, variable = "bootstrap_uncertainty", name = "units", type = "NC_CHAR", value = "bootstrap from 1 to 10")

# Method dimension
dim.def.nc(nc, dimname = "hill_diversity", dimlength = 20)
var.def.nc(nc, varname = "hill_diversity", vartype = "NC_FLOAT", dimensions = "hill_diversity")
att.put.nc(nc, variable = "hill_diversity", name = "standard_name", type = "NC_CHAR", value = "hill_diversity")
att.put.nc(nc, variable = "hill_diversity", name = "long_name", type = "NC_CHAR", value = "hill scaling factor diversity")
att.put.nc(nc, variable = "hill_diversity", name = "units", type = "NC_CHAR", value = "hill scaling factor from 0 to 5, steps of 0.25")

# Type dimension
dim.def.nc(nc, dimname = "observation_type", dimlength = 4)
var.def.nc(nc, varname = "observation_type", vartype = "NC_FLOAT", dimensions = "observation_type")
att.put.nc(nc, variable = "observation_type", name = "standard_name", type = "NC_CHAR", value = "observation_type")
att.put.nc(nc, variable = "observation_type", name = "long_name", type = "NC_CHAR", value = "observation_type")
att.put.nc(nc, variable = "observation_type", name = "units", type = "NC_CHAR", value = "observation_type, including metagenomic, abundance, biomass and occurrence")

# --- 4.4. Variables
# CRS variable
var.def.nc(nc, varname = "crs", vartype = "NC_CHAR", dimensions = NA)
att.put.nc(nc, variable = "crs", name = "grid_mapping_name", type = "NC_CHAR", value = "latitude_longitude")
att.put.nc(nc, variable = "crs", name = "long_name", type = "NC_CHAR", value = "CRS definition")
att.put.nc(nc, variable = "crs", name = "longitude_of_prime_meridian", type = "NC_DOUBLE", value = 0.)
att.put.nc(nc, variable = "crs", name = "semi_major_axis", type = "NC_DOUBLE", value = 6378137.)
att.put.nc(nc, variable = "crs", name = "inverse_flattening", type = "NC_DOUBLE", value = 298.257223563)
att.put.nc(nc, variable = "crs", name = "spatial_ref", type = "NC_CHAR", value = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AXIS[\"Latitude\",NORTH],AXIS[\"Longitude\",EAST],AUTHORITY[\"EPSG\",\"4326\"]]')
att.put.nc(nc, variable = "crs", name = "GeoTransform", type = "NC_CHAR", value = '-180 0.08333333333333333 0 90 0 -0.08333333333333333 ')

# --- 4.5. Diversity
# Any random variable
var.def.nc(nc, varname = "diversity", vartype = "NC_FLOAT", dimensions = c("lon", "lat", "time", "hill_diversity","bootstrap_uncertainty","observation_type"))
att.put.nc(nc, variable = "diversity", name = "long_name", type = "NC_CHAR", value = "Modeled diversity")
att.put.nc(nc, variable = "diversity", name = "units", type = "NC_CHAR", value = "Quantile scale from 1 to 100")

# Other recommended attributes
att.put.nc(nc, variable = "diversity", name = "grid_mapping", type = "NC_CHAR", value = "crs")
att.put.nc(nc, variable = "diversity", name = "_FillValue", type = "NC_FLOAT", value = -9999.9)

# Close and check
sync.nc(nc)
print.nc(nc)

# Write values
var.put.nc(nc, "diversity", all_div)

# Close and check
sync.nc(nc)
print.nc(nc)

close.nc(nc)





