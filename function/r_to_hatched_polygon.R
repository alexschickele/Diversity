#' @name r_to_hatched_polygon
#' @description This function converts a raster object into hatched polygons based on a specified threshold.
#' @param RASTER A raster object representing quantitative values.
#' @param THRESHOLD A numeric value indicating the threshold above which to create hatched polygons.
#' @param COL The color to fill the polygons. Default is "black".
#' @param DENSITY The density of the hatch lines. Default is 20.
#' @param BORDER The color of the polygon borders. Default is NA (no border).
#'
#' @return Hatched polygons representing areas where raster values are above the threshold.

r_to_hatched_polygon <- function(RASTER,
                                 THRESHOLD,
                                 COL = "black",
                                 DENSITY = 20,
                                 BORDER = NA){
  
  # Add a box around the raster at the minimum value - to avoid polygon outside the raster limit
  tmp <- as.matrix(RASTER)
  tmp[1,] <- min(tmp, na.rm = TRUE)
  tmp[180,] <- min(tmp, na.rm = TRUE)
  tmp[,1] <- min(tmp, na.rm = TRUE)
  tmp[,360] <- min(tmp, na.rm = TRUE)
  RASTER <- setValues(RASTER, tmp)
  
  # Create a contour line at the threshold
  rc0 <- raster::rasterToContour(RASTER > THRESHOLD, nlevels = 1)
  
  # Close the lines (to form a polygon) and plot it
  all_lines <- rc0@lines[[1]]@Lines
  all_poly <- lapply(1:length(all_lines), function(id){
    poly_bounds <- rbind(all_lines[[id]]@coords, all_lines[[id]]@coords[1,])
    polygon(x = poly_bounds[,1], y = poly_bounds[,2], col = COL, density = DENSITY, border = BORDER)
    
  }) # end polygon plot loop
  
} # end function
