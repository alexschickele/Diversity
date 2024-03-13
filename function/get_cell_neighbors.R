#' @name get_cell_neighbors
#' @description little utility function to identify neighboring cells from a user
#' defined one. Useful for beta diversity and local averages... for example.
#' @param NX number of X columns in the spatial grid
#' @param NY number of Y columns in the spatial grid
#' @param ID cell identified (i.e., as a unique cell index in a matrix; 1:NX*NY)
#' @param BUFFER radius in which we extract the neighbors (i.e., in number of cells)
#' @return a list containing the ID the center and neighbors cells

get_cell_neighbors <- function(NX = 360,
                               NY = 180,
                               ID = 1,
                               BUFFER = 1){
  
  # --- 1. Initialize
  # --- 1.1. Libraries
  require(raster)
  
  # --- 1.2. Create base grid - in raster format
  # We do not care about the extent or coordinates. 
  r <- raster(ncols = NX, nrows = NY)
  
  # --- 2. Get center row and col
  center_id <- ID
  center_rc <- rowColFromCell(r, ID)
  
  # --- 3. Get neighbors ID
  # --- 3.1. Extract all cells from the area
  neighbors_id <- cellFromRowColCombine(r, 
                                     col = (center_rc[2]-BUFFER):(center_rc[2]+BUFFER),
                                     row = (center_rc[1]-BUFFER):(center_rc[1]+BUFFER))
  # --- 3.2. Remove center cell and NA
  # NA cell are outside the RASTER extent
  neighbors_id <- neighbors_id[which(neighbors_id != ID & !is.na(neighbors_id))] 
  
  # --- 4. Wrap up and save
  out <- list(center_id = center_id,
              neighbors_id = neighbors_id)
  
  return(out)
} # END FUNCTION
