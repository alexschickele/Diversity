#' @description
#' Function to test whether all components of a map are significantly different on the
#' factor representing the Z dimension.
#' 
#' @example Which geographical cell present significant distribution in their diversity 
#' resulting from different data types ?
#' 
#' @param A longitude x latitude x factor array
#' @param R resolution of the window
#' @param THRESHOLD for the p-value
#' @return A lon x lat two column data frame of the centroids of window that are not significant

moving_window_anova <- function(A, R, THRESHOLD){

  
  # --- A. Define the window resolution
  xy_id <- expand.grid(lon = seq(1, dim(A)[[1]], R),
                       lat = seq(1, dim(A)[[2]], R))
  
  # --- B. Loop over to apply ANOVA
  pval <- lapply(1:nrow(xy_id), function(x){
    #' @param W Extracted moving window
    
    # --- Extract, flatten and remove NA
    W <- A[xy_id[x,1]:(xy_id[x,1]+R-1), xy_id[x,2]:(xy_id[x,2]+R-1), ] %>% 
      apply(3, as.vector)
    id <- which(!is.na(apply(W, 1, sum)))
    W <- W[id, ] # Window with na.rm = T
    
    # --- If enough values perform ANOVA
    if(sum(W) == 0){out <- NA}
    else if(length(id) < 3){out <- 1} else {
      
      val = as.vector(W)
      group = rep(1:ncol(W), each = nrow(W)) %>% as.factor()
      out <- aov(val ~ group) %>% summary() %>% .[[1]] %>% .[1, "Pr(>F)"]
    } # end if
    return(out)
    
  }) %>% unlist() # end lapply
  # --- C. Concatenate in a dataframe and return
  pval_result <- data.frame(lon = xy_id$lon - 180 + (R/2),
                            lat = (xy_id$lat - 90 + (R/2)) %>% rev(),
                            pval = pval) %>% dplyr::filter(pval > THRESHOLD)
  return(pval_result)
} # end function
