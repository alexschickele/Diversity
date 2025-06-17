#' =============================================================================
#' @name diversity_plot
#' @description Function to generate diversity plots from the output of "diversity_assessment.R"
#' 
#' @param INPUT_FOLDER name of the CEPHALOPOD run folder we want the data from (i.e., for the CALL object)
#' @param OUTPUT_FOLDER The path to the folder where the diversity data is stored (i.e., for the DIVERSITY.RData object).
#' @param PROFILE Name of the environmental variable to compute the profile from (sensu CALL$ENV_VAR). Set to NULL for no profile plots.
#' @return Generates diversity plots and saves them in a PDF file.


diversity_plot <- function(INPUT_FOLDER,
                           OUTPUT_FOLDER,
                           PROFILE = NULL){
  
  # --- 1. Initialize
  # --- 1.1. Dummy raster
  r0 <- raster::raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, 
                       vals=NA)
  
  # --- 1.2. Load diversity object
  load(paste0(OUTPUT_FOLDER, "/DIVERSITY.RData"))
  HILL_NB <- DIVERSITY$alpha %>% dimnames() %>% .[[2]]
  
  # --- 1.3. Load CALL object - for environmental variables
  load(paste0(INPUT_FOLDER, "/CALL.RData"))
  
  # --- 1.4. Graphical data object
  GRAPH_DATA <- list()
  
  # --- 2. Compute the data tables
  # --- 2.1. For alpha diversity
  # --- 2.1.1. Global average
  div_avg <- apply(DIVERSITY$alpha, c(1,2), function(x)(x = mean(x, na.rm = TRUE)))
  div_avg_r <- lapply(1:dim(div_avg)[[2]], function(x){
    r <- setValues(r0, div_avg[,x])
  }) %>% stack()
  div_avg_r[is.na(div_avg_r)] <- 0
  GRAPH_DATA[["alpha"]][["avg"]] <- div_avg_r
  
  # --- 2.1.2. Standard deviation
  div_avg <- apply(DIVERSITY$alpha, c(1,2), function(x)(x = sd(x, na.rm = TRUE)/max(x, na.rm = TRUE)))
  div_avg_r <- lapply(1:dim(div_avg)[[2]], function(x){
    r <- setValues(r0, div_avg[,x])
  }) %>% stack()
  div_avg_r[is.na(div_avg_r)] <- 0
  GRAPH_DATA[["alpha"]][["sd"]] <- div_avg_r
  
  # --- 2.2. For beta diversity
  if(is.null(DIVERSITY$beta) == FALSE){
    # --- 2.2.1. Global average
    div_avg <- apply(DIVERSITY$beta, c(1,2), function(x)(x = mean(x, na.rm = TRUE)))
    div_avg_r <- lapply(1:dim(div_avg)[[2]], function(x){
      r <- setValues(r0, div_avg[,x])
    }) %>% stack()
    div_avg_r[is.na(div_avg_r)] <- 0
    GRAPH_DATA[["beta"]][["avg"]] <- div_avg_r
    
    # --- 2.2.2. Standard deviation
    div_avg <- apply(DIVERSITY$beta, c(1,2), function(x)(x = sd(x, na.rm = TRUE)/max(x, na.rm = TRUE)))
    div_avg_r <- lapply(1:dim(div_avg)[[2]], function(x){
      r <- setValues(r0, div_avg[,x])
    }) %>% stack()
    div_avg_r[is.na(div_avg_r)] <- 0
    GRAPH_DATA[["beta"]][["sd"]] <- div_avg_r
  } # end if beta
  
  # --- 2.3. MESS
  MESS_DATA <- apply(DIVERSITY$mess, 1, function(x){x[x > 0] <- 1
                                                    x[x < 0] <- -1
                                                    x <-  median(x, na.rm = TRUE)})
  MESS_R <- setValues(r0, -MESS_DATA) # we take the opposite because the threshold of the plot is "above"
  
  # --- 2.4. Land
  land <- GRAPH_DATA[[1]][[1]][[1]]
  land[land == 0] <- 9999
  land[land != 9999] <- NA
  
  # --- 2.5. Color palette
  pal <- colorRampPalette(c("white",brewer.pal(9, "GnBu"),"black"))(100)

  # --- 3. Mapping
  # --- 3.1. Initialize
  pdf(paste0(OUTPUT_FOLDER, "/diversity_maps.pdf"))
  par(mfrow = c(4,2), mar = c(2,3,2,3))
  
  # --- 3.2. Plot iteratively
  for(i in 1:length(GRAPH_DATA)){
    for(h in 1:length(HILL_NB)){
      # --- 3.2.1. Diversity value
      plot(GRAPH_DATA[[i]][["avg"]][[h]], col = pal,
           main = paste(names(GRAPH_DATA)[i], "div. ( Hill =", HILL_NB[h],")"), cex.axis = 0.6, cex.main = 0.6)
      plot(raster::rasterToContour(GRAPH_DATA[[i]][["avg"]][[h]], nlevels = 4), lwd = 0.5, add = TRUE, col = "gray20")
      box("figure", col="black", lwd = 1)
      
      # --- 3.2.2. Standard deviation
      threshold <- 0.25
      if(max(getValues(GRAPH_DATA[[i]][["sd"]][[h]]), na.rm = TRUE) > threshold){
        r_to_hatched_polygon(RASTER = GRAPH_DATA[[i]][["sd"]][[h]],
                             THRESHOLD = threshold,
                             COL = "red3",
                             DENSITY = 40,
                             BORDER = NA)
      } # if enough SD
      
      # --- 3.2.3. MESS data
      threshold <- 0
      if(max(getValues(MESS_R), na.rm = TRUE) > threshold){
        r_to_hatched_polygon(RASTER = MESS_R,
                             THRESHOLD = threshold,
                             COL = "gray80",
                             DENSITY = 20,
                             BORDER = NA)
      } # if enough SD
      
      # --- 3.2.4. Land
      plot(land, col = "antiquewhite4", legend = FALSE, add = TRUE)
      
    } # for h hill
  } # for i
  
  # --- 3.3. Wrap up
  dev.off()
  
  # --- 4. Profiles
  # --- 4.1. Extract the environmental data
  tmp<- lapply(CALL$ENV_DATA, function(x){
    out <- x[[PROFILE]] %>% getValues()
  }) %>% abind(along = 2)
  PROFILE_DATA <- array(tmp, dim = c(dim(tmp), dim(DIVERSITY[[1]])[[3]])) %>% aperm(c(1,3,2)) %>% as.integer()
  
  # --- 4.2. Plot the profiles
  pdf(paste0(OUTPUT_FOLDER, "/diversity_profiles_",PROFILE,".pdf"))
  par(mfrow = c(3,2), mar = c(6,4,3,2))
  
  for(i in 1:length(GRAPH_DATA)){
    for(h in 1:length(HILL_NB)){
      # --- 4.2.1 Build a frequency table diversity x env.
      freq_m <- table(PROFILE_DATA, 
                      DIVERSITY[[i]][,h,,] %>% as.integer())
    
    # --- 4.2.2 Plot iteratively
    image(freq_m, col = pal, axes = FALSE, xlab = PROFILE, ylab = paste("Hill = ", HILL_NB[h]), 
          main = "Hill diversity profile", sub = "Color scale represent the log(frequency)")
    axis(side = 1, at = seq(0,1,0.1), as.integer(seq(min(PROFILE_DATA, na.rm = TRUE), max(PROFILE_DATA, na.rm = TRUE), length.out = 11)))
    axis(side = 2, at = seq(0,1,0.1), as.integer(seq(min(DIVERSITY[[i]][,h,,], na.rm = TRUE), max(DIVERSITY[[i]][,h,,], na.rm = TRUE), length.out = 11)), las = 2)
    box()
    abline(h = seq(0,1,0.02), v = seq(0,1,0.02), col = scales::alpha("gray20", 0.2))
    abline(h = seq(0,1,0.1), v = seq(0,1,0.1), col = scales::alpha("gray20", 0.5))
    box("figure", col="black", lwd = 1)
    } # end for h hill
  } # end for i div
  
  # --- 4.3. Wrap up
  dev.off()
  
} # END FUNCTION
