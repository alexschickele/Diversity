#' =============================================================================
#' @name diversity_assessment
#' @description this functions sources the output folders of a CEPHALOPOD run in
#' order to perform a diversity assessment. The initializing by setting a run
#' name and sourcing the config file is necessary as in the start of CEPHALOPOD.
#' 
#' @param INPUT_FOLDER name of the CEPHALOPOD run folder we want the data from
#' @param OUTPUT_FOLDER name of the DIVERSITY folder to write the files in
#' @param HILL_NB Hill numbers to be computed for the diversity
#' @param BETA_DIV Boolean, if TRUE, computes the beta diversity as well
#' (the latter can be quite computing intensive)
#' @param FOCAL_SIZE radius (in cells) of the focal used to compute gamma
#' diversity (i.e., regional diversity)
#' @param MAX_CLUSTERS the maximum number of clusters for parallel computing
#' (overwrites the one in .config.R for Euler use essentially)
#' 
#' @return an RData containing all diversities as output

diversity_assessment <- function(INPUT_FOLDER,
                                 OUTPUT_FOLDER,
                                 HILL_NB = 0:3,
                                 BETA_DIV = FALSE,
                                 FOCAL_SIZE = 2,
                                 MAX_CLUSTERS = 10){
  
  # --- 1. Initialize & setup
  # --- 1.1. Load CALL
  load(paste0(INPUT_FOLDER,"/CALL.RData"))
  
  # --- 1.2. Initialize DIVERSITY save object
  DIVERSITY <- list()
  
  
  # --- 1.3. Dummy raster
  r0 <- raster::raster(nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90, 
                       vals=NA)
  
  
  # --- 1.4. Binarize function for richness
  binarize_hsi <- function(y_hat, q = seq(max(0, min(y_hat, na.rm = TRUE)), # security if min > 0.25
                                          min(1, max(y_hat, na.rm = TRUE)), # security if max < 0.75
                                          length.out = 50)){
    
    # Remove the extremes of the sequence (to avoid only 0's or 1's)
    q <- q[2:(length(q)-1)]
    
    # security if there is not enough values to have a sequence
    if(length(unique(q)) == 1){
      threshold <- unique(q)
      y_hat_bin <- y_hat
    } else {
      spearman <- 0
      threshold <- 0
      
      for(j in q){
        y_hat_bin <- y_hat
        y_hat_bin[y_hat_bin <= j] <- 0
        y_hat_bin[y_hat_bin > j] <- 1
        tmp <- cor(y_hat, y_hat_bin, method = "spearman", use = "pairwise.complete.obs")
        print(tmp)
        
        if(tmp > spearman & !is.na(tmp)){
          spearman <- tmp
          threshold <- j
        } # update
      } # for threshold
    } # end if
    
    y_hat_bin <- y_hat # reset
    y_hat_bin[y_hat_bin <= threshold] <- 0
    y_hat_bin[y_hat_bin > threshold] <- 1
    return(y_hat_bin)
    
  } # end function
  
  
  # --- 2. Build the ensemble(s)
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - START"))
  
  # --- 2.1. Extract file information
  # Which subfolder list and which model in the ensemble
  all_files <- list.files(INPUT_FOLDER, recursive = TRUE)
  model_files <- unique(dirname(all_files[grepl("MODEL.RData", all_files)])) #%>% .[1:30]
  ensemble_files <- mclapply(model_files, function(x){
    memory_cleanup() # low memory use
    
    load(paste0(INPUT_FOLDER,"/", x, "/MODEL.RData"))
    if(length(MODEL$MODEL_LIST) >= 1){
      load(paste0(INPUT_FOLDER,"/", x, "/QUERY.RData"))
      return(list(SUBFOLDER_NAME = x, MODEL_LIST = MODEL$MODEL_LIST, Y = QUERY$Y, MESS = QUERY$MESS))
    } else {
      return(NULL)
    } # if model list
  }, mc.cores = MAX_CLUSTERS) %>% .[lengths(.) != 0]
  
  # --- 2.2. Loop over the files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - loop over files"))
  tmp <- mclapply(ensemble_files, function(x){
    memory_cleanup() # low memory use
    
    # --- 2.2.1. Load MODEL files
    load(paste0(INPUT_FOLDER,"/", x$SUBFOLDER_NAME, "/MODEL.RData"))
    
    # --- 2.2.2. Extract projections in a matrix
    # If there is more than 1 algorithm, we extract and average across algorithm
    # Output matrix is cell x bootstrap x month
    if(length(x$MODEL_LIST) > 1){
      m <- lapply(x$MODEL_LIST, function(y){
        MODEL[[y]][["proj"]]$y_hat
      }) %>% abind(along = 4) %>% apply(c(1,2,3), function(z)(z = mean(z, na.rm = TRUE)))
    } else {
      m <- MODEL[[x$MODEL_LIST]][["proj"]]$y_hat
    } # end if
  }, mc.cores = MAX_CLUSTERS, mc.cleanup = TRUE)
  
  # --- 2.3. Stack in a cell x species x bootstrap x month matrix
  # --- 2.3.1. Re-arrange the array
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - format to array"))
  data <- tmp %>% 
    abind(along = 4) %>% 
    aperm(c(1,4,2,3))
  # --- 2.3.2. Pretty dimensions
  dimnames(data)[[2]] <- lapply(ensemble_files, function(x){out <- x$SUBFOLDER_NAME}) %>% unlist() %>% as.character()
  dimnames(data)[[4]] <- 1:12 %>% as.character()
  
  # --- 2.4. Memory cleanup
  rm(tmp)
  gc() # clean garbage and temporary files
  message(paste0(Sys.time(), "--- DIVERSITY: build the ensembles - DONE"))
  
  # --- 4. Binarize
  # --- 4.1. Computing the binary dataset
  message(paste0(Sys.time(), "--- DIVERSITY: binarize data - START"))
  data_bin <- mclapply(1:dim(data)[[3]], function(b){
    memory_cleanup() # low memory use
    
    data_bootstrap <- data[,,b,]
    out <- lapply(1:dim(data)[[4]], function(m){
      memory_cleanup() # low memory use
      
      df0 <- data_bootstrap[,,m] %>% as.data.frame()
      df0_bin <- apply(df0, 2, function(z)(z = binarize_hsi(y_hat = z)))
      return(df0_bin)
    }) %>% abind(along = 3)
    return(out)
  }, mc.cores = round(MAX_CLUSTERS/4, 0), mc.cleanup = TRUE) %>% abind(along = 4) %>% aperm(c(1,2,4,3)) # end lapply
  
  # --- 4.2. Pretty dimensions
  dimnames(data_bin) <- dimnames(data)
  
  # --- 4.3. Memory cleanup
  gc() # clean garbage and temporary files
  message(paste0(Sys.time(), "--- DIVERSITY: binarize data - DONE"))
  
  # --- 5. Initialize diversity computing
  full_cell_id <- which(!is.na(data[,1,1,1]))
  cell_vector <- getValues(r0) %>% as.numeric()
  loop_over <- expand.grid(HILL_NB, dimnames(data)[[3]], dimnames(data)[[4]])
  names(loop_over) <- c("Hill_value","Bootstrap","Month")
  
  # --- 6. Alpha diversities
  # --- 6.1. Compute the diversities
  message(paste0(Sys.time(), "--- DIVERSITY: compute alpha - START"))
  alpha_div_list <- mclapply(1:nrow(loop_over), function(x){
    memory_cleanup() # low memory use
    
    # --- 6.1.1. Taking the binary dataset if Hill equals 0
    if(loop_over$Hill_value[x] == 0 & CALL$DATA_TYPE == "binary"){
      df0 <- data_bin[full_cell_id,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
    } else{
      df0 <- data[full_cell_id,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
    } # end if Hill == 0
    
    # --- 6.1.2. Compute the diversity
    div0 <- hill_taxa(comm = df0, q = loop_over$Hill_value[x]) # compute diversity
    
    # --- 6.1.3. Assign cells back
    div <- cell_vector # base
    div[full_cell_id] <- div0 # fill non empty cells
    return(div)
    
  }, mc.cores = round(MAX_CLUSTERS/2, 0)) # limited due to memory use
  message(paste0(Sys.time(), "--- DIVERSITY: compute alpha - FORMATING"))
  
  # --- 6.2. Back transformation to array
  alpha_div_data <- array(data = NA,
                          dim = c(length(cell_vector), length(HILL_NB), dim(data)[[3]], dim(data)[[4]]),
                          dimnames = list(NULL,
                                          as.character(HILL_NB),
                                          dimnames(data)[[3]], 
                                          dimnames(data)[[4]]))
  
  for(x in 1:nrow(loop_over)){
    alpha_div_data[, as.character(loop_over$Hill_value[x]), loop_over$Bootstrap[x], loop_over$Month[x]] <- alpha_div_list[[x]]
  } # end for
  
  alpha_div_data[is.infinite(alpha_div_data)] <- NA # security
  message(paste0(Sys.time(), "--- DIVERSITY: compute alpha - DONE"))
  
  # --- 6.3. Save into diversity object
  DIVERSITY[["alpha"]] <- alpha_div_data
  message(paste0(Sys.time(), "--- DIVERSITY: save alpha - DONE"))
  
  # --- 7. Beta diversity
  if(BETA_DIV == TRUE){
    # --- 7.1. Compute the cell indices to average - for beta diversity
    message(paste0(Sys.time(), "--- DIVERSITY: get cells for beta - START"))
    beta_cell_id <- mclapply(1:length(full_cell_id), function(x){
      cell_id <- get_cell_neighbors(ID = full_cell_id[x], BUFFER = 2) %>% unlist()
      return(cell_id)
    }, mc.cores = 30) %>% bind_rows()
    message(paste0(Sys.time(), "--- DIVERSITY: get cells for beta - DONE"))
    
    # --- 7.2. Compute the beta diversity
    # We average the composition in the focal and compute the species overlap with the alpha
    message(paste0(Sys.time(), "--- DIVERSITY: compute beta - START"))
    beta_div_list <- mclapply(1:nrow(loop_over), function(x){
      memory_cleanup() # low memory use
      s = Sys.time()
      # --- 7.2.1. Build the focal community matrix
      df0_focal <- lapply(1:ncol(beta_cell_id[,-1]), function(z){
        memory_cleanup() # low memory use
        tmp <- beta_cell_id[,z] %>% unlist() %>% as.vector()
        
        # Taking the binary dataset if Hill equals 0
        if(loop_over$Hill_value[x] == 0 & CALL$DATA_TYPE == "binary"){
          df0_focal <- data_bin[tmp,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
        } else{
          df0_focal <- data[tmp,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
        } # end if Hill == 0
        
      }) %>% abind(along = 3) %>% 
        apply(c(1,2), function(zz)(zz = mean(zz, na.rm = TRUE))) # loop over the neighbor cells
      
      # --- 7.2.2. Build the center community matrix
      # Taking the binary dataset if Hill equals 0
      tmp <- beta_cell_id[,1] %>% unlist() %>% as.vector()
      if(loop_over$Hill_value[x] == 0 & CALL$DATA_TYPE == "binary"){
        df0_center <- data_bin[tmp,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
      } else{
        df0_center <- data[tmp,,loop_over$Bootstrap[x], loop_over$Month[x]] %>% as.data.frame() # subset full cells
      } # end if Hill == 0
      
      # --- 7.2.3. Compute diversity
      div0 <- lapply(1:nrow(df0_center), function(z){
        memory_cleanup() # low memory use
        df0 <- rbind(df0_center[z,], df0_focal[z,])
        # Security if communities are equal or composed of absences
        if(nrow(distinct(df0)) == 1 | sum(df0[1,]) == 0 | sum(df0[2,]) == 0){out <- 0} else{
          out <- hill_taxa_parti(comm = df0, q = loop_over$Hill_value[x], check_data = FALSE, show_warning = FALSE) %>% .[,"region_similarity"]
        } # end if security
        
      }) %>% unlist() # loop over cells
      
      # --- 7.2.4. Assign cells back
      div <- cell_vector # base
      div[full_cell_id] <- div0 # fill non empty cells
      Sys.time()-s
      return(div)
      
    }, mc.cores = round(MAX_CLUSTERS/2, 0)) # take car of the memory use !
    message(paste0(Sys.time(), "--- DIVERSITY: compute beta - FORMATTING"))
    
    # --- 7.3. Back transformation to array
    beta_div_data <- array(data = NA,
                           dim = c(length(cell_vector), length(HILL_NB), dim(data)[[3]], dim(data)[[4]]),
                           dimnames = list(NULL,
                                           as.character(HILL_NB),
                                           dimnames(data)[[3]], 
                                           dimnames(data)[[4]]))
    
    for(x in 1:nrow(loop_over)){
      beta_div_data[, as.character(loop_over$Hill_value[x]), loop_over$Bootstrap[x], loop_over$Month[x]] <- beta_div_list[[x]]
    } # end for
    
    beta_div_data[is.infinite(beta_div_data)] <- NA # security
    message(paste0(Sys.time(), "--- DIVERSITY: compute beta - DONE"))
    
    # --- 7.4. Save into diversity object
    DIVERSITY[["beta"]] <- beta_div_data
    message(paste0(Sys.time(), "--- DIVERSITY: save beta - DONE"))
    
  } # end if beta
  
  # --- 8. Compute MESS
  message(paste0(Sys.time(), "--- DIVERSITY: extract mess - START"))
  mess_data <- lapply(ensemble_files, function(x){
    mess <- getValues(x$MESS)
  }) %>% abind(along = 3) %>% aperm(c(1,3,2))
  DIVERSITY[["mess"]] <- mess_data
  message(paste0(Sys.time(), "--- DIVERSITY: extract mess - DONE"))
  
  # --- 9. Save diversity file
  message(paste0(Sys.time(), "--- DIVERSITY: save all - START"))
  save(DIVERSITY, file = paste0(OUTPUT_FOLDER, "/DIVERSITY.RData"),
       compress = "gzip", compression_level = 6)
  message(paste0(Sys.time(), "--- DIVERSITY: save all - DONE"))
  
} # END FUNCTION



