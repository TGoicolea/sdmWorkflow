


.getBestComb <- function(species,preds,n=1,comb,bgn=1000) {
  # The following two lines are not needed if we sure the input is presence-only with "species" column
  # but to make sure, I just re-generated the SpatVector!!
  .sp <-vect(crds(species)) # make a SpatVector with coordinates of the presence locations
  .sp$species <- 1 # add species column
  
  .auc <- rep(0,length(comb))
  #------
  for (i in seq_along(comb)) {
    .w <- comb[[i]] # which variables/layers combination in item i of the comb list?
    .prs <- preds[[.w]] # selected variable combinations!
    d <- sdmData(species~., .sp, .prs, bg=list(method='gRandom',n=bgn))
    m <- sdm(species~., d, methods = c('glmp','brt','svm','bioclim.dismo','mda','maxent','rf','mlp','cart'), #,'mda'maxent','rf','mlp''cart',
             n=n,replication='boot')
    e <- .getEval(m,n='species',setting=list(method='weighted',stat='auc')) # from model -> get Evaluation based on Ensemble
    .auc[i] <- e@statistics$AUC
  }
  
  cat('\n The top 10 combinations (sorted): ',paste(order(.auc,decreasing = T),collapse = ', '))
  
  .auc
}
#-------


# The core function in Workflow (can be used in a more efficient parallelisation!)
.coreWF <- function(i,preds,species,n=10,comb,bgn=1000,path='.',filename,spn,...) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  # The following two lines are not needed if we sure the input is presence-only with "species" column
  # but to make sure, I just re-generated the SpatVector!!
  .sp <- vect(crds(species)) # make a SpatVector with coordinates of the presence locations
  .sp$species <- 1 # add species column
  
  .w <- comb[[i]] # which variables/layers combination in item i of the comb list?
  .prs <- preds[[.w]] # selected variable combinations!
  d <- sdmData(species~., .sp, .prs, bg=list(method = 'gRandom',n=bgn))
  m <- sdm(species~., d, methods = c('glmp','brt','rf','svm','cart','maxent','bioclim.dismo','mlp','mda'),
           n=n,replication='boot')
  #----
  if (path == '.') write.sdm(m, filename = paste0('models/M_',i,'_',spn,'_',filename,'.sdm'),overwrite=TRUE)
  else write.sdm(m, filename = paste0(path,'/models/M_',i,'_',spn,'_',filename,'.sdm'),overwrite=TRUE)
  #----------
  # Save predictions:
  if (!dir.exists(file.path(path, "predicts"))) dir.create(file.path(path, "predicts"), recursive = TRUE, showWarnings = FALSE)
  if (path == '.') .pr <- predict(m, .prs, filename = paste0('predicts/pr_COMB-',i,'_',spn,'_',filename,'.tif'),overwrite=T)
  else .pr <- predict(m, .prs, filename = paste0(path,'/predicts/pr_COMB-',i,'_',spn,'_',filename,'.tif'),overwrite=T)
  
  #------
  if (!dir.exists(file.path(path, "ensembles"))) dir.create(file.path(path, "ensembles"), recursive = TRUE, showWarnings = FALSE)
  if (path == '.') en <- ensemble(m, .pr, setting=list(method='weighted',stat='auc'),filename = paste0('ensembles/en_COMB-',i,'_',spn,'_',filename,'.tif'),overwrite=T)
  else en <- ensemble(m, .prs, setting=list(method='weighted',stat='auc'),filename = paste0(path,'/ensembles/en_COMB-',i,'_',spn,'_',filename,'.tif'),overwrite=T) #originally was .prs, but I changed to .prs
  #---------
  rm(m); gc()
  TRUE
}
#---------
# error-handling of WF:
### .wf_ER <- function(i,...) {
#   e <- try(.coreWF(i,...),silent=TRUE)
#   if (inherits(e,'try-error')) FALSE
#   else TRUE
# }
#second version
.wf_ER <- function(i, preds_path, species_path, n, comb, bgn, path, filename, spn) {
  species <- vect(species_path)  # Load species data reading the path instead of sending species data as function argument (which might be too large and was causing problems), only pass file paths and read them inside .wf_ER
  preds <- rast(preds_path)  # Load predictors reading the path instead of sending predictors as function arguments (which might be too large), only pass file paths and read them inside .wf_ER
  e <- try(.coreWF(i, preds, species, n, comb, bgn, path, filename, spn), silent = TRUE)
  if (inherits(e, 'try-error')) FALSE else TRUE
}
#--------------



runWF <- function(species,preds,n=10,comb,bgn=1000,path='.',filename,ncore=8,spn) {
  if (!is.null(ncore)) {
    ncore <- c(parallel::detectCores()-2,ncore)
    ncore <- ncore[which.min(ncore)]
  }
  
  # The following two lines are not needed if we sure the input is presence-only with "species" column
  # but to make sure, I just re-generated the SpatVector!!
  .sp <- vect(crds(species)) # make a SpatVector with coordinates of the presence locations
  .sp$species <- 1 # add species column
  
  if (!is.null(ncore) && ncore > 1) {
    # we need to save the files, and read them in parallel sessions (revise the code if you already have files)
    species_path<-paste0('temporal/_temp_species_file_',spn,'_',filename,'.shp')
    preds_path<-paste0('temporal/_temp_predictors_file_',spn,'_',filename,'.tif')
    writeVector(.sp,species_path,overwrite=T)
    writeRaster(preds,filename=preds_path,overwrite=T)
    
    library(parallel)
    cl <- makeCluster(ncore)
    #cl <- parallel::makeForkCluster(ncore)
    ###cl <- parallel::makePSOCKcluster(ncore)  # Explicit and avoids confusion
    
    clusterEvalQ(cl, {
      library(terra)
      library(sdm)
      getmethodNames()
      # .sp <- vect(paste0('temporal/_temp_species_file_',spn,'_',filename,'.shp'))
      # preds <- rast(paste0('temporal/_temp_predictors_file_',spn,'_',filename,'.tif'))
      NULL
    })
    ###clusterExport(cl, c("path","filename","spn","comb","n","bgn",".coreWF"))
    ###2    clusterExport(cl, varlist = c("path", "filename", "comb", "n", "bgn", ".coreWF", "spn", "preds", "species"), envir = environment())
    clusterExport(cl, varlist = c("path", "filename", "comb", "n", "bgn", ".coreWF", "spn",".wf_ER","species_path","preds_path"), envir = environment())
    
    
    ### o <- parLapply(cl,1:length(comb),.wf_ER)
    o <- parLapply(cl, 1:length(comb), function(i) .wf_ER(i, preds_path,species_path, n, comb, bgn, path, filename, spn))
    
    stopCluster(cl)
    gc()
    
    return(unlist(o))
    
  } else {
    ###sapply(cl,1:length(comb),.wf_ER)
    sapply(1:length(comb), function(i) .wf_ER(i, preds_path, species_path, n, comb, bgn, path, filename, spn))
    
  }
  
}

##################
#project the models
##################

.corePWF <- function(i, pr_new_list, path, filename, spn) {
  print(paste0(">>> Running .corePWF for combination i = ", i))
  
  model_file <- paste0(path, '/models/M_', i, '_', spn, '_', filename, '.sdm')
  if (!file.exists(model_file)) {
    print(paste0("❌ Model file NOT FOUND for COMB-", i, " → Skipping!"))
    return(list(AUC = NA, en = NULL))  # Return empty values if model is missing
  }
  
  m <- read.sdm(model_file)
  print(paste0("Model loaded for COMB-", i))
  
  stat_file <- paste0(path, "/ensembles/statistics/st_COMB-", i, "_", spn, "_", filename, ".rds")
  if (!file.exists(stat_file)) {
    print("Calculating statistics")
    e <- .getEval(m, n = 'species', setting = list(method = 'weighted', stat = 'auc'))
    if (!dir.exists(paste0(path, "/ensembles/statistics/"))) {
      dir.create(paste0(path, "/ensembles/statistics/st_COMB-", i, "_", spn, "_", filename, ".rds"))}
    saveRDS(e, paste0(path, "/ensembles/statistics/st2_COMB-", i, "_", spn, "_", filename, ".rds"))
  } else {
    e <- readRDS(stat_file)
  }
  print(paste0("✅ Model statistics for COMB-", i))
  .AUC <- e@statistics$AUC  # Extract AUC
  
  ### for (s in seq_along(pr_new_list)) {
  s <- 1
  scenario_name <- names(pr_new_list)
  .pr_new <- rast(pr_new_list[[s]])
  
  prediction_file <- paste0(path, "/predicts/projections/proj_pr_", scenario_name, "_COMB-", i, "_", spn, "_", filename, ".tif")
  print(paste0("🔹 Processing: ", scenario_name, " → COMB-", i))
  #if (file.exists(prediction_file)) {
  #print(paste0("⚠️ File already exists, skipping: ", prediction_file))
  #} else {
  tryCatch({
    print("Running prediction...")
    .pr<-predict(m, .pr_new, filename = prediction_file, overwrite = TRUE)
    print(paste0("✅ Prediction complete for ", scenario_name, " (COMB-", i, ")"))
  }, error = function(e) {
    print(paste0("❌ ERROR in COMB-", i, " (Scenario: ", scenario_name, "): ", e$message))
  })
  # }
  
  ensemble_file <- paste0(path, "/ensembles/projections/proj_en_", scenario_name, "_COMB-", i, "_", spn, "_", filename, ".tif")
  en <- NULL
  #if (file.exists(ensemble_file)) {en<-rast(ensemble_file)} else{
  tryCatch({
    print("Running ensemble...")
    en <- ensemble(m, .pr, setting = list(method = 'weighted', stat = 'auc'), filename = ensemble_file, overwrite = TRUE)
    print(paste0("✅ Ensemble complete for ", scenario_name, " (COMB-", i, ")"))
  }, error = function(e) {
    print(paste0("❌ ERROR in COMB-", i, " (Scenario: ", scenario_name, "): ", e$message))
  })
  #}
  ### }
  
  rm(m)
  gc()
  
  return(list(AUC = .AUC, en = en))  # Return AUC and ensemble raster
}

.pwf_ER <- function(i, pr_new_list, comb, path, filename, spn) {
  e <- try(.corePWF(i, pr_new_list, path, filename, spn), silent = TRUE)
  if (inherits(e, 'try-error')) return(list(AUC = NA, en = NULL))
  return(e)
}

runPWF <- function(pr_new_list, comb, path = '.', filename, ncore = 8, spn) {
  if (!is.null(ncore)) {
    ncore <- min(parallel::detectCores() - 2, ncore)
  }
  
  .auc <- rep(NA, length(comb))  # Create empty AUC vector
  e_stack <- list()  # List to store ensemble rasters
  
  if (!is.null(ncore) && ncore > 1) {
    library(parallel)
    cl <- makeCluster(ncore)
    
    clusterEvalQ(cl, {
      library(terra)
      library(sdm)
      NULL
    })
    
    clusterExport(cl, varlist = c("path", "filename", "spn", "comb", "pr_new_list", ".corePWF", ".pwf_ER"), envir = environment())
    
    results <- parLapply(cl, 1:length(comb), function(i) .pwf_ER(i, pr_new_list, comb, path, filename, spn))
    
    stopCluster(cl)
    gc()
    
  } else {
    results <- lapply(1:length(comb), function(i) .pwf_ER(i, pr_new_list, comb, path, filename, spn))
  }
  
  # Collect AUC values and ensemble rasters
  for (i in seq_along(comb)) {
    .auc[i] <- results[[i]]$AUC
    if (!is.null(results[[i]]$en)) {
      e_stack[[i]] <- results[[i]]$en
    }
  }
  
  # Convert e_stack list to a SpatRaster stack
  e_stack <- rast(e_stack)
  
  # Normalize AUC values to sum to 1
  .auc <- .auc / sum(.auc, na.rm = TRUE)
  
  # Compute final ensemble weighted mean
  enFinal <- weighted.mean(e_stack, .auc, na.rm = TRUE)
  
  # Save final ensemble raster
  ensemble_file <- paste0(path, "/ensembles/projections/proj_en_", scenario_name, "_COMB-all_", spn, "_", filename, ".tif")
  writeRaster(enFinal, ensemble_file, overwrite = TRUE)
  
  return(list(AUC = .auc, enFinal = enFinal))
}


