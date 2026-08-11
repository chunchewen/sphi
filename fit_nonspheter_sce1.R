#===============================================================================
# Scenario 1: No Spatial Dependence
#
# Fit the nonspatial Heterogeneity-Index model to 100 simulated datasets.
#
# Input:
#   simdata_scenario1/sim_data_001.rds, ..., sim_data_100.rds
#
# Output:
#   mcmc_nonsp_scenario1/mcmc_fit_001.rds, ..., mcmc_fit_100.rds
#===============================================================================

rm(list = ls())


#-------------------------------------------------------------------------------
# 1. Define file paths
#-------------------------------------------------------------------------------

simulation_dir <- paste0(
  "C:/Users/chech/OneDrive - Dartmouth College/",
  "Research/HeterIndex/Paper Simulation/"
)

data_dir <- file.path(
  simulation_dir,
  "simdata_scenario1"
)

mcmc_dir <- file.path(
  simulation_dir,
  "mcmc_nonsp_scenario1"
)

# Create the output directory if it does not already exist
dir.create(
  mcmc_dir,
  showWarnings = FALSE,
  recursive    = TRUE
)


#-------------------------------------------------------------------------------
# 2. Load the nonspatial MCMC function
#-------------------------------------------------------------------------------

source(
  file.path(
    simulation_dir,
    "mcmc_fn_heter.R"
  )
)


#-------------------------------------------------------------------------------
# 3. Specify the simulation replicates and MCMC settings
#-------------------------------------------------------------------------------

nrep <- 1:100

nsim <- 5000
burn <- 2500
thin <- 5

# Number of retained posterior draws
n_save <- floor((nsim - burn) / thin)

message(
  "Each fitted model will retain approximately ",
  n_save,
  " posterior draws."
)


#-------------------------------------------------------------------------------
# 4. Fit the nonspatial model to each simulated dataset
#-------------------------------------------------------------------------------

for (s in nrep) {
  
  cat("\n")
  cat("============================================================\n")
  cat("Fitting Scenario 1 replicate ", s, " of ", length(nrep), "\n", sep = "")
  cat("============================================================\n")
  
  input_file <- file.path(
    data_dir,
    sprintf("sim_data_%03d.rds", s)
  )
  
  output_file <- file.path(
    mcmc_dir,
    sprintf("mcmc_fit_%03d.rds", s)
  )
  
  #---------------------------------------------------------------------------
  # Optional: skip replicates that have already been fitted
  #---------------------------------------------------------------------------
  
  if (file.exists(output_file)) {
    
    message(
      "Replicate ",
      s,
      " has already been fitted. Skipping."
    )
    
    next
  }
  
  #---------------------------------------------------------------------------
  # Check whether the simulated dataset exists
  #---------------------------------------------------------------------------
  
  if (!file.exists(input_file)) {
    
    warning(
      "Input file does not exist for replicate ",
      s,
      ": ",
      input_file
    )
    
    next
  }
  
  #---------------------------------------------------------------------------
  # Fit the model
  #---------------------------------------------------------------------------
  
  start_time <- Sys.time()
  
  tryCatch({
    
    # Read the sth simulated dataset
    sim_s <- readRDS(input_file)
    
    # Set a reproducible seed for the MCMC sampler
    # This is distinct from the seed used to generate the dataset.
    set.seed(10000 + s)
    
    # Run the nonspatial Heterogeneity-Index model
    fit_s <- run_heter_mcmc(
      sim_s$dat,
      nsim = nsim,
      burn = burn,
      thin = thin
    )
    
    # Save the complete fitted-model object
    saveRDS(
      fit_s,
      file = output_file
    )
    
    elapsed_time <- difftime(
      Sys.time(),
      start_time,
      units = "mins"
    )
    
    message(
      "Successfully saved replicate ",
      s,
      " to:\n",
      output_file
    )
    
    message(
      "Elapsed time: ",
      round(as.numeric(elapsed_time), 2),
      " minutes."
    )
    
    # Remove large objects before proceeding to the next replicate
    rm(sim_s, fit_s)
    invisible(gc())
    
  }, error = function(e) {
    
    warning(
      "Replicate ",
      s,
      " failed: ",
      conditionMessage(e)
    )
    
  })
}


#===============================================================================
# 5. Check the completed MCMC fits
#===============================================================================

mcmc_files <- list.files(
  mcmc_dir,
  pattern   = "^mcmc_fit_[0-9]{3}\\.rds$",
  full.names = TRUE
)

message(
  "Number of completed MCMC fits: ",
  length(mcmc_files),
  " out of ",
  length(nrep)
)


#-------------------------------------------------------------------------------
# 6. Load selected fitted models for inspection
#-------------------------------------------------------------------------------

mcmcdat1 <- readRDS(
  file.path(
    mcmc_dir,
    "mcmc_fit_001.rds"
  )
)

mcmcdat3 <- readRDS(
  file.path(
    mcmc_dir,
    "mcmc_fit_003.rds"
  )
)