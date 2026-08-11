#===============================================================================
# Scenario 2: Low Spatial Dependence
#
# Data-generating model: Spatial SpHI model
# Fitted model:          Spatial SpHI model
#
# Purpose:
# Evaluate the performance of the correctly specified spatial model under
# relatively low ZCTA-level spatial variation.
#===============================================================================

rm(list = ls())


#-------------------------------------------------------------------------------
# 1. Define directories
#-------------------------------------------------------------------------------

simulation_dir <- paste0(
  "C:/Users/chech/OneDrive - Dartmouth College/",
  "Research/HeterIndex/Paper Simulation/"
)

data_dir <- file.path(
  simulation_dir,
  "simdata_scenario2"
)

mcmc_dir <- file.path(
  simulation_dir,
  "mcmc_sp_scenario2"
)

dir.create(
  mcmc_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


#-------------------------------------------------------------------------------
# 2. Load the spatial-model MCMC function
#-------------------------------------------------------------------------------

source(
  file.path(
    simulation_dir,
    "mcmc_fn_spheter.R"
  )
)


#-------------------------------------------------------------------------------
# 3. Specify replicates and MCMC settings
#-------------------------------------------------------------------------------

nrep <- 1:100

nsim <- 5000
burn <- 2500
thin <- 5

# Number of retained posterior draws
n_save <- floor((nsim - burn) / thin)

message(
  "Each fit will retain ",
  n_save,
  " posterior draws."
)


#-------------------------------------------------------------------------------
# 4. Fit the spatial SpHI model
#-------------------------------------------------------------------------------

for (s in nrep) {
  
  cat("\n")
  cat("============================================================\n")
  cat("Running Scenario 2 replicate ", s, " of 100\n", sep = "")
  cat("============================================================\n")
  
  input_file <- file.path(
    data_dir,
    sprintf("sim_data_%03d.rds", s)
  )
  
  output_file <- file.path(
    mcmc_dir,
    sprintf("mcmc_fit_%03d.rds", s)
  )
  
  # Skip completed replicates
  if (file.exists(output_file)) {
    message(
      "Replicate ",
      s,
      " has already been fitted. Skipping."
    )
    next
  }
  
  # Check that the simulated dataset exists
  if (!file.exists(input_file)) {
    warning(
      "Input dataset was not found for replicate ",
      s,
      ": ",
      input_file
    )
    next
  }
  
  start_time <- Sys.time()
  
  tryCatch({
    
    # Read data generated under low spatial dependence
    sim_s <- readRDS(input_file)
    
    # Use a reproducible MCMC seed
    set.seed(50000 + s)
    
    # Fit the spatial Heterogeneity-Index model
    fit_s <- run_spheter_mcmc(
      sim_s$dat,
      nsim = nsim,
      burn = burn,
      thin = thin
    )
    
    # Save the fitted-model object
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
      ". Elapsed time: ",
      round(as.numeric(elapsed_time), 2),
      " minutes."
    )
    
    rm(sim_s, fit_s)
    invisible(gc())
    
  }, error = function(e) {
    
    warning(
      "Replicate ",
      s,
      " failed: ",
      conditionMessage(e)
    )
    
    if (exists("sim_s", inherits = FALSE)) {
      rm(sim_s)
    }
    
    if (exists("fit_s", inherits = FALSE)) {
      rm(fit_s)
    }
    
    invisible(gc())
  })
}


#===============================================================================
# 5. Check completed fits
#===============================================================================

mcmc_files <- list.files(
  mcmc_dir,
  pattern = "^mcmc_fit_[0-9]{3}\\.rds$",
  full.names = TRUE
)

message(
  "Completed Scenario 2 spatial-model fits: ",
  length(mcmc_files),
  " out of 100."
)


#-------------------------------------------------------------------------------
# 6. Load selected fits for inspection
#-------------------------------------------------------------------------------

mcmcdat1 <- readRDS(
  file.path(mcmc_dir, "mcmc_fit_001.rds")
)

mcmcdat3 <- readRDS(
  file.path(mcmc_dir, "mcmc_fit_003.rds")
)

mcmcdat5 <- readRDS(
  file.path(mcmc_dir, "mcmc_fit_005.rds")
)