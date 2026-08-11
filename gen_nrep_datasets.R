#===============================================================================
# Generate 100 Replicated Datasets for Three Simulation Scenarios
#
# Scenario 1: No ZCTA-level spatial variation
# Scenario 2: Low ZCTA-level spatial variation
# Scenario 3: High ZCTA-level spatial variation
#===============================================================================

rm(list = ls())

#-------------------------------------------------------------------------------
# File paths
#-------------------------------------------------------------------------------

simulation_dir <- paste0(
  "C:/Users/chech/OneDrive - Dartmouth College/",
  "Research/HeterIndex/Paper Simulation/"
)

# Load the spatial and nonspatial data-generation functions
source(
  file.path(simulation_dir, "spheter_data_gen_fn.R")
)

source(
  file.path(simulation_dir, "heter_data_gen_fn.R")
)

# Replicate identifiers
nrep <- 1:100


#===============================================================================
# Scenario 1: No spatial variation
#===============================================================================

out_dir_sce1 <- file.path(
  simulation_dir,
  "simdata_scenario1"
)

# Create the output directory if it does not already exist
dir.create(
  out_dir_sce1,
  recursive = TRUE,
  showWarnings = FALSE
)

for (s in nrep) {
  
  # Generate the sth replicated dataset
  sim_s <- generate_heter_data(
    seed = s
  )
  
  # Construct the output filename
  output_file <- file.path(
    out_dir_sce1,
    sprintf("sim_data_%03d.rds", s)
  )
  
  # Save the complete simulated-data object
  saveRDS(
    sim_s,
    file = output_file
  )
  
  message(
    "Scenario 1: saved replicate ",
    s,
    " of ",
    length(nrep),
    " to ",
    output_file
  )
}


#===============================================================================
# Scenario 2: Low spatial variation
#
# The spatial random effects are generated as:
#
#   delta ~ MVN{0, tau_delta^(-1) R(rho)}
#
# Here, tau_delta = 0.50 corresponds to a marginal spatial variance of
# approximately 1 / 0.50 = 2, before applying the sum-to-zero constraint.
#===============================================================================

out_dir_sce2 <- file.path(
  simulation_dir,
  "simdata_scenario2"
)

dir.create(
  out_dir_sce2,
  recursive = TRUE,
  showWarnings = FALSE
)

for (s in nrep) {
  
  sim_s <- generate_spatial_heter_data(
    seed   = s,
    rho_km = 30,
    taud   = 0.50
  )
  
  output_file <- file.path(
    out_dir_sce2,
    sprintf("sim_data_%03d.rds", s)
  )
  
  saveRDS(
    sim_s,
    file = output_file
  )
  
  message(
    "Scenario 2: saved replicate ",
    s,
    " of ",
    length(nrep),
    " to ",
    output_file
  )
}


#===============================================================================
# Scenario 3: High spatial variation
#
# Here, tau_delta = 0.25 corresponds to a marginal spatial variance of
# approximately 1 / 0.25 = 4, before applying the sum-to-zero constraint.
# The larger rho also produces more persistent spatial correlation.
#===============================================================================

out_dir_sce3 <- file.path(
  simulation_dir,
  "simdata_scenario3"
)

dir.create(
  out_dir_sce3,
  recursive = TRUE,
  showWarnings = FALSE
)

for (s in nrep) {
  
  sim_s <- generate_spatial_heter_data(
    seed   = s,
    rho_km = 50,
    taud   = 0.25
  )
  
  output_file <- file.path(
    out_dir_sce3,
    sprintf("sim_data_%03d.rds", s)
  )
  
  saveRDS(
    sim_s,
    file = output_file
  )
  
  message(
    "Scenario 3: saved replicate ",
    s,
    " of ",
    length(nrep),
    " to ",
    output_file
  )
}


#===============================================================================
# Optional checks
#===============================================================================

# Load replicate 60 from each scenario
simdat_sce1 <- readRDS(
  file.path(out_dir_sce1, "sim_data_060.rds")
)

simdat_sce2 <- readRDS(
  file.path(out_dir_sce2, "sim_data_060.rds")
)

simdat_sce3 <- readRDS(
  file.path(out_dir_sce3, "sim_data_060.rds")
)

# Examine the true parameter values
simdat_sce1$truth
simdat_sce2$truth
simdat_sce3$truth