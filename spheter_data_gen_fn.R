#===============================================================================
# Generate Data from the Spatial Heterogeneity-Index Model
#
# This function simulates binary patient-level outcomes under a hierarchical
# logistic regression model with:
#
#   1. Patient-level fixed effects
#   2. HRR-level multivariate random effects
#   3. ZCTA-level spatial random effects
#
# The function also calculates the true HRR-specific heterogeneity index
# across insurance groups.
#===============================================================================

generate_spatial_heter_data <- function(
    seed             = NULL,
    dir_hrr          = paste0(
      "C:/Users/chech/OneDrive - Dartmouth College/",
      "Research/HeterIndex/Simulation/"
    ),
    J                = 19,      # Number of HRRs
    K                = 983,     # Number of ZCTAs
    L                = 4,       # Number of insurance groups
    nis_range        = c(500, 1000),  # Range of HRR-specific sample sizes
    insurance_prob   = c(0.35, 0.30, 0.20, 0.15),
    
    # Fixed effects associated with X:
    # intercept, x2, x3, x4, w1, x2*w1, x3*w1, x4*w1
    beta = c(
      -1.00, -0.25, -0.50, -0.45,
      0.15,  0.25,  0.35,  0.25
    ),
    
    # Fixed effects associated with the additional covariates:
    # w2, w3, and w4
    alpha = c(0.25, 0.50, 0.25),
    
    # Covariance matrix for the HRR-level random effects:
    # theta1 = shared random intercept
    # theta2-theta4 = insurance-specific deviations from the reference group
    sigma = matrix(
      c(
        1.0, -0.3, -0.6, -0.2,
        -0.3,  0.6,  0.2,  0.6,
        -0.6,  0.2,  1.2,  0.5,
        -0.2,  0.6,  0.5,  1.5
      ),
      nrow = L,
      ncol = L,
      byrow = TRUE
    ),
    
    rho_km             = 80,    # Spatial correlation range in kilometers
    taud               = 0.25,  # Spatial precision parameter
    year_map           = 2010,
    state_map          = "FL",
    return_map_objects = FALSE
) {
  
  #-----------------------------------------------------------------------------
  # Load required packages
  #-----------------------------------------------------------------------------
  
  require(mvtnorm)
  require(tigris)
  require(dplyr)
  require(sf)
  require(geosphere)
  
  # Set the random-number seed to ensure reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  #-----------------------------------------------------------------------------
  # Obtain the Florida ZCTA map
  #-----------------------------------------------------------------------------
  
  zipcode <- zctas(
    state = state_map,
    year  = year_map
  )
  
  # Sort ZCTAs by their numeric codes and create sequential internal identifiers
  zipcode <- zipcode[order(as.numeric(zipcode$ZCTA5CE10)), ]
  zipcode$zipcode    <- as.numeric(zipcode$ZCTA5CE10)
  zipcode$zipcode_id <- seq_len(nrow(zipcode))
  
  #-----------------------------------------------------------------------------
  # Link each ZCTA to its Dartmouth Atlas HRR
  #-----------------------------------------------------------------------------
  
  dhrr <- read.csv(
    file.path(dir_hrr, "hrr_fl.csv")
  )
  
  zipcode <- left_join(
    zipcode,
    dhrr[, c("zipcode", "hrrnum_id")],
    by = "zipcode"
  )
  
  # Retain only ZCTAs that can be assigned to an HRR
  zipcode <- zipcode[!is.na(zipcode$hrrnum_id), ]
  zipcode <- zipcode[order(zipcode$zipcode_id), ]
  
  # Update K if the number of mapped ZCTAs differs from the supplied value
  K_use <- nrow(zipcode)
  
  if (K_use != K) {
    message(
      "K in function call = ", K,
      ", but mapped ZCTAs = ", K_use,
      ". Using K = ", K_use, " from merged map."
    )
    
    K <- K_use
  }
  
  #-----------------------------------------------------------------------------
  # Generate HRR-specific sample sizes
  #-----------------------------------------------------------------------------
  
  # nis[j] is the number of simulated patients in HRR j
  nis <- sample(
    nis_range[1]:nis_range[2],
    size    = J,
    replace = TRUE
  )
  
  # Total sample size
  N <- sum(nis)
  
  # HRR identifier for each patient
  id <- rep(seq_len(J), nis)
  
  #-----------------------------------------------------------------------------
  # Generate insurance-group membership
  #-----------------------------------------------------------------------------
  
  # Insurance groups:
  #   1 = FFS Medicare (reference)
  #   2 = Medicare Advantage
  #   3 = Private insurance
  #   4 = Uninsured
  xtype <- sample(
    seq_len(L),
    size    = N,
    replace = TRUE,
    prob    = insurance_prob
  )
  
  # Dummy variables for the non-reference insurance groups
  x2 <- as.integer(xtype == 2)
  x3 <- as.integer(xtype == 3)
  x4 <- as.integer(xtype == 4)
  
  #-----------------------------------------------------------------------------
  # Generate additional patient-level covariates
  #-----------------------------------------------------------------------------
  
  # Continuous covariate
  w1 <- rnorm(N)
  
  # Binary covariate
  w2 <- rbinom(N, size = 1, prob = 0.5)
  
  # Three-category covariate with category A as the reference group
  wx <- sample(
    c("A", "B", "C"),
    size    = N,
    replace = TRUE
  )
  
  w3 <- as.integer(wx == "B")
  w4 <- as.integer(wx == "C")
  
  #-----------------------------------------------------------------------------
  # Construct the fixed-effect design matrices
  #-----------------------------------------------------------------------------
  
  # X contains:
  #   1. Intercept
  #   2. Insurance-group indicators
  #   3. Continuous covariate w1
  #   4. Insurance-by-w1 interaction terms
  X <- cbind(
    1,
    x2,
    x3,
    x4,
    w1,
    x2 * w1,
    x3 * w1,
    x4 * w1
  )
  
  # W contains the remaining patient-level covariates
  W <- cbind(
    w2,
    w3,
    w4
  )
  
  p <- ncol(X)  # Number of columns in X: 8
  q <- ncol(W)  # Number of columns in W: 3
  
  #-----------------------------------------------------------------------------
  # Generate HRR-level multivariate random effects
  #
  # theta_mat is a J x L matrix:
  #
  #   Column 1: Shared HRR random intercept
  #   Column 2: HRR-specific deviation for insurance group 2
  #   Column 3: HRR-specific deviation for insurance group 3
  #   Column 4: HRR-specific deviation for insurance group 4
  #-----------------------------------------------------------------------------
  
  theta_mat <- mvtnorm::rmvnorm(
    n     = J,
    sigma = sigma
  )
  
  # Impose a sum-to-zero constraint separately on each random-effect component
  theta_mat <- sweep(
    theta_mat,
    MARGIN = 2,
    STATS  = colMeans(theta_mat),
    FUN    = "-"
  )
  
  theta1 <- theta_mat[, 1]
  theta2 <- theta_mat[, 2]
  theta3 <- theta_mat[, 3]
  theta4 <- theta_mat[, 4]
  
  # Expand HRR-level random effects to the patient level
  Theta1 <- theta1[id]
  Theta2 <- theta2[id]
  Theta3 <- theta3[id]
  Theta4 <- theta4[id]
  
  #-----------------------------------------------------------------------------
  # Assign patients to ZCTAs nested within their HRRs
  #-----------------------------------------------------------------------------
  
  # Create a list containing the ZCTAs within each HRR
  zctas_by_hrr <- split(
    zipcode$zipcode_id,
    zipcode$hrrnum_id
  )
  
  # Initialize the patient-level ZCTA identifiers
  zcta_id_i <- integer(N)
  
  # Randomly assign each patient to a ZCTA within the patient's HRR
  for (j in seq_len(J)) {
    
    idx <- which(id == j)
    
    zcta_id_i[idx] <- sample(
      zctas_by_hrr[[as.character(j)]],
      size    = length(idx),
      replace = TRUE
    )
  }
  
  # Number of simulated patients in each ZCTA
  nis.z <- tabulate(zcta_id_i, nbins = K)
  
  #-----------------------------------------------------------------------------
  # Generate ZCTA-level spatial random effects
  #-----------------------------------------------------------------------------
  
  # Obtain the geographic centroid of each ZCTA
  cent   <- suppressWarnings(st_centroid(zipcode))
  coords <- st_coordinates(cent)
  
  # Calculate pairwise great-circle distances between ZCTA centroids
  dist_km <- geosphere::distm(
    coords,
    fun = geosphere::distHaversine
  ) / 1000
  
  # Exponential spatial correlation matrix:
  # R[k,m] = exp(-distance[k,m] / rho)
  #
  # A small diagonal term is added for numerical stability
  R <- exp(-dist_km / rho_km) + diag(0.0001, K)
  
  # Generate the ZCTA-level spatial random effects:
  # delta ~ MVN(0, tau_delta^(-1) R)
  delta_gen <- c(
    mvtnorm::rmvnorm(
      n     = 1,
      sigma = (1 / taud) * R
    )
  )
  
  # Impose a sum-to-zero constraint
  delta <- delta_gen - mean(delta_gen)
  
  #-----------------------------------------------------------------------------
  # Generate the binary outcome
  #
  # For patient i, the linear predictor is:
  #
  # eta_i =
  #   X_i beta + W_i alpha
  #   + Theta1_i
  #   + Theta2_i x2_i
  #   + Theta3_i x3_i
  #   + Theta4_i x4_i
  #   + delta_{k(i)}
  #
  # where k(i) denotes the ZCTA containing patient i.
  #-----------------------------------------------------------------------------
  
  eta <- as.vector(
    X %*% beta +
      W %*% alpha +
      Theta1 +
      Theta2 * x2 +
      Theta3 * x3 +
      Theta4 * x4 +
      delta[zcta_id_i]
  )
  
  # Convert the linear predictor to a probability using the logistic link
  pr <- plogis(eta)
  
  # Generate the binary outcome
  y <- rbinom(
    n    = N,
    size = 1,
    prob = pr
  )
  
  #-----------------------------------------------------------------------------
  # Calculate the implied covariance among insurance-specific HRR effects
  #
  # The original random effects are:
  #
  #   theta_j = (theta_j1, theta_j2, theta_j3, theta_j4)'
  #
  # The insurance-specific HRR effects are:
  #
  #   FFS       = theta_j1
  #   MA        = theta_j1 + theta_j2
  #   Private   = theta_j1 + theta_j3
  #   Uninsured = theta_j1 + theta_j4
  #
  # Therefore:
  #
  #   Sigma* = A Sigma A'
  #-----------------------------------------------------------------------------
  
  A_mat <- matrix(
    c(
      1, 0, 0, 0,
      1, 1, 0, 0,
      1, 0, 1, 0,
      1, 0, 0, 1
    ),
    nrow  = L,
    ncol  = L,
    byrow = TRUE
  )
  
  sigma_star  <- A_mat %*% sigma %*% t(A_mat)
  cormat_star <- cov2cor(sigma_star)
  
  insurance_names <- c(
    "FFS",
    "MA",
    "Private",
    "Uninsured"
  )
  
  rownames(sigma_star) <- colnames(sigma_star) <- insurance_names
  rownames(cormat_star) <- colnames(cormat_star) <- insurance_names
  
  #-----------------------------------------------------------------------------
  # Calculate the true HRR- and insurance-specific linear predictors
  #
  # For HRR j:
  #
  #   FFS:
  #     Delta_j1 = beta_1 + theta_j1
  #
  #   MA:
  #     Delta_j2 = beta_1 + beta_2 + theta_j1 + theta_j2
  #
  #   Private:
  #     Delta_j3 = beta_1 + beta_3 + theta_j1 + theta_j3
  #
  #   Uninsured:
  #     Delta_j4 = beta_1 + beta_4 + theta_j1 + theta_j4
  #-----------------------------------------------------------------------------
  
  Delta_true <- matrix(
    NA_real_,
    nrow = J,
    ncol = L
  )
  
  Delta_true[, 1] <-
    beta[1] +
    theta_mat[, 1]
  
  Delta_true[, 2] <-
    beta[1] + beta[2] +
    theta_mat[, 1] + theta_mat[, 2]
  
  Delta_true[, 3] <-
    beta[1] + beta[3] +
    theta_mat[, 1] + theta_mat[, 3]
  
  Delta_true[, 4] <-
    beta[1] + beta[4] +
    theta_mat[, 1] + theta_mat[, 4]
  
  colnames(Delta_true) <- insurance_names
  
  # The true H-index is the standard deviation of the insurance-specific
  # linear predictors within each HRR
  H_true <- apply(
    Delta_true,
    MARGIN = 1,
    FUN    = sd
  )
  
  # Rank HRRs from the largest to the smallest H-index
  rank_H_true <- rank(
    -H_true,
    ties.method = "average"
  )
  
  #-----------------------------------------------------------------------------
  # Assemble the returned object
  #-----------------------------------------------------------------------------
  
  out <- list(
    
    # Simulated data and dimension information
    dat = list(
      y           = y,
      X           = X,
      W           = W,
      id          = id,
      xtype       = xtype,
      zcta_id_i   = zcta_id_i,
      nis         = nis,
      nis.z       = nis.z,
      N           = N,
      J           = J,
      K           = K,
      L           = L,
      p           = p,
      q           = q
    ),
    
    # True parameter values used to generate the data
    truth = list(
      beta         = beta,
      alpha        = alpha,
      sigma        = sigma,
      sigma_star   = sigma_star,
      cormat_star  = cormat_star,
      theta_mat    = theta_mat,
      delta        = delta,
      taud         = taud,
      rho_km       = rho_km,
      Delta_true   = Delta_true,
      H_true       = H_true,
      rank_H_true  = rank_H_true
    )
  )
  
  # Optionally return the spatial map and distance objects
  if (return_map_objects) {
    out$map <- list(
      zipcode = zipcode,
      coords  = coords,
      dist_km = dist_km,
      R       = R
    )
  }
  
  return(out)
}