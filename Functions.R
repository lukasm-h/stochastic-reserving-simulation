############# helper functions ################

## Technical helper

# Apply absolute value if needed (for comparison)

toggle_comp <- function(x, a, stat){
  if (a == "absolute") {
    x <- abs(x)
  }
  
  if(stat == "mean") {
    x <- mean(x, na.rm = TRUE)
  }
  
  if(stat == "sd") {
    x <- sd(x, na.rm = TRUE)
  }
  x
}

# Integrate Finney-Function (1941, bias correction for LN fitting)
g_k <- function(t,k){
  value = 0
 n <- 0:9
 
 product_terms <- cumprod(k+ 2*n)
 summands <- ((k^(2*n) * (k + 2*n)) / ((k+1)^n * product_terms)) * (t^n / factorial(n))
 
 sum(summands)
}


## Simulation_routine: picking the right function by index

bootstrap_result <- function(triangle, model = c("Mack", "ODP", "LN"), B){
  switch(
    model,
    Mack = bootstrap_mack(triangle, B),
    ODP = bootstrap_odp(triangle, B),
    LN = bootstrap_LN(triangle, B)
  )
}
  
analytical_result <- function(triangle, model = c("Mack", "ODP", "LN")){
  switch(
    model,
    Mack = analytic_mack(triangle),
    ODP = {
      inc <- cum2incr(as.triangle(triangle))
      
      # Count where the fallback is needed
      fallback_columns <- which(colSums(inc < 0, na.rm = TRUE) > 0)
      
      if(length(fallback_columns) > 0) {
       # message("ODP fallback was used.")
        result <- analytic_odp_fallback(triangle)
      } else {
        result <- analytic_odp(triangle)}
      result$fallback_columns <- fallback_columns
      result
      },
    LN = analytic_ln(triangle)
    )
}

known_parameters_result <- function(triangle, parameters, model = c("Mack", "ODP", "LN"), nb_sims){
  switch(
    model,
    Mack = results_mack_known_parameters(triangle, parameters$cl_factors, parameters$cl_sigmas, parameters$residuals, nb_sims),
    ODP = results_odp_known_parameters(parameters$alpha_hat, parameters$beta_hat, parameters$dispersion, nb_sims),
    LN = results_ln_known_parameters(triangle, parameters$mus, parameters$sigma_sqs, nb_sims)
  )
}

sim_triangle <- function(model, parameters){
  switch(
    model,
    Mack = sim_triangle_Mack(parameters$cl_factors, parameters$cl_sigmas, parameters$residuals, parameters$start_values)$rectangle,
    ODP = sim_triangle_ODP(parameters$alpha_hat, parameters$beta_hat, parameters$dispersion),
    LN = sim_triangle_LN(parameters$mus, parameters$sigma_sqs, parameters$start_values)
  )
}
## Empirical-Distributional work

# Get standardized residuals

standardize_residuals <- function(residuals){
  
  empirical_sd <- sqrt(mean((residuals - mean(residuals))^2))
  standardized_residuals <- (residuals - mean(residuals)) / empirical_sd
  
  standardized_residuals
}


# Summarize Reserves results

reserve_results <- function(reserves){
  list(
    mean = mean(reserves),
    sd = sd(reserves),
    q50 = quantile(reserves, 0.5),
    q65 = quantile(reserves, 0.65),
    q80 = quantile(reserves, 0.8),
    q95 = quantile(reserves, 0.95)
  )
}




## Triangle work

# Hide future part of triangle
hide_future <- function(data){
  ## Checks
  if (nrow(data) != ncol(data)) {
    stop("Matrix must be symmetric!")
  } 
  
  data[row(data) + col(data) > nrow(data) + 1] <- NA_real_
  data
}

## Mack-helpers

# Get the Mack-Residuals with variance correction
mack_residuals <- function(triangle){
  # Check
  if (any(triangle <= 0, na.rm = TRUE)) {
    stop("Values must be positive!")
  }
  
  N <- nrow(triangle)
  mack_model <- MackChainLadder(triangle, est.sigma="Mack")
  
  # Only up to N-2 defined residuals, index all the relevant terms to 1 to N-2
  periods <- seq_len(N-2)
  elements_per_period <- N - periods
  
  
  # Create the relevant indices from the full triangle to avoid slower matrix calculation
  relevant_AY <- sequence(elements_per_period)
  relevant_dev <- rep.int(periods, elements_per_period)
  
  # Construct C_i,j and Ci,j+1 (shifted) vectors from the relevant matrix values
  cumulative <- triangle[cbind(relevant_AY, relevant_dev)]
  cumulative_shifted <- triangle[cbind(relevant_AY, relevant_dev+1)]
  
  # Get right number of replications of the relevant Mack-sizes for the residual calculation
  f_rep <- mack_model$f[relevant_dev]
  sigma_rep <- mack_model$sigma[relevant_dev]
  elements_per_period_rep <- elements_per_period[relevant_dev]
  
  residuals <- sqrt(elements_per_period_rep / (elements_per_period_rep -1)) * 
    (cumulative_shifted - f_rep * cumulative) / (sigma_rep * sqrt(cumulative))
}


################# 0. Analytical results #################

## 1. Mack - Analytical

# Analytic calculation of Mack estimates
analytic_mack <- function(triangle){
  if(any(triangle <= 0, na.rm = TRUE)) {
    warning("Analytic Mack failed, non-positive value in triangle.")
    return(
      list(
        mean = NA_real_,
        msep = NA_real_
      )
    )
  }
  
  mack_model <- MackChainLadder(triangle, est.sigma = "Mack")
  
  list(
    mean = sum(summary(mack_model)$ByOrigin$IBNR),
    msep = mack_model$Total.Mack.S.E
  )
}

# For simulation of Mack-reserves conditional on some upper triangle with known parameters
results_mack_known_parameters <- function(triangle, cl_factors, cl_sigmas, residuals, nb_sims){
  N <- nrow(triangle)
  
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  to_ultimate_factors <- c(1, cumprod(rev(cl_factors)))
  
  conditional_mean_ultimate <- sum(diagonal * to_ultimate_factors)
  
  # Preparing for conditional variance (see Kor. 11)
  after_k_cl_factors <- rev(to_ultimate_factors[1:(N-1)]) 
  sigma_sum_simplified <- c(0, cumsum(rev(cl_sigmas^2 * after_k_cl_factors / cl_factors)))
  
  conditional_variance <- sum(diagonal * to_ultimate_factors * sigma_sum_simplified)
  
  # Empirical Quantiles from simulation with known parameters
  reserves_sim <- numeric(nb_sims)
  
  #Also count resamples, if they were needed (should basically never happen, because models are congruent for this function)
  res_resample_future_values_total <- numeric(N-1)
  
  for (i in seq_len(nb_sims)){
    sim_reserve <- sim_reserve_mack(diagonal, cl_factors, cl_sigmas, residuals)
    reserves_sim[i] <- sim_reserve$reserve_obs 
    res_resample_future_values_total <- res_resample_future_values_total + sim_reserve$res_resample_future_values_bs
  }
  
  sim_results <- reserve_results(reserves_sim)
  
  list(
    mean = conditional_mean_ultimate - sum(diagonal),
    sd = sqrt(conditional_variance),
    MC_mean = sim_results$mean,
    MC_sd = sim_results$sd,
    q50 = sim_results$q50,
    q65 = sim_results$q65,
    q80 = sim_results$q80,
    q95 = sim_results$q95,
    res_resample_future_values = res_resample_future_values_total
  )  
}


## 2. ODP - Analytical

# Analytical ODP results on mean and MSEP
analytic_odp <- function(triangle){
  odp_model <- glmReserve(as.triangle(triangle), mse.method = "formula")
  
  list(
    mean = odp_model$summary["total", "IBNR"],
    msep = odp_model$summary["total", "S.E"],
    failed = integer(0)
  )
}

# ODP - Analytical_fallback
#GLM reserve does only work for positive increments, but the quasi ODP-score should allow for some relaxing as long as positive columns/rows are there. So we need a fallback.
analytic_odp_fallback <- function(triangle){
  
  ## We essentially reproduce glmReserve()-function from the ChainLadder-package, but with the quasi-ODP
  
  # Get incremental values
  inc <- cum2incr(as.triangle(triangle))
  
  # Get dataframe
  data <- as.LongTriangle(inc, varnames = c("origin", "dev"), na.rm = FALSE)
  
  data_past <- data[!is.na(data$value),]
  data_future <- data[is.na(data$value),]
  
  # Check if system of equations is actually solvable (so quasi-score working?)
  bad_rows <- which(tapply(data_past$value, data_past$origin, sum) <= 0)
  bad_columns <- which(tapply(data_past$value, data_past$dev, sum) <= 0)
  
  if(length(bad_rows) > 0) {
    stop("Full row negative!")
  }
  
  if(length(bad_columns) > 0) {
    return(
      list(
        mean = NA_real_,
        msep = NA_real_,
        failed = bad_columns
      )
    )
  }
  
  # Define quasi-ODP
  
  quasi_odp <- quasi(
    link = "log",
    variance = list(
      varfun = function(mu) mu,
      validmu = function(mu){all(is.finite(mu)) && all(mu >0)
        },
      dev.resids = function(y,mu,weight){
        (weight * (y-mu)^2) / mu
      },
      initialize = expression({
        n <- rep.int(1, nobs)
        if (mean(y) <= 0) {
          stop("With mean of the data being non-postive, no ODP fitting is possible.")
        }
        mustart <- rep(mean(y), length(y)) # this is the alternative to the glmReserve-fitting-start, which can crash under the relaxed restrictions
      }),
      name = "Quasi-ODP"
    )
  )
  
  # Now actually fit GLM
  fit <- glm(
    value ~ factor(origin) + factor(dev),
    family = quasi_odp,
    data = data_past
  )
  
  # Get future values
  future_increments <- predict(fit, newdata = data_future, type = "response")
  future_linked_values <- predict(fit, newdata = data_future, type = "link")
  
  # Calculate reserve estimator
  mean_reserve <- sum(future_increments)
  
  # Now for MSEP, we just mirror glmReserve
  
  phi <- with(fit, sum(weights * residuals^2)/df.residual)
  
  msep_process_error <- phi * mean_reserve
  
  design_future <- model.matrix(
    delete.response(terms(fit)),
    data_future,
    xlev = fit$xlevels
  )
  
  Cov <- design_future %*% vcov(fit) %*% t(design_future)
  
  delta <- fit$family$mu.eta(future_linked_values)
  
  msep_estimation_error <- as.numeric(t(delta) %*% Cov %*% delta)
  
  msep <- sqrt(msep_process_error + msep_estimation_error)
  
  # Return as before
  list(
    mean = mean_reserve,
    msep = msep,
    failed = integer(0)
  )
}

# For simulation of (quasi)-ODP-reserves conditional on known parameters (and potentially conditional on upper triangle as well, but that is not important for simulation, since incremental claims are independent)
results_odp_known_parameters <- function(alphas, betas, dispersion, nb_sims){
  
  # Checks
  if (dispersion <= 1) {
    stop("Dispersion must be greater 1 for this simulation function!")
  }
  if (length(alphas) != length(betas)) {
    stop("Must have symmetric parameter number")
  }
  
  N <- length(alphas)
  
  # Since only future is relevant, extract only parameter pair that refer to the future
  
  # How often does each AY contribute to the future
  AY_count <- seq_len(N) - 1
  
  # Replicate the AY for the above calculated times
  future_AY <- rep(
    seq_len(N),
    times = AY_count
  )
  
  # For each of the AYs in the vector above, give the future dev_period
  future_dev <- (N+1) - sequence(1:(N-1))
  
  future_means <- alphas[future_AY] * betas[future_dev]
  future_sizes <- future_means / (dispersion - 1)
  future_variances <- future_means * dispersion
  
  
  # Simulate reserves for empirical quantiles 
  reserves_sim <- numeric(nb_sims)
  for (i in seq_len(nb_sims)){
    
    future_claims <- rnbinom(
      n = length(future_means),
      mu = future_means,
      size = future_sizes
    )
    
    reserves_sim[i] <- sum(future_claims) 
  }
  
  sim_results <- reserve_results(reserves_sim)
  
  list(
    mean = sum(future_means),
    sd = sqrt(sum(future_variances)),
    MC_mean = sim_results$mean,
    MC_sd = sim_results$sd,
    q50 = sim_results$q50,
    q65 = sim_results$q65,
    q80 = sim_results$q80,
    q95 = sim_results$q95
  )  
  
  
}


## 3. LN - Analytical

# Analytical LN results on mean and MSEP
analytic_ln <- function(triangle){
  N <- nrow(triangle)
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  
  # Fit parameters
  parameters_est <- fit_triangle_LN(triangle)
  mus <- parameters_est$mus
  sigma_sqs <- parameters_est$sigma_sqs
  
  # Calculate dev factor estimators with Finney correction
  theta_hat <-  numeric(N-1)
  
  # Create the zero-indexed maximum for easy conversion to written formula
  n_zero <- N-1
  j <- 0:(n_zero-2)
  k_s <- n_zero-j-1
  
  theta_hat[1:(N-2)] <- exp(mus[-length(mus)]) * mapply(FUN = g_k, t = 0.5 * sigma_sqs[-length(sigma_sqs)], k = k_s)
  theta_hat[N-1] <- triangle[1, N] / triangle[1, N-1]
  
  # Calculate msep for individual years
  # Define the estimators first
  theta_hat_2 <- numeric(N-1)
  
  theta_hat_2[1:(N-2)] <- exp(2*mus[-length(mus)]) * mapply(FUN = g_k, t = 2 * sigma_sqs[-length(sigma_sqs)], k = k_s)
  theta_hat_2[N-1] <- (triangle[1, N] / triangle[1, N-1])^2
  
  theta_2_tilde <- numeric(N-1)
  
  theta_2_tilde[1:(N-2)] <- exp(2*mus[-length(mus)]) * mapply(FUN = g_k, t = ((n_zero-j-2)/(n_zero-j-1)) * sigma_sqs[-length(sigma_sqs)], k = k_s)
  theta_2_tilde[N-1] <- exp(-sigma_sqs[length(sigma_sqs)]) * (triangle[1, N] / triangle[1, N-1])^2
  
  # Calculate the actual msep for the individual years
  msep_ind <- numeric(N)
  msep_ind[2:N] <- diagonal[2:N]^2 * (cumprod(rev(theta_hat_2)) - 2* cumprod(rev(theta_2_tilde)) + cumprod(rev(theta_hat^2)))
  
  # Calculate the second (mixed) sum of the msep for all AYs combined
  theta_hat_product <- c(1, cumprod(rev(theta_hat)))
  theta_2_tilde_product <- c(1, cumprod(rev(theta_2_tilde)))

  ultimate_values <- diagonal * theta_hat_product
  larger_i_sum <- rev(cumsum(rev(ultimate_values))) - ultimate_values
  
  mixed_sum <- sum(diagonal * ((theta_hat_product^2 - theta_2_tilde_product) / theta_hat_product) * larger_i_sum)
  
  msep2 <- sum(msep_ind) + 2*mixed_sum
  
  # Return
  list(
    mean = sum(ultimate_values - diagonal),
    msep = sqrt(msep2)
  )
}

# For simulation of LN-reserves conditional on some upper triangle with known parameters
results_ln_known_parameters <- function(triangle, mus, sigma_sqs, nb_sims){
  
  # Preparing relevant quantities
  N <- nrow(triangle)
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  
  # Calculate exp. dev. factors with known parameters
  mean_dev_factors <- exp(mus + 0.5 * sigma_sqs)  
  to_ultimate_factors <- c(1, cumprod(rev(mean_dev_factors)))
  
  conditional_mean_ultimate <- sum(diagonal * to_ultimate_factors)
  
  # Calculate variance of dev factors and thus process variance with known parameters (see Lemma LN 9a)
  var_dev_factors <- exp(2 * mus + sigma_sqs) * (exp(sigma_sqs) - 1)
  
  # See Lemma LN9b
  to_ultimate_var_product_1 <- c(1, cumprod(rev(var_dev_factors + mean_dev_factors^2)))
  to_ultimate_var_product_2 <- c(1, cumprod(rev(mean_dev_factors^2)))
  to_ultimate_var <- to_ultimate_var_product_1 - to_ultimate_var_product_2
  
  conditional_variance_ultimate <- sum(diagonal^2 * to_ultimate_var)
  
  # Empirical Quantiles from simulation with known parameters
  reserves_sim <- numeric(nb_sims)
  for (i in seq_len(nb_sims)){
    reserves_sim[i] <- sim_reserve_LN(diagonal, mus, sigma_sqs) 
  }
  
  sim_results <- reserve_results(reserves_sim)
  
  list(
    mean = conditional_mean_ultimate - sum(diagonal),
    sd = sqrt(conditional_variance_ultimate),
    MC_mean = sim_results$mean,
    MC_sd = sim_results$sd,
    q50 = sim_results$q50,
    q65 = sim_results$q65,
    q80 = sim_results$q80,
    q95 = sim_results$q95
  )  
}

################# A. Fitting on a whole triangle ################

## 1. Mack Triangle

fit_triangle_mack <- function(triangle, emp_norm = FALSE){
  mack_model <- MackChainLadder(triangle, est.sigma = "Mack")
  
  # Get correct residuals
  if (emp_norm) {
    residuals <- standardize_residuals(mack_residuals(triangle))
  } else {
    residuals <- mack_residuals(triangle) # optionally can normalize by empirical moments
  }
  
  # Return
  list(
    cl_factors = mack_model$f[-length(mack_model$f)],
    cl_sigmas = mack_model$sigma,
    residuals = residuals,
    start_values = triangle[,1]
  )
}


# 2. ODP

fit_triangle_ODP <- function(triangle, cum = TRUE){
  ## Checks
  
  if (nrow(triangle) != ncol(triangle)) {
    stop("Matrix must be symmetric!")
  } 
  
  if (any(colSums(cum2incr(triangle), na.rm =TRUE) <= 0)) {
    stop("Values must be non-negative!")
  }
  
  # Fit GLM-model
  N <- nrow(triangle)
  glm_model <- glmReserve(triangle, cum = cum, mse.method = "formula")$model
  
  # Retrieve parameters in log-form
  ay_names <- paste0("factor(origin)", glm_model$xlevels[["factor(origin)"]][-1])
  dev_names <- paste0("factor(dev)", glm_model$xlevels[["factor(dev)"]][-1])
  
  glm_parameters_ay <- numeric(length = N)
  glm_parameters_dev <- numeric(length = N)
  
  glm_parameters_ay[-1] <- glm_model$coefficients[ay_names]
  
  glm_parameters_dev[-1] <- glm_model$coefficients[dev_names]
  
  # Calculate actual parameters 
  alpha_hat <- exp(glm_model$coefficients["(Intercept)"] + glm_parameters_ay)
  beta_hat <- exp(glm_parameters_dev)
  dispersion <- summary(glm_model)$dispersion
  
  # Return
  list(
    alpha_hat = alpha_hat,
    beta_hat = beta_hat,
    dispersion = dispersion
  )
}

## 3. Log-Normal-Triangle 

fit_triangle_LN <- function(triangle) {
  ## Checks
  if (any(triangle <= 0, na.rm = TRUE)) {
    stop("Values must be positive!")
  }
  
  # Calculate individual dev_factors
  log_dev_factors <- log(triangle[, -1, drop = FALSE] 
                         / triangle[,-ncol(triangle), drop = FALSE])
  
  # Calculate parameter estimators
  mus <- colMeans(log_dev_factors, na.rm = TRUE)
  sigma_sqs <- apply(log_dev_factors, MARGIN = 2, FUN = var, na.rm = TRUE)
  
  # Set last variance parameter via Mack approximation
  sigma_sqs[length(sigma_sqs)] <- min(c(sigma_sqs[length(sigma_sqs)-1]^2/sigma_sqs[length(sigma_sqs)-2], sigma_sqs[length(sigma_sqs)-1], sigma_sqs[length(sigma_sqs)-2]))
  
  list(mus = mus, sigma_sqs = sigma_sqs, start_values = triangle[,1])
}



################# B. Simulating ################

### 1. Mack-Model

sim_triangle_Mack <- function(cl_factors, cl_sigmas, residuals, start_values){
  N <- length(start_values)
  
  ## Checks
  if (length(cl_factors) != length(cl_sigmas)) {
    stop("cl_factors and sigmas must same length.")
  }
  
  if (any(cl_sigmas < 0)) {
    stop("sigmas must be non-negative.")
  }
  
  if (any(start_values < 0)) {
    stop("start_values must be non-negative.")
  }
  if (N-1 != length(cl_factors)) {
    stop("Triangle must be symmetric!")
  }
  
  # Sample residuals for all cells after time 0
  sample_residuals <- matrix(
    sample(
      residuals, 
      size = N * (N-1),
      replace = TRUE
      ),
    nrow = N,
    ncol = N-1)
  
  # Create empty Cij-Matrix
  cumulative <- matrix(
    NA_real_,
    nrow = N,
    ncol = N
  )
  
  # Set starting values
  cumulative[, 1] <- start_values
  
  # Prepare counter for potential resampling of triangles
  resample_count_mack_simulation <- numeric(N-1)
  
  # Calculate full claims data
  for (j in seq_len(N-1)) {
    
    cumulative[, j+1] = cumulative[,j] * cl_factors[j] + 
      cl_sigmas[j] * sqrt(cumulative[,j]) * sample_residuals[,j]
    
    # Check for negative cumulative values
    negative <- cumulative[,j+1] <= 0
    
    attempt <- 0
    
    while (any(negative)) {
    
      attempt <- attempt + 1
      
      if(attempt > 50){
        stop(
          "Could not obtain positive cumulative Mack claims after 50 tries."
        )
      }
      
      # Add to the counter
      resample_count_mack_simulation[j] <- resample_count_mack_simulation[j] + sum(negative)
    
      # Redraw exactly the residuals that are responsible for the negative cumulative claims (happens seldom, so just a quick fix)
      sample_residuals[negative, j] <- sample(
        residuals,
        size = sum(negative),
        replace = TRUE
      )
      
      # Then recalculate
      cumulative[negative, j+1] = cumulative[negative,j] * cl_factors[j] + 
        cl_sigmas[j] * sqrt(cumulative[negative,j]) * sample_residuals[negative,j]
      
      negative <- cumulative[,j+1] <= 0
      
      }
  }
  
  # Return
  list(
    rectangle = cumulative,
    resample_count_mack_simulation = resample_count_mack_simulation
    )
}

sim_and_refit_devfactors_mack <- function(means, factors, weights, residuals){
  N <- nrow(means)
  
  # Create and fill matrix with resampled residuals
  resample_residuals <- matrix(
    NA_real_,
    nrow = N,
    ncol = N-1
  )
  
  resample_residuals[row(resample_residuals) + col(resample_residuals) <= N] <- sample(
    residuals, 
    size = (N * (N-1)) / 2,
    replace = TRUE
  )
  
  resample_dev_factors <- means + factors * resample_residuals
  
#  if (any(resample_dev_factors <= 0, na.rm = TRUE)) {
#    warning("Negative individual dev factor detected.")
#  }
  
  cl_factors_refit <- colSums(weights * resample_dev_factors, na.rm = TRUE)
  
  # If there is ever a non-positive fitted CL-factor: Count it and retry:
  
  res_resample_dev_factor_bs <- integer(N - 1)
  
  attempt <- 0
  
  while(any(cl_factors_refit <= 0)){
    attempt <- attempt + 1
    
    if(attempt > 50){
      stop(
        "Could not obtain positive CL-Mack-factor after 50 tries."
      )
    }
    
    bad_cl_dev_periods <- which(cl_factors_refit <= 0)
    
    # Counter
    res_resample_dev_factor_bs[bad_cl_dev_periods] <- res_resample_dev_factor_bs[bad_cl_dev_periods] + 1
    
    # warning(
    #   sprintf(
    #     "Attempt %d: Nonpositive CL-Mack-factor in development periods: %s. Resampling those periods.",
    #     attempt,
    #     paste(bad_cl_dev_periods, collapse = ", ")
    #   )
    # )
    
    for (j in bad_cl_dev_periods){
      individual_dev_factor_count <- seq_len(N-j)
      
      resample_residuals[individual_dev_factor_count, j] <- sample(
        residuals,
        size = length(individual_dev_factor_count),
        replace = TRUE
      )
      
      resample_dev_factors[individual_dev_factor_count, j] <- means[individual_dev_factor_count, j] + factors[individual_dev_factor_count, j] * resample_residuals[individual_dev_factor_count, j]
      
      cl_factors_refit[j] <- sum(weights[individual_dev_factor_count, j] * resample_dev_factors[individual_dev_factor_count, j], na.rm = TRUE)
    }
  }
  
  if (any(cl_factors_refit <= 0)){
    stop("Non-positive fitted CL-Mack-factor detected!")
  }
  # Return
  list(
    cl_factors = cl_factors_refit,
    res_resample_dev_factor_bs = res_resample_dev_factor_bs
  )
  
}

# Sim only lower triangle and return reserves from given diagonal, residuals and cl_factors
sim_reserve_mack <- function(diagonal, cl_factors, cl_sigmas, residuals){
  N <- length(diagonal)
  
  # First draw new residuals for lower triangle
  future_residuals <- sample(
    residuals,
    size = (N * (N-1)) / 2,
    replace = TRUE
  )
  
  # Set start values for recursion on diagonal and residual index for correct sampling
  C_star <- diagonal
  res_index <- 1
  
  # Counter for resample attempts
  res_resample_future_values_bs <- integer(N - 1)
  
  # Loop for all development periods, starting with period 1
  for (l in seq_len(N-1)) {
    # Check which AYs are developed for this period -> we do it period by period!
    dev_AYs <- (N-l+1):N
    
    # Get correct number of residuals for this dev period
    res_for_current_period <- future_residuals[res_index:(res_index+l-1)]
    
    # Only develop the AY-values that need to be developed in this step!
    C_star_minus <- C_star[dev_AYs]
    
    C_star_new <- cl_factors[l] * C_star_minus + cl_sigmas[l] * sqrt(C_star_minus) * res_for_current_period
    
    # Make sure that no negative cumulative claims are generated, so resample then (very rare)
    negative <- C_star_new < 0
    
    attempt <- 0
    
    while (any(negative)) {
      
      attempt <- attempt + 1
    
      if(attempt > 50){
        stop(
          "Could not obtain positive cumulative Mack claims after 50 tries."
        )
      }
      
      # Add to the counter
      res_resample_future_values_bs[l] <- res_resample_future_values_bs[l] + sum(negative)
      
      #warning(sum(negative), "negative cumulative claims in dev period", l, "; residuals resampled there.")
        
      # Redraw exactly the residuals that are responsible for the negative cumulative claims (happens seldom, so just a quick fix)
      res_for_current_period[negative] <- sample(
        residuals,
        size = sum(negative),
        replace = TRUE
      )
      
      # Then recalculate
      C_star_new[negative] <- cl_factors[l] * C_star_minus[negative] + cl_sigmas[l] * sqrt(C_star_minus[negative]) * res_for_current_period[negative]
      
      negative <- C_star_new < 0
    }
    
    C_star[dev_AYs] <- C_star_new
    
    res_index <- res_index + l
  }
  # Return reserve
  list(
    reserve_obs = sum(C_star - diagonal),
    res_resample_future_values_bs = res_resample_future_values_bs
  )
}

### 2. ODP-Model (actually using NB as a simplification)
sim_triangle_ODP <- function(alphas, betas, dispersion){
  # Checks
  if (dispersion <= 1) {
    stop("Dispersion must be greater 1 for this simulation function!")
  }
  if (length(alphas) != length(betas)) {
    stop("Must have symmetric parameter number")
  }
  
  
  # Create matrix of means and size for NB
  mean_matrix <- outer(alphas, betas)
  size_matrix <- mean_matrix / (dispersion - 1)
  
  # Simulate outcomes from each parameter pair for full rectangle
  
  incremental <- matrix(
    rnbinom(
      length(mean_matrix),
      mu = as.vector(mean_matrix),
      size = as.vector(size_matrix)
    ),
    nrow = nrow(mean_matrix),
    ncol = ncol(mean_matrix)
  )
  
  cumulative <- t(apply(incremental, MARGIN = 1, FUN = cumsum))
  
  # Return
  cumulative
}




### 3. LN-Model

sim_and_refit_devfactors_LN <- function(mus, sigma_sqs){
  ## Checks
  if (length(mus) != length(sigma_sqs)) {
    stop("mus and sigma_sqs must same length.")
  }
  
  if (any(sigma_sqs < 0)) {
    stop("sigma_sqs must be non-negative.")
  }
  
  # Simulate dev_factors 
  
    # Save element count per column (because each development period in the triangle has different number of factor observations)
    N <- length(mus) + 1
    elements_per_period <- N - seq_len(N-1)
  
  dev_factors <- matrix(
    NA_real_,
    nrow = N,
    ncol = N-1
  )
  
  upper_id <- row(dev_factors) + col(dev_factors) <= N
  
  dev_factors[upper_id] <- rlnorm(
      sum(elements_per_period),
      meanlog = rep(mus, times = elements_per_period),
      sdlog = rep(sqrt(sigma_sqs), times = elements_per_period)
    )
  
  # Calculate parameter estimators
  mus_refit <- colMeans(log(dev_factors), na.rm = TRUE)
  sigma_sqs_refit <- apply(log(dev_factors), MARGIN = 2, FUN = var, na.rm = TRUE)
  
  # Set last variance parameter via Mack approximation
  sigma_sqs_refit[length(sigma_sqs_refit)] <- min(c(sigma_sqs_refit[length(sigma_sqs_refit)-1]^2/sigma_sqs_refit[length(sigma_sqs_refit)-2], sigma_sqs_refit[length(sigma_sqs_refit)-1], sigma_sqs_refit[length(sigma_sqs_refit)-2]))
  
  list(mus_refit = mus_refit, sigma_sqs_refit = sigma_sqs_refit)
}


sim_triangle_LN <- function(mus, sigma_sqs, start_values){
  ## Checks
  if (length(mus) != length(sigma_sqs)) {
    stop("mus and sigma_sqs must same length.")
  }
  
  if (any(sigma_sqs < 0)) {
    stop("sigma_sqs must be non-negative.")
  }
  
  if (any(start_values < 0)) {
    stop("start_values must be non-negative.")
  }
  if (length(start_values)-1 != length(mus)) {
    stop("Triangle must be symmetric!")
  }
  
  N <- length(start_values)
  
  # Simulate dev_factors 
  dev_factors <- matrix(
    rlnorm(
      N * (N-1),
      meanlog = rep(mus, each = N),
      sdlog = rep(sqrt(sigma_sqs), each = N)
    ),
    nrow = N,
    ncol = N-1
  )
  
  # Create empty Cij-Matrix
  cumulative <- matrix(
    NA_real_,
    nrow = N,
    ncol = N
  )
  
  # Set starting values
  cumulative[, 1] <- start_values
  
  # Calculate full claims data
  for (j in seq_len(N-1)) {
    cumulative[, j+1] = cumulative[,j] * dev_factors[,j]
  }
  
  # Return
  cumulative
}

# Simulating empirical (conditional) reserves from parameters and diagonal

sim_reserve_LN <- function(diagonal, mus, sigma_sqs){
  N <- length(diagonal)
  # Keep only resp. sums of parameters because of LN-distribution of latest-to-ultimate-factor
  parameters <- list(latest_to_ultimate_mus = c(0, cumsum(rev(mus))),
                           latest_to_ultimate_sigma_sqs = c(0, cumsum(rev(sigma_sqs))))
  
  # Simulate lower triangle from it
  latest_to_ultimate_factors <- rlnorm(
    N,
    meanlog = parameters$latest_to_ultimate_mus,
    sdlog = sqrt(parameters$latest_to_ultimate_sigma_sqs)
  )
  
  ultimate <- diagonal * latest_to_ultimate_factors
  
  # Return total reserve
  sum(ultimate - diagonal)
}

################# C. Bootstrap Routine ################

### 1. Mack-Bootstrap
bootstrap_mack <- function(triangle, B){
  # Checks
  
  if (nrow(triangle) != ncol(triangle)) {
    stop("Matrix must be symmetric!")
  } 
  
  if (any(triangle <= 0, na.rm = TRUE)) {
      warning("Bootstrap Mack failed, non-positive value in triangle.")
      result <-   list(
        mean = NA_real_,
        sd = NA_real_,
        q50 = NA_real_,
        q65 = NA_real_,
        q80 = NA_real_,
        q95 = NA_real_
      )
      result$res_resample_dev_factor <- rep(NA_real_, N - 1)
      result$res_resample_future_values <- rep(NA_real_, N - 1)
    return(result)
    }
  
  # Step 0: Save diagonal of triangle (latest values)
  N <- nrow(triangle)
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  
  # Step 1+2: Fitting for original estimators and residuals
  parameters_est <- fit_triangle_mack(triangle, emp_norm = TRUE)
  
  
  ### LOOP: Prepare empty reserve vector and constant quantities, then bootstrap for observations
  reserve_bs_obs <- numeric(B)
  
    # Before: Prepare relevant quantities for pseudo-data generation (repating them in a matrix to create the matrix of resampled dev factors)
      shifted_triangle <- triangle[,-N, drop = FALSE]
  
      mean_rep <- matrix(
      rep(parameters_est$cl_factors, each = N),
      nrow = N,
      ncol = N-1
    )
    
    factor_rep <- sweep(
      1 / sqrt(shifted_triangle),
      MARGIN = 2,
      STATS = parameters_est$cl_sigmas,
      FUN = "*"
    )
    
    cell_weights <- matrix(
      0,
      nrow = N,
      ncol = N-1
    )
    
    cell_weights[row(cell_weights) + col(cell_weights) <= N] <- shifted_triangle[row(shifted_triangle) + col(shifted_triangle) <= N]
    
    cell_weights <- sweep(
      cell_weights,
      MARGIN = 2,
      STATS = colSums(cell_weights),
      FUN = "/"
    )
    
    # Count resamples
    res_resample_dev_factor_total <- numeric(N-1)
    res_resample_future_values_total <- numeric(N-1)
    
  for (i in seq_len(B)){
    # Step 3.1 and 3.2 Simulate new dev factors and get refitted parameters 
    refit_cl_factors <- sim_and_refit_devfactors_mack(mean_rep, factor_rep, cell_weights, parameters_est$residuals)
    
    res_resample_dev_factor_total <- res_resample_dev_factor_total + refit_cl_factors$res_resample_dev_factor_bs
    
    # Step 3.3: Simulate lower triangle recursively from it and return reserves
    sim_reserve <- sim_reserve_mack(diagonal, refit_cl_factors$cl_factors, parameters_est$cl_sigmas, parameters_est$residuals)
    reserve_bs_obs[i] <- sim_reserve$reserve_obs
    
    res_resample_future_values_total <- res_resample_future_values_total + sim_reserve$res_resample_future_values_bs
    
  }
    # Step 4: Give results on reserve
    result <- reserve_results(reserve_bs_obs)
    result$res_resample_dev_factor <- res_resample_dev_factor_total
    result$res_resample_future_values <- res_resample_future_values_total
    result
}

### 2. ODP-Bootstrap
bootstrap_odp <- function(triangle, B){
  
  # Step 1-3 Do whole bootstrap with integrated packages function
  bootstrap_model <- BootChainLadder(
    triangle,
    R = B,
    process.distr = "od.pois",
  )
  
    # Get bootstrap reserve realizations
  reserve_bs_obs <- bootstrap_model$IBNR.Total
  
  # Step 4: Give results on reserve
  reserve_results(reserve_bs_obs)
}

### 3. LN-Bootstrap
bootstrap_LN <- function(triangle, B){
  # Checks
  
  if (nrow(triangle) != ncol(triangle)) {
    stop("Matrix must be symmetric!")
  } 
  if (any(triangle <= 0, na.rm = TRUE)) {
    stop("Values must be positive!")
  }
  
  # Step 0: Save diagonal of triangle (latest values)
  N <- nrow(triangle)
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  
  # Step 1: Fitting for original estimators
  parameters_est <- fit_triangle_LN(triangle)
  
  ### LOOP: Prepare empty reserve vector, then bootstrap for observations
  reserve_bs_obs <- numeric(B)
  
  for (i in seq_len(B)){
  # Step 3.1+3.2: Simulate new dev factors and get refitted parameters 
    parameters_refit <- sim_and_refit_devfactors_LN(parameters_est$mus, parameters_est$sigma_sqs)
    
  # Step 3.3 Simulate total reserve
    reserve_bs_obs[i] <- sim_reserve_LN(diagonal, parameters_refit$mus_refit, parameters_refit$sigma_sqs_refit)
  }
  
  # Step 4: Give results on reserve
  reserve_results(reserve_bs_obs)
}

############################ APPENDIX
### LN-Simulation only using the paramaters once estimated (so process error only)

simple_ln_simulation <- function(triangle, B){
  # Step 0: Save diagonal of triangle (latest values)
  N <- nrow(triangle)
  diagonal <- triangle[cbind(seq_len(N), N:1)]
  
  # Step 1: Fitting for original estimators
  parameters_est <- fit_triangle_LN(triangle)
  
  ### LOOP: Prepare empty reserve vector, then simulate new observations
  reserve_bs_obs <- numeric(B)
  
  for (i in seq_len(B)){
    # Keep only resp. sums of parameters because of LN-distribution of latest-to-ultimate-factor
    parameters_to_ult <- list(latest_to_ultimate_mus = c(0, cumsum(rev(parameters_est$mus))),
                             latest_to_ultimate_sigma_sqs = c(0, cumsum(rev(parameters_est$sigma_sqs))))
    
    # Step 3.3: Simulate lower triangle from it
    latest_to_ultimate_factors <- rlnorm(
      N,
      meanlog = parameters_to_ult$latest_to_ultimate_mus,
      sdlog = sqrt(parameters_to_ult$latest_to_ultimate_sigma_sqs)
    )
    
    ultimate <- diagonal * latest_to_ultimate_factors
    
    # Return total reserve
    reserve_bs_obs[i] <- sum(ultimate - diagonal)
  }
  
  # Step 4: Give results on reserve
  reserve_results(reserve_bs_obs)
}