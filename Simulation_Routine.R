# Main simulation routine. This is meant to be run line by line, so you can actually inspect the process and the results.
# Section 0 has a lot of manual input, so please put in your specific paths, triangle_names etc. 

library(ChainLadder)
library(cli)

############### 0. Read parameter set and prepare simulation ##########
name <- "[NAME OF PARAMETER SET]"
parameter_set_path <- "[PATH OF PARAMETER SET]" #this expects the parameters in the format of the output of "Parameter_calibration.R"

sim_parameters <- readRDS(parameter_set_path)

set.seed(1512)

# Define the number of triangles for each generating model
nb_triangle_sim <- 10000
# Define how many bootstrap runs for the fitting section
B <- 1000

################ 0.1 Simulation preparation with no manual input required ######
source("Functions.R")

generating_models <- c("Mack", "ODP", "LN")


############# 1 . Simulate triangles from given parameter set for each generating model (ONLY ONCE PER PARAMETER_SET) ##########
# Select triangle (symmetric) size
N <- length(sim_parameters$ODP$alpha_hat)

# Create array for the simulated triangles
rectangles <- array(
  NA_real_,
  dim = c(
    N,
    N,
    nb_triangle_sim,
    length(generating_models)
  ),
  dimnames = list(
    AY = NULL,
    dev = NULL,
    triangle_nb = NULL,
    generating_model = generating_models
  )
)

# Count the errors that happen in Mack simulation because of negative cumulative claims
mack_simulation_errors <- integer(N-1)

# Simulate for each generating model, using the parameters above
for (gen_m in generating_models){
  for (l in seq_len(nb_triangle_sim)){
    if (gen_m == "Mack"){
      Mack_sims <- sim_triangle_Mack(sim_parameters[[gen_m]]$cl_factors, sim_parameters[[gen_m]]$cl_sigmas, sim_parameters[[gen_m]]$residuals, sim_parameters[[gen_m]]$start_values)
      rectangles[,,l,gen_m] <- Mack_sims$rectangle
      mack_simulation_errors <- mack_simulation_errors + Mack_sims$resample_count_mack_simulation
      next
    }
    
    rectangles[,,l,gen_m] <- sim_triangle(gen_m, sim_parameters[[gen_m]])
  }
}

saveRDS(rectangles, file = paste0("rectangles_",name,"_simulations.rds")) ## WATCH OUT; THESE ARE RECTANGLES!

# Hide lower part of rectanlge, so only triangles
claims <- rectangles

for (gen_m in generating_models){
  for (l in seq_len(nb_triangle_sim)){
    claims[,,l,gen_m] <- hide_future(claims[,,l,gen_m])
  }
}

saveRDS(claims, file = paste0("claims_",name,"_simulations.rds")) ## WATCH OUT; THESE ARE TRIANGLES NOW!

############## 1.1 Calculate 'true' ultimate, MSEP and Quantile ###########

# Access the triangles, if not loaded in already
#claims <- readRDS("X.rds")

# Create array for the data that needs to be saved
fitting_models <- c("true", "Mack", "ODP", "LN")
statistics <- c("mean_bs", "mean_ana", "msep_bs", "msep_ana", "q50", "q65", "q80", "q95")

full_results <- array(
  NA_real_,
  dim = c(
    length(generating_models),
    nb_triangle_sim,
    length(fitting_models),
    length(statistics)
  ),
  dimnames = list(
    generating_model = generating_models,
    triangle_nb = NULL,
    fitted_model = fitting_models,
    statistic = statistics
  )
)

# Create array for error analysis (debugging when triangles are not compatible)

error_types <- c("ODP_full_column_neg", "ODP_negative_increment", "Mack_dev_factor_resample", "Mack_future_resample")

error_analysis <- array(
  0,
  dim = c(
    length(generating_models),
    nb_triangle_sim,
    length(fitting_models),
    length(error_types),
    N
  ),
  dimnames = list(
    gen_m = generating_models,
    triangle = NULL,
    fit_m = fitting_models,
    error_type = error_types,
    dev = NULL
  )
)


### Actually calculate 'true' values of ultimate, MSEP, quantiles etc. for the triangles with the known parameters

# Create progress bar
total_ana_fits <- (length(generating_models) -1) * nb_triangle_sim

progress_bar_1 <- cli_progress_bar(
  name = "Analytical fitting of simulated triangles",
  total = total_ana_fits, 
  type = "iterator",
  auto_terminate = FALSE,
  clear = FALSE,
  format = paste0(
    "{pb_bar} {pb_percent} / ",
    "{pb_current} / {pb_total} completed / ",
    "{pb_total - pb_current} remaining / ",
    "{pb_eta_str} / {pb_status}"
  )
)

# Mack and LN have for each triangle conditional 'true' results, so they have to be calculated for each triangle
for (gen_m in c("Mack", "LN")){
  for (l in seq_len(nb_triangle_sim)){
    # Show current position
    cli_progress_update(
      inc = 0,
      id = progress_bar_1,
      status = sprintf(
        "gen_m: %s / triangle: %d/%d",
        gen_m,
        l,
        nb_triangle_sim
      )
    )
    
    
    # Calculate full result
    result <- known_parameters_result(claims[,,l,gen_m], sim_parameters[[gen_m]], gen_m, nb_triangle_sim)
    
    # Map to right matrix elements
    full_results[gen_m,l,"true","mean_bs"] <- result$MC_mean
    full_results[gen_m,l,"true","mean_ana"] <- result$mean
    full_results[gen_m,l,"true","msep_bs"] <- result$MC_sd #which does not mean msep here, but rather sd of the real reserve distribution
    full_results[gen_m,l,"true","msep_ana"] <- result$sd #which does not mean msep here, but rather sd of the real reserve distribution
    
    full_results[gen_m,l,"true","q50"] <- result$q50
    full_results[gen_m,l,"true","q65"] <- result$q65
    full_results[gen_m,l,"true","q80"] <- result$q80
    full_results[gen_m,l,"true","q95"] <- result$q95
    
    if(gen_m == "Mack") {
      error_analysis[gen_m, l, "true", "Mack_future_resample", seq_len(N-1)] <- result$res_resample_future_values
    }
    
    # Update index of progress bar
    cli_progress_update(
      id = progress_bar_1,
      inc = 1)
  }
}
# ODP-'true'-results only depend on the parameters, so only needed to be calculated once
odp_result <- results_odp_known_parameters(sim_parameters[["ODP"]]$alpha_hat, sim_parameters[["ODP"]]$beta_hat, sim_parameters[["ODP"]]$dispersion, nb_triangle_sim)
  
for (l in seq_len(nb_triangle_sim)){
  full_results["ODP",l,"true","mean_bs"] <- odp_result$MC_mean
  full_results["ODP",l,"true","mean_ana"] <- odp_result$mean
  full_results["ODP",l,"true","msep_bs"] <- odp_result$MC_sd #which does not mean msep here, but rather sd of the real reserve distribution
  full_results["ODP",l,"true","msep_ana"] <- odp_result$sd #which does not mean msep here, but rather sd of the real reserve distribution
  
  full_results["ODP",l,"true","q50"] <- odp_result$q50
  full_results["ODP",l,"true","q65"] <- odp_result$q65
  full_results["ODP",l,"true","q80"] <- odp_result$q80
  full_results["ODP",l,"true","q95"] <- odp_result$q95
}

############## 2. Fitting the simulated triangles in different model settings ##############

# Create progress bar
total_tasks <- length(generating_models) * nb_triangle_sim * (length(fitting_models)-1)

progress_bar2 <- cli_progress_bar(
  name = "Fitting of simulated triangles",
  total = total_tasks, 
  type = "iterator",
  auto_terminate = FALSE,
  clear = FALSE,
  format = paste0(
    "{pb_bar} {pb_percent} / ",
    "{pb_current} / {pb_total} completed / ",
    "{pb_total - pb_current} remaining / ",
    "{pb_eta_str} / {pb_status}"
  )
)

# Loop over each generating model, for each triangle generated by that and each fitting model
for (gen_m in generating_models){
  for (l in seq_len(nb_triangle_sim)){
    triangle <- claims[,,l,gen_m]
    for (fit_m in fitting_models[-1]){
      # Show current position
      cli_progress_update(
        inc = 0,
        id = progress_bar2,
        status = sprintf(
          "gen_m: %s / triangle: %d/%d / fit_m: %s",
          gen_m,
          l,
          nb_triangle_sim,
          fit_m
        )
      )
      
      #Skip combinations already completed
      # if (!is.na(full_results[gen_m, l, fit_m, "q95"])) {
      #   cli_progress_update(
      #     id = progress_bar2,
      #     inc = 1
      #   )
      #   next
      # }
      
      # Skip if any value is non-positive
      if (any(triangle <= 0, na.rm = TRUE)) {
        warning("Triangle with non-positive cum. value encoutered. Skipped.")
        cli_progress_update(
          id = progress_bar2,
          inc = 1)
        next
      }
      
     # # For debugging: Print current triangle
    #  message(
    #    sprintf("gen_m: %s | triangle: %d", gen_m, l)
    #  )
      
      
      # Fit the model analytical
      ana_result <- analytical_result(triangle, fit_m)
      
      # Save analytical irregularities
      
      if (fit_m == "ODP") {
        error_analysis[gen_m, l, fit_m, "ODP_full_column_neg", ana_result$failed] <- 1
        
        # Skip triangle for ODP fit on negative incremental columns, no reasonable values from BS as well then
        if (length(ana_result$failed)>0) {
          cli_progress_update(
            id = progress_bar2,
            inc = 1)
          next
        }
        
        error_analysis[gen_m, l, fit_m, "ODP_negative_increment", ana_result$fallback_columns] <- 1
      }
      
      
      # Fit the model bootstrap
      bs_result <- bootstrap_result(triangle, fit_m, B)
      
      

      # Save BS-irregularities      
      if(fit_m == "Mack") {
        error_analysis[gen_m, l, fit_m, "Mack_dev_factor_resample", seq_len(N-1)] <- bs_result$res_resample_dev_factor
        error_analysis[gen_m, l, fit_m, "Mack_future_resample", seq_len(N-1)] <- bs_result$res_resample_future_values
      }
      
     
      
      
      # Save results
      full_results[gen_m, l, fit_m, "mean_ana"] <- ana_result$mean
      full_results[gen_m, l, fit_m, "msep_ana"] <- ana_result$msep
      full_results[gen_m, l, fit_m, "mean_bs"] <- bs_result$mean
      full_results[gen_m, l, fit_m, "msep_bs"] <- bs_result$sd
      full_results[gen_m, l, fit_m, "q50"] <- bs_result$q50
      full_results[gen_m, l, fit_m, "q65"] <- bs_result$q65
      full_results[gen_m, l, fit_m, "q80"] <- bs_result$q80
      full_results[gen_m, l, fit_m, "q95"] <- bs_result$q95
    
      
      
      # Update index of progress bar
      cli_progress_update(
        id = progress_bar2,
        inc = 1
      )
    }
  }
}


# Save results
saveRDS(full_results, paste0("full_results_",name,".rds"))


################## 3. Some error checks #########################

  # 3.0.1 Check if there were many errors in Mack simulation
  mack_simulation_errors
  
  # 3.0.2 Check if while calculating 'true' Mack result, any resampling was neccessary
  apply(error_analysis["Mack",,"true","Mack_future_resample",], 
        c(2),
        sum)
  
  
  # 3.0.3 How many triangles had negative incremental values?
  
  triangle_count_neg_incremental_values <- colSums(apply(
    claims,
    c(3,4),
    function(x) any(cum2incr(x) <= 0, na.rm = TRUE)
  ))
  
  triangle_count_neg_incremental_values
  
  
  ## 3.1 How many triangles were skipped due to ODP_negative columns?
  ODP_skipped <- apply(error_analysis[, ,"ODP","ODP_full_column_neg",], c(1,3), sum, na.rm = TRUE)
  ODP_skipped
  
  ## 3.2 How many triangles should have been skipped for all models because of negative cumulative values?
  triangle_skipped <- colSums(apply(
      claims,
      c(3,4),
      function(x) any(x <= 0, na.rm = TRUE)
    ))
  
  triangle_skipped
  
  # 3.3 Check if only these triangles were actually skipped
  skipped <- apply(
    full_results[,,,"mean_bs"],
    c(1,3),
    function(x) sum((is.na(x)))
  )
  
  skipped
  
  # 3.4.1 Mack check for number of dev_factor_resamples
  apply(error_analysis[,,"Mack","Mack_dev_factor_resample",], 
        c(1,3),
        sum)
  
  # 3.4.2 Mack check for number of future cumulative claims generation resample
  apply(error_analysis[,,"Mack","Mack_future_resample",], 
        c(1,3),
        sum)
  
  # 3.5 Were Analytic and BS result close?
  
  # Mean
  apply(full_results[,,,"mean_bs"] / full_results[,,,"mean_ana"] -1,
        3,
        mean,
        na.rm = TRUE
        )
  
  apply(full_results[,,,"mean_bs"] / full_results[,,,"mean_ana"] -1,
        3,
        median,
        na.rm = TRUE
  )
  
    # more granular
    apply(full_results[,,,"mean_bs"] / full_results[,,,"mean_ana"] -1,
          c(1,3),
          mean,
          na.rm = TRUE
    )
  
    # and its sd (sd of the mean-diff):
    apply(full_results[,,,"mean_bs"] / full_results[,,,"mean_ana"] -1,
          3,
          sd,
          na.rm = TRUE
    )
  
  
  # MSEP
  apply(full_results[,,,"msep_bs"] / full_results[,,,"msep_ana"] -1,
        3,
        mean,
        na.rm = TRUE
        )
  
  apply(full_results[,,,"msep_bs"] / full_results[,,,"msep_ana"] -1,
        3,
        median,
        na.rm = TRUE
  )
      
    # more granular
      apply(full_results[,,,"msep_bs"] / full_results[,,,"msep_ana"] -1,
            c(1,3),
            mean,
            na.rm = TRUE
      )
  
    # and its sd (sd of the MSEP-diff)
    apply(full_results[,,,"msep_bs"] / full_results[,,,"msep_ana"] -1,
          3,
          sd,
          na.rm = TRUE
    )


#################### 4. Comparison of model specifications ######################


# Define relevant quantities for comparison
comparison_statistics <- c("mean_diff_bs", 
                           "mean_diff_rel_bs", 
                           "mean_diff_ana", 
                           "mean_diff_rel_ana", 
                           "msep_diff_bs", 
                           "msep_diff_rel_bs",
                           "msep_diff_ana", 
                           "msep_diff_rel_ana",
                           "q50_diff",
                           "q50_diff_rel",
                           "q65_diff",
                           "q65_diff_rel",
                           "q80_diff",
                           "q80_diff_rel",
                           "q95_diff",
                           "q95_diff_rel",
                           "q50_diff_rel_bs_norm",
                           "q65_diff_rel_bs_norm",
                           "q80_diff_rel_bs_norm",
                           "q95_diff_rel_bs_norm",
                           "q50_diff_rel_ana_norm",
                           "q65_diff_rel_ana_norm",
                           "q80_diff_rel_ana_norm",
                           "q95_diff_rel_ana_norm")

# Define comparison types: observed diffs or regular diffs
comparison_types <- c("absolute", "regular")

# Define moments to be compared: mean and sd (for constructing MC_confidence)
moments <- c("mean", "sd")

# Create comparsion array to be filled by loop afterwards
comparison_array <- array(
  NA_real_,
  dim = c(
    length(generating_models),
    length(fitting_models),
    length(fitting_models),
    length(comparison_statistics),
    length(moments),
    length(comparison_types)
  ),
  dimnames = list(
    generating_model = generating_models,
    fitted_model_1 = fitting_models,
    fitted_model_2 = fitting_models,
    comparison_statistic = comparison_statistics,
    moments = moments,
    comparison_type = comparison_types
  )
)

# Fill comparison matrix

# Fix generating model
for (gen_m in generating_models){
  # Loop over all comparison-pairs
  for (fit_m1 in fitting_models){
    for(fit_m2 in fitting_models){
      for(moment in moments){
        for (comp in comparison_types){
          # Calculate empirical mean diffs
          comparison_array[gen_m, fit_m1, fit_m2, "mean_diff_bs", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "mean_bs"]-full_results[gen_m, , fit_m2, "mean_bs"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "mean_diff_rel_bs", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "mean_bs"]/full_results[gen_m, , fit_m2, "mean_bs"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "mean_diff_ana", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "mean_ana"]-full_results[gen_m, , fit_m2, "mean_ana"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "mean_diff_rel_ana", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "mean_ana"]/full_results[gen_m, , fit_m2, "mean_ana"] - 1 , comp, moment)
          
          # Calculate empirical MSEP diffs
          comparison_array[gen_m, fit_m1, fit_m2, "msep_diff_bs", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "msep_bs"]-full_results[gen_m, , fit_m2, "msep_bs"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "msep_diff_rel_bs", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "msep_bs"]/full_results[gen_m, , fit_m2, "msep_bs"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "msep_diff_ana", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "msep_ana"]-full_results[gen_m, , fit_m2, "msep_ana"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "msep_diff_rel_ana", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "msep_ana"]/full_results[gen_m, , fit_m2, "msep_ana"] - 1 , comp, moment)
          
          # Calculate empirical Quantile diffs
          comparison_array[gen_m, fit_m1, fit_m2, "q50_diff", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q50"]-full_results[gen_m, , fit_m2, "q50"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q50_diff_rel", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q50"]/full_results[gen_m, , fit_m2, "q50"] - 1 , comp, moment)
          
          comparison_array[gen_m, fit_m1, fit_m2, "q65_diff", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q65"]-full_results[gen_m, , fit_m2, "q65"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q65_diff_rel", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q65"]/full_results[gen_m, , fit_m2, "q65"] - 1 , comp, moment)
          
          comparison_array[gen_m, fit_m1, fit_m2, "q80_diff", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q80"]-full_results[gen_m, , fit_m2, "q80"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q80_diff_rel", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q80"]/full_results[gen_m, , fit_m2, "q80"] - 1 , comp, moment)
          
          comparison_array[gen_m, fit_m1, fit_m2, "q95_diff", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q95"]-full_results[gen_m, , fit_m2, "q95"], comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q95_diff_rel", moment, comp] <- toggle_comp(full_results[gen_m, , fit_m1, "q95"]/full_results[gen_m, , fit_m2, "q95"] - 1 , comp, moment)
          
          # Calculate normalized empiricial Quantile diffs
          comparison_array[gen_m, fit_m1, fit_m2, "q50_diff_rel_bs_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q50"]-full_results[gen_m, , fit_m1, "mean_bs"]+full_results[gen_m, , fit_m2, "mean_bs"])/full_results[gen_m, , fit_m2, "q50"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q65_diff_rel_bs_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q65"]-full_results[gen_m, , fit_m1, "mean_bs"]+full_results[gen_m, , fit_m2, "mean_bs"])/full_results[gen_m, , fit_m2, "q65"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q80_diff_rel_bs_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q80"]-full_results[gen_m, , fit_m1, "mean_bs"]+full_results[gen_m, , fit_m2, "mean_bs"])/full_results[gen_m, , fit_m2, "q80"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q95_diff_rel_bs_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q95"]-full_results[gen_m, , fit_m1, "mean_bs"]+full_results[gen_m, , fit_m2, "mean_bs"])/full_results[gen_m, , fit_m2, "q95"] - 1 , comp, moment)
          
          # Analytically normalized empirical Quantile diffs as well
          comparison_array[gen_m, fit_m1, fit_m2, "q50_diff_rel_ana_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q50"]-full_results[gen_m, , fit_m1, "mean_ana"]+full_results[gen_m, , fit_m2, "mean_ana"])/full_results[gen_m, , fit_m2, "q50"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q65_diff_rel_ana_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q65"]-full_results[gen_m, , fit_m1, "mean_ana"]+full_results[gen_m, , fit_m2, "mean_ana"])/full_results[gen_m, , fit_m2, "q65"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q80_diff_rel_ana_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q80"]-full_results[gen_m, , fit_m1, "mean_ana"]+full_results[gen_m, , fit_m2, "mean_ana"])/full_results[gen_m, , fit_m2, "q80"] - 1 , comp, moment)
          comparison_array[gen_m, fit_m1, fit_m2, "q95_diff_rel_ana_norm", moment, comp] <- toggle_comp((full_results[gen_m, , fit_m1, "q95"]-full_results[gen_m, , fit_m1, "mean_ana"]+full_results[gen_m, , fit_m2, "mean_ana"])/full_results[gen_m, , fit_m2, "q95"] - 1 , comp, moment)
        }
      }
    }
  }
}
# Save results
# NB: diff means diff of model 1 to model 2!!! (how much larger is fit_model1 quantity compared to fit_model2 quantity)

saveRDS(comparison_array, file = paste0("comparison_",name,".rds"))

# Save Error overview
saveRDS(error_analysis, file = paste0("erroranalysis_",name,".rds"))

# Export the comparison results to CSV as well
source("Export_Statistics.R")
export_statistics_to_csv(comparison_array, name)
