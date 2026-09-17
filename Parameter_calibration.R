#### Getting parameter from calibration triangle

library(ChainLadder)
library(readxl)

calibrate_parameters <- function(file_path, triangle_name){
  
  ############### Read calibration triangle
  
  triangle_data <- read_excel(file_path)
  
  triangle <- as.matrix(triangle_data[,-1])
  
  ########### Retrieve parameters #########
  generating_models <- c("Mack", "ODP", "LN")
  
  source("Functions.R")
  
  ###### Now parameters calibration
  # Create array for parameters
  sim_parameters <- setNames(
    vector("list", length(generating_models)),
    generating_models)
  
  sim_parameters[["Mack"]] <- fit_triangle_mack(triangle, emp_norm = TRUE)
  sim_parameters[["ODP"]] <- fit_triangle_ODP(as.triangle(triangle))
  sim_parameters[["LN"]] <- fit_triangle_LN(triangle)
  
  ## Save parameters fitted based on this calibration triangle
  
  saveRDS(sim_parameters, paste0(triangle_name, "calibrated_parameters.rds"))
}
