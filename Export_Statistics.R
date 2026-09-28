### Export relevant statistics to csv (this expects the format of the output 'comparison_array' of "Simulation_Routine.R")

export_statistics_to_csv <- function(comparison_data, name){

  ## Some info about the structure
  generating_models <- c("Mack", "ODP", "LN")
  comparison_statistics <- c("mean_diff_bs", 
                             "mean_diff_rel_bs", 
                             "mean_diff_ana", 
                             "mean_diff_rel_ana", 
                             "msep_diff_bs", 
                             "msep_diff_rel_bs",
                             "msep_diff_ana", 
                             "msep_diff_rel_ana",
                             "process_err_diff_rel_ana",
                             "parameter_err_diff_rel_ana",
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
  moments <- c("mean", "sd")
  
  ## Actual export array
  comparison_export <- do.call(
    rbind,
    lapply(comparison_statistics, function(stat) {
      
      do.call(
        rbind,
        lapply(generating_models, function(gen_m) {
          
          do.call(
            rbind,
            lapply(moments, function(moment) {
              values <- comparison_data[
                gen_m,
                ,
                gen_m,
                stat,
                moment,
                "regular"
              ]
              
              data.frame(
                statistic = stat,
                moment = moment,
                gen_m = gen_m,
                t(values),
                check.names = FALSE,
                row.names = NULL
              )
            })
          )
        })
      )
    })
  )
  
  write.csv(
    comparison_export,
    paste0("comparison_results_",name,".csv"),
    row.names = FALSE
  )
}