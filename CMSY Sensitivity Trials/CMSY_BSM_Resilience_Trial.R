library (datalimited2)

#####  Cod data
### CMSY with BSM (includes survey index)

###### Create Results Table ###################
write.results = function(results, df){
  output = data.frame(resilience,
                      df[1,2], df[1,3], df[1,4],
                      df[2,2], df[2,3], df[2,4],
                      df[3,2], df[3,3], df[3,4],
                      df[4,2], df[4,3], df[4,4],
                      df[5,2], df[5,3], df[5,4])
  rbind (results, output)
}

results_resilience_BSM = data.frame(matrix(ncol = 16, nrow = 0))

###############################################


### Very Low Resilience
resilience = "Very low"
output_vlow <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                     btype="biomass", resilience = resilience)
plot_dlm(output_vlow)
ref_vlow <- output_vlow[["ref_pts"]]
results_resilience_BSM = write.results(results_resilience_BSM, ref_vlow)

### Low Resilience
resilience = "Low"
output_low <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                  btype="biomass", resilience = resilience)
plot_dlm(output_low)
ref_low <- output_low[["ref_pts"]]
results_resilience_BSM = write.results(results_resilience_BSM, ref_low)

### Med Resilience
resilience = "Medium"
output_med <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                  btype="biomass", resilience = resilience)
plot_dlm(output_med)
ref_med <- output_med[["ref_pts"]]
results_resilience_BSM = write.results(results_resilience_BSM, ref_med)

### High Resilience
resilience = "High"
output_high <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                   btype="biomass", resilience = resilience)
plot_dlm(output_high)
ref_high <- output_high[["ref_pts"]]
results_resilience_BSM = write.results(results_resilience_BSM, ref_high)

results_resilience_BSM = setNames(results_resilience_BSM, c("Resilience", "r", "r.low", "r.high",
                                                    "k", "k.low", "k.high", 
                                                    "msy", "msy.low", "msy.high",
                                                    "bmsy", "bmsy.low", "bmsy.high",
                                                    "fmsy", "fmsy.low", "fmsy.high"))

write.csv(results_resilience_BSM, "CMSY Sensitivity Trials\\Results\\results_resilience_BSM.csv")
