library (datalimited2)

#####  Cod data
### No abundance indices. Just catch data

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

results_resilience = data.frame(matrix(ncol = 16, nrow = 0))

###############################################


### Very Low Resilience
resilience = "Very low"
output_vlow <- cmsy2(year=cod.data$year, catch=cod.data$catch, resilience = resilience)
plot_dlm(output_vlow)
ref_vlow <- output_vlow[["ref_pts"]]
results_resilience = write.results(results_resilience, ref_vlow)

### Low Resilience
resilience = "Low"
output_low <- cmsy2(year=cod.data$year, catch=cod.data$catch, resilience = "Low")
plot_dlm(output_low)
ref_low <- output_low[["ref_pts"]]
results_resilience = write.results(results_resilience, ref_low)

### Med Resilience
resilience = "Medium"
output_med <- cmsy2(year=cod.data$year, catch=cod.data$catch, resilience = "Medium")
plot_dlm(output_med)
ref_med <- output_med[["ref_pts"]]
results_resilience = write.results(results_resilience, ref_med)

### High Resilience
resilience = "High"
output_high <- cmsy2(year=cod.data$year, catch=cod.data$catch, resilience = "High")
plot_dlm(output_high)
ref_high <- output_high[["ref_pts"]]
results_resilience = write.results(results_resilience, ref_high)

results_resilience = setNames(results_resilience, c("Resilience", "r", "r.low", "r.high",
                                               "k", "k.low", "k.high", 
                                               "msy", "msy.low", "msy.high",
                                               "fmsy", "fmsy.low", "fmsy.high",
                                               "bmsy", "bmsy.low", "bmsy.high"))

write.csv(results_resilience, "CMSY Sensitivity Trials\\Results\\results_resilience_CMSY.csv")
