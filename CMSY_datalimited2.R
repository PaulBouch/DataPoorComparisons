library (datalimited2)

#####  Cod data
### No abundance indices. Just catch data

output <- cmsy2(year=cod.data$year, catch=cod.data$catch, resilience = "Low")
plot_dlm(output)

ref_pts <- output[["ref_pts"]]
ref_ts <- output[["ref_ts"]]


#### CMSY with BSM

output <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
              btype="biomass", resilience = "Low")
plot_dlm(output)

# Extract reference points and time series from output
ref_pts <- output[["ref_pts"]]
ref_ts <- output[["ref_ts"]]