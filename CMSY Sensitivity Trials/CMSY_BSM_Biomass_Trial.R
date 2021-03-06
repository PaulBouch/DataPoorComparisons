library (datalimited2)

#####  Cod data
### CMSY with BSM (includes survey index)

cod.data = read.csv("Data//Cod_7_ek//Cod.data.csv")

###### Create Results Table ###################
write.results = function(results, df){
  output = data.frame(Profile, resilience,
                      df[1,2], df[1,3], df[1,4],
                      df[2,2], df[2,3], df[2,4],
                      df[3,2], df[3,3], df[3,4],
                      df[4,2], df[4,3], df[4,4],
                      df[5,2], df[5,3], df[5,4])
  rbind (results, output)
}

results_biomass_BSM = data.frame(matrix(ncol = 17, nrow = 0))

###############################################
#### Test different combinations of starting and final prior ranges
#### Using Medium Resilience, which is assumed to be correct

resilience = "Medium"

### H to H Profile
Profile ="H_H"
output_HH <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                  btype="CPUE", resilience = resilience, 
                  stb.low = 0.5,
                  stb.hi = 0.9,
                  endb.low = 0.5,
                  endb.hi = 0.9)
plot_dlm(output_HH)
ref_HH <- output_HH[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_HH)

### H to M Profile
Profile ="H_M"
output_HM <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.5,
                 stb.hi = 0.9,
                 endb.low = 0.2,
                 endb.hi = 0.6)
plot_dlm(output_HM)
ref_HM <- output_HM[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_HM)

### H to L Profile
Profile ="H_L"
output_HL <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.5,
                 stb.hi = 0.9,
                 endb.low = 0.01,
                 endb.hi = 0.4)
plot_dlm(output_HL)
ref_HL <- output_HL[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_HL)

### M to H Profile
Profile ="M_H"
output_MH <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.2,
                 stb.hi = 0.6,
                 endb.low = 0.5,
                 endb.hi = 0.9)
plot_dlm(output_MH)
ref_MH <- output_MH[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_MH)

### M to M Profile
Profile ="M_M"
output_MM <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.2,
                 stb.hi = 0.6,
                 endb.low = 0.2,
                 endb.hi = 0.6)
plot_dlm(output_MM)
ref_MM <- output_MM[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_MM)

### M to L Profile
Profile ="M_L"
output_ML <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.2,
                 stb.hi = 0.6,
                 endb.low = 0.01,
                 endb.hi = 0.4)
plot_dlm(output_ML)
ref_ML <- output_ML[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_ML)

### L to H Profile
Profile ="L_H"
output_LH <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.01,
                 stb.hi = 0.4,
                 endb.low = 0.5,
                 endb.hi = 0.9)
plot_dlm(output_LH)
ref_LH <- output_LH[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_LH)

### L to M Profile
Profile ="L_M"
output_LM <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.01,
                 stb.hi = 0.4,
                 endb.low = 0.2,
                 endb.hi = 0.6)
plot_dlm(output_LM)
ref_LM <- output_LM[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_LM)

### L to L Profile
Profile ="L_L"
output_LL <- bsm(year=cod.data$year, catch=cod.data$catch, biomass=cod.data$index1,
                 btype="CPUE", resilience = resilience, 
                 stb.low = 0.01,
                 stb.hi = 0.4,
                 endb.low = 0.01,
                 endb.hi = 0.4)
plot_dlm(output_LL)
ref_LL <- output_LL[["ref_pts"]]
results_biomass_BSM = write.results(results_biomass_BSM, ref_LL)

results_biomass_BSM = setNames(results_biomass_BSM, c("Profile", "Resilience", "r", "r.low", "r.high",
                                                            "k", "k.low", "k.high", 
                                                            "msy", "msy.low", "msy.high",
                                                            "bmsy", "bmsy.low", "bmsy.high",
                                                            "fmsy", "fmsy.low", "fmsy.high"))

write.csv(results_biomass_BSM, "CMSY Sensitivity Trials\\Results\\results_biomass_BSM.csv")
