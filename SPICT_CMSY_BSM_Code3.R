library (datalimited2)
library(spict)
options(scipen = 999)


data = read.csv("Data//Whg_7_a//Whg.a.Data.csv")
data$catch = data$catch / 1000

Stock = "Whg_7_a"

No.Index = ncol(data) - 3

###############################################################################
# Prior information for resilience (r) of stock ---------------------------

# Can either classify as High, Medium, Low or Very low or manually set the boundaries
pb.resilience = "Medium"
## Lower r estimate
pb.r.low = NA
## Upper r estimate
pb.r.hi = NA

if (is.na(pb.r.low))  {if (pb.resilience == "High"){pb.r.low = 0.6; pb.r.hi = 1.5}
  else if (pb.resilience == "Medium"){pb.r.low = 0.2; pb.r.hi = 0.8}
  else if (pb.resilience == "Low"){pb.r.low = 0.05; pb.r.hi = 0.5}
  else if (pb.resilience == "Very low"){pb.r.low = 0.015; pb.r.hi = 0.1}}


# Function to calculate moving average
ma <- function(x){
  x.1 <- stats::filter(x,rep(1/3,3),sides=1)
  x.1[1] <- x[1]
  x.1[2] <- (x[1]+x[2])/2
  return(x.1)
}


# Priors ------------------------------------------------------------------

###################################################
### Prior for r
###################################################
pb.start.r = c(pb.r.low, pb.r.hi)

### Mean of r range specified by the resilience
pb.log.r = log(mean(pb.start.r))
pb.log.r.sd = abs ((pb.log.r - log(pb.r.low))/2)

nyr = length(data$year)
ct <- ma(data$catch)

######################################################
### k Prior and biomass priors for CMSY
######################################################

pb.endbio_prior <- function(nyr, ct.raw, ct){
  rawct.ratio=ct.raw[nyr]/max(ct)
    endbio  <- if(ct[nyr]/max(ct) > 0.8) {c(0.4,0.8)} else if(rawct.ratio < 0.5) {c(0.01,0.4)} else {c(0.2,0.6)}
    
    # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
    if(endbio[2]==0.4){
      if(rawct.ratio< 0.05) {endbio[2] <- 0.1} else
        if(rawct.ratio< 0.15) {endbio[2] <- 0.2} else
          if(rawct.ratio< 0.35) {endbio[2] <- 0.3} else {endbio[2] <- 0.4}
    }
  return(endbio)
}

pb.endbio = pb.endbio_prior(nyr, data$catch, ct)

pb.k_prior <- function(endbio, r.low, r.hi, ct){
  # initial prior range of k values, assuming min k will be larger than max catch / prior for r
  if(mean(endbio) <= 0.5){
    start.k <- c(max(ct)/r.hi,4*max(ct)/r.low)
  }else{
    start.k <- c(2*max(ct)/r.hi, 12*max(ct)/r.low)
  }
  return(start.k)
}

pb.k.prior = pb.k_prior(pb.endbio, pb.r.low, pb.r.hi, ct)

pb.k.low = pb.k.prior[1]
pb.k.hi =  pb.k.prior[2] 

### Mean of k.low and k.high when defined by default for CMSY
pb.log.k =  log(mean(pb.k.prior))
pb.log.k.sd = (pb.log.k - log(pb.k.low))/4

###############################################################
###  prior for q
###############################################################


lyr               <- ifelse(mean(pb.start.r)>=0.5,5,10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish
mean.last.ct      <-mean(data$catch[(nyr-lyr):nyr],na.rm=T) # get mean catch of last years
mean.last.cpue    <-mean(data$index1[(nyr-lyr):nyr],na.rm=T) # get mean of CPUE of last years

gm.start.r      <- exp(mean(log(pb.start.r))) # get geometric mean of prior r range
if(mean(pb.endbio) >= 0.5) {  # if biomass is high
  q.1           <- mean.last.cpue*0.25*gm.start.r/mean.last.ct
  q.2           <- mean.last.cpue*0.5*pb.start.r[2]/mean.last.ct
} else {
  q.1           <- mean.last.cpue*0.5*gm.start.r/mean.last.ct
  q.2           <- mean.last.cpue*pb.start.r[2]/mean.last.ct
  }

pb.q.prior         <- c(q.1,q.2)
pb.q.low = pb.q.prior[1]
pb.q.hi = pb.q.prior[2]
pb.log.q = log(mean(pb.q.prior))
pb.log.q.sd = (pb.log.q -log(pb.q.low))/4

###########################################################################################
###########################################################################################
###########################################################################################
###                         Spict Function                                              ###
###########################################################################################
###########################################################################################
###########################################################################################

pb.spict = function(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                    pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                    pb.log.q, pb.log.q.sd, ...){
  inp <- list()
  inp$obsC <- data$catch
  inp$timeC <- data$year

  ### Use one abundance index or all
  if (Use.Index == 1){
  ### Use only one abundance index
  index = data[complete.cases(data$index1), ]
  inp$obsI <- index$index1
  inp$timeI <- index$year
  } else {
  if (No.Index == 2){
    ### Use 2 abundance index
  index1 = data[complete.cases(data$index1), ]
  index2 = data[complete.cases(data$index2), ]
  inp$obsI[[1]] <- index1$index1
  inp$obsI[[2]] <- index2$index2
  inp$timeI <- list(index1$year, index2$year)
  }
  if (No.Index == 3){
    ### Use 3 abundance index
    index1 = data[complete.cases(data$index1), ]
    index2 = data[complete.cases(data$index2), ]
    index3 = data[complete.cases(data$index3), ]
    inp$obsI[[1]] <- index1$index1
    inp$obsI[[2]] <- index2$index2
    inp$obsI[[3]] <- index3$index3
    inp$timeI <- list(index1$year, index2$year, index3$year)
  } else {
    ### Use 4 abundance index
    index1 = data[complete.cases(data$index1), ]
    index2 = data[complete.cases(data$index2), ]
    index3 = data[complete.cases(data$index3), ]
    index4 = data[complete.cases(data$index4), ]
    inp$obsI[[1]] <- index1$index1
    inp$obsI[[2]] <- index2$index2
    inp$obsI[[3]] <- index3$index3
    inp$obsI[[3]] <- index3$index3
    inp$timeI <- list(index1$year, index2$year, index3$year, index4$year) 
  }}

  ### fix n=2
  if (Use.n.prior == "Y"){inp$priors$logn <- c(log(2), 1e-3)}
  ### provide r prior
  if (Use.r.prior == "Y"){inp$priors$logr <- c(pb.log.r, pb.log.r.sd, 1)}
  ### Does not seem to like having k prior
  if (Use.k.prior == "Y"){inp$priors$logK <- c(pb.log.k, pb.log.k.sd, 1)}
  ### provide prior for q
  if (Use.q.prior == "Y"){inp$priors$logq <- c(pb.log.q, pb.log.q.sd, 1)}
  
  assign("last.warning", NULL, envir = baseenv())
  inp <- check.inp(inp)
  
  Prior_warning = warnings()
  
  
  
  #plots the Catch and index 
  plotspict.data(inp)
  
  assign("last.warning", NULL, envir = baseenv())
  # Run the model
  res <- fit.spict(inp)
  
  Fit_warning = warnings()
  
  # What are the calculated values
  summary(res)
  
  Convergence = capture.output(summary(res))[1]
  
  # plot everything
  # plot(res)
  # # Calculate the residuals
  # res <- calc.osa.resid(res)
  # # are the residuals ok?
  # plotspict.diagnostic(res)
  # sumspict.diagnostics(res)
  
  ####   Produce the reference value table
  tab1 <- sumspict.srefpoints(res);
  tab1 = as.data.frame(tab1)
  tab2 <- sumspict.states(res);
  tab2 = as.data.frame(tab2)
  spict.r = get.par('r', res)
  spict.k = get.par('K', res)
  spictrk = rbind(spict.r, spict.k)
  
  results_spict = setNames(data.frame(matrix(ncol = 48, nrow = 0)), c("Stock", "Method", "Index","FixN2",
                                                                       "Use.r", "Use.k", "Use.q",
                                                                       "Resilience", "Prior.r.low", "Prior.r.hi", 
                                                                       "Prior.k.low", "Prior.k.hi", 
                                                                       "Prior.log.k", "Prior.log.k.sd", 
                                                                       "Prior.q.low", "Prior.q.hi", "Prior.log.q", "Prior.log.q.sd", 
                                                                       "r", "r.low", "r.high",
                                                                       "k", "k.low", "k.high", 
                                                                       "msy", "msy.low", "msy.high",
                                                                       "fmsy", "fmsy.low", "fmsy.high",
                                                                       "bmsy", "bmsy.low", "bmsy.high",
                                                                       "b.end", "b.end.low", "b.end.hi",
                                                                       "bbmsy", "bbmsy.low", "bbmsy.hi",
                                                                       "f.end", "f.end.low", "f.end.hi",
                                                                       "ffmsy", "ffmsy.low", "ffmsy.hi",
                                                                       "Converge", "Prior.Warn", "Fit.Warn"))
  
  results_spict [1,1:18] = c(Stock, "SPICT", Use.Index, Use.n.prior, Use.r.prior,
                              Use.k.prior, Use.q.prior,
                              pb.resilience, pb.r.low, pb.r.hi, 
                              pb.k.low, pb.k.hi, pb.log.k, pb.log.k.sd,
                              pb.q.low, pb.q.hi, pb.log.q, pb.log.q.sd)
  
  results_spict [1,46:48] = c(Convergence[1],                               
                              if (length(Prior_warning)== 0){""} else  {names(Prior_warning)[1]},
                              if (length(Fit_warning)== 0){""} else  {names(Fit_warning)[1]}) 
  
  results_spict [1,19:45] = c(spictrk[1,2], spictrk[1,1], spictrk[1,3],
                               spictrk[2,2], spictrk[2,1], spictrk[2,3],
                               tab1[3,1], tab1[3,2], tab1[3,3], #msy
                               tab1[2,1], tab1[2,2], tab1[2,3], #fmsy
                               tab1[1,1], tab1[1,2], tab1[1,3], #bmsy
                               tab2[1,1], tab2[1,2], tab2[1,3], # b end
                               tab2[3,1], tab2[3,2], tab2[3,3], # bbmsy
                               tab2[2,1], tab2[2,2], tab2[2,3], # f end
                               tab2[4,1], tab2[4,2], tab2[4,3]) # ffmsy

  
  
  
  results = rbind (results, results_spict)
  return (results)
}



###########################################################################################
###########################################################################################
###########################################################################################
###                          Assessments                                                ###
###########################################################################################
###########################################################################################
###########################################################################################


# CMSY --------------------------------------------------------------------
cmsy = cmsy2(year = data$year, catch = data$catch, r.low = pb.r.low, r.hi = pb.r.hi)
plot_dlm (cmsy)
cmsy_ref = cmsy[["ref_pts"]]
cmsy_est = cmsy[["ref_ts"]]


results_cmsy = data.frame(Stock, "CMSY", 1, "Y", "Y", "Y", "Y",pb.resilience, 
                          pb.r.low, pb.r.hi,  pb.k.low, pb.k.hi, 
                          pb.log.k, pb.log.k.sd, 
                          pb.q.low, pb.q.hi, pb.log.q, pb.log.q.sd, 
                          cmsy_ref[1,2], cmsy_ref[1,3], cmsy_ref[1,4],
                          cmsy_ref[2,2], cmsy_ref[2,3], cmsy_ref[2,4],
                          cmsy_ref[3,2], cmsy_ref[3,3], cmsy_ref[3,4],
                          cmsy_ref[4,2], cmsy_ref[4,3], cmsy_ref[4,4],
                          cmsy_ref[5,2], cmsy_ref[5,3], cmsy_ref[5,4],
                          cmsy_est [ nrow(cmsy_est), c(4:9, 13:15, 19:21)],
                          "",  "", "")

results_cmsy = setNames(results_cmsy, c("Stock", "Method", "Index","FixN2",
                                        "Use.r", "Use.k", "Use.q",
                                        "Resilience", "Prior.r.low", "Prior.r.hi", 
                                        "Prior.k.low", "Prior.k.hi", 
                                        "Prior.log.k", "Prior.log.k.sd", 
                                        "Prior.q.low", "Prior.q.hi", "Prior.log.q", "Prior.log.q.sd", 
                             "r", "r.low", "r.high",
                             "k", "k.low", "k.high", 
                             "msy", "msy.low", "msy.high",
                             "fmsy", "fmsy.low", "fmsy.high",
                             "bmsy", "bmsy.low", "bmsy.high",
                             "b.end", "b.end.low", "b.end.hi",
                             "bbmsy", "bbmsy.low", "bbmsy.hi",
                             "f.end", "f.end.low", "f.end.hi",
                             "ffmsy", "ffmsy.low", "ffmsy.hi",
                             "Converge", "Prior.Warn", "Fit.Warn"))


# BSM ---------------------------------------------------------------------
### This will be using the r prior and q prior
### Not certain how it uses the k prior? Seems to be used as a start point?

bsm = bsm(year = data$year, catch = data$catch,
              biomass=data$index1, btype="CPUE",
              r.low = pb.r.low, r.hi = pb.r.hi)

#plot_dlm (bsm)
bsm_ref = bsm[["ref_pts"]]
bsm_est = bsm[["ref_ts"]]



results_bsm = data.frame(Stock, "BSM", 1, "Y", "Y", "Y", "Y",pb.resilience, 
                          pb.r.low, pb.r.hi,  pb.k.low, pb.k.hi, 
                          pb.log.k, pb.log.k.sd, 
                          pb.q.low, pb.q.hi, pb.log.q, pb.log.q.sd, 
                          bsm_ref[1,2], bsm_ref[1,3], bsm_ref[1,4],
                          bsm_ref[2,2], bsm_ref[2,3], bsm_ref[2,4],
                          bsm_ref[3,2], bsm_ref[3,3], bsm_ref[3,4],
                          bsm_ref[4,2], bsm_ref[4,3], bsm_ref[4,4],
                          bsm_ref[5,2], bsm_ref[5,3], bsm_ref[5,4],
                          bsm_est [ nrow(bsm_est), c(4:9, 13:15, 19:21)],
                         "", "","")

results_bsm = setNames(results_bsm, c("Stock", "Method", "Index","FixN2",
                                      "Use.r", "Use.k", "Use.q",
                                      "Resilience", "Prior.r.low", "Prior.r.hi", 
                                       "Prior.k.low", "Prior.k.hi", 
                                       "Prior.log.k", "Prior.log.k.sd", 
                                       "Prior.q.low", "Prior.q.hi", "Prior.log.q", "Prior.log.q.sd", 
                                        "r", "r.low", "r.high",
                                        "k", "k.low", "k.high", 
                                        "msy", "msy.low", "msy.high",
                                        "bmsy", "bmsy.low", "bmsy.high",
                                        "fmsy", "fmsy.low", "fmsy.high",
                                        "b.end", "b.end.low", "b.end.hi",
                                        "bbmsy", "bbmsy.low", "bbmsy.hi",
                                        "f.end", "f.end.low", "f.end.hi",
                                        "ffmsy", "ffmsy.low", "ffmsy.hi",
                                      "Converge", "Prior.Warn", "Fit.Warn"))

results = rbind (results_cmsy, results_bsm)


# SPICT Option 1 -------------------------------------------------------------------
### No priors and n not fixed at 2 - one abundance index
Use.Index = 1
Use.n.prior = "N"
Use.r.prior = "N"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 2 -------------------------------------------------------------------
### No priors but n = 2
Use.Index = 1
Use.n.prior = "Y"
Use.r.prior = "N"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 3 -------------------------------------------------------------------
### r priors and n = 2
Use.Index = 1
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 4 -------------------------------------------------------------------
### r and k priors and n = 2
Use.Index = 1
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 5 -------------------------------------------------------------------
### r and q priors and n = 2
Use.Index = 1
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 6 -------------------------------------------------------------------
### r, k, and q priors and n = 2
Use.Index = 1
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 7 -------------------------------------------------------------------
Use.Index = 1
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 8 -------------------------------------------------------------------
Use.Index = 1
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 9 -------------------------------------------------------------------
Use.Index = 1
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 10 -------------------------------------------------------------------
Use.Index = 1
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
#### Options using all available indices
# SPICT Option 1 -------------------------------------------------------------------
### No priors and n not fixed at 2 - one abundance index
Use.Index = "All"
Use.n.prior = "N"
Use.r.prior = "N"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 2 -------------------------------------------------------------------
### No priors but n = 2
Use.Index = "All"
Use.n.prior = "Y"
Use.r.prior = "N"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 3 -------------------------------------------------------------------
### r priors and n = 2
Use.Index = "All"
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 4 -------------------------------------------------------------------
### r and k priors and n = 2
Use.Index ="All"
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 5 -------------------------------------------------------------------
### r and q priors and n = 2
Use.Index = "All"
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 6 -------------------------------------------------------------------
### r, k, and q priors and n = 2
Use.Index = "All"
Use.n.prior = "Y"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 7 -------------------------------------------------------------------
Use.Index = "All"
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 8 -------------------------------------------------------------------
Use.Index = "All"
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "N"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)


#####################################################################################
# SPICT Option 9 -------------------------------------------------------------------
Use.Index = "All"
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "N"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

#####################################################################################
# SPICT Option 10 -------------------------------------------------------------------
Use.Index = "All"
Use.n.prior = "N"
Use.r.prior = "Y"
Use.k.prior = "Y"
Use.q.prior = "Y"

results = pb.spict(data, No.Index, Use.Index, Use.n.prior, Use.r.prior, Use.k.prior, Use.q.prior,
                   pb.log.r, pb.log.r.sd, pb.log.k, pb.log.k.sd,
                   pb.log.q, pb.log.q.sd)

######################################################################################
### Save results ######################################
write.csv(results, (paste("Results\\", format(Sys.Date(),format="%d%m%y"), "Results", Stock,  ".csv")))

          