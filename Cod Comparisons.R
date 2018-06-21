library (datalimited2)
library(spict)
options(scipen = 999)


cod.data = read.csv("Data//Cod_7_ek//Cod.data.csv")
cod.data$catch = cod.data$catch / 1000

# Calculate moving average
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

pb.resilience = "Medium"
pb.r.low = 0.2
pb.r.hi = 0.8
pb.start.r = c(pb.r.low, pb.r.hi)

### Mean of r range specified by the resilience
pb.log.r = log(0.5)
pb.log.r.sd = abs ((pb.log.r - log(pb.r.low))/2)


nyr = length(cod.data$year)
ct <- ma(ct.raw)

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

pb.endbio = pb.endbio_prior(nyr, cod.data$catch, ct)

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
mean.last.ct      <-mean(cod.data$catch[(nyr-lyr):nyr],na.rm=T) # get mean catch of last years
mean.last.cpue    <-mean(cod.data$index1[(nyr-lyr):nyr],na.rm=T) # get mean of CPUE of last years

gm.start.r      <- exp(mean(log(pb.start.r))) # get geometric mean of prior r range
if(mean(pb.endbio) >= 0.5) {  # if biomass is high
  q.1           <- mean.last.cpue*0.25*gm.start.r/mean.last.ct
  q.2           <- mean.last.cpue*0.5*start.r[2]/mean.last.ct
} else {
  q.1           <- mean.last.cpue*0.5*gm.start.r/mean.last.ct
  q.2           <- mean.last.cpue*start.r[2]/mean.last.ct
  }

pb.q.prior         <- c(q.1,q.2)
pb.q.low = pb.q.prior[1]
pb.q.hi = pb.q.prior[2]
pb.log.q = log(mean(pb.q.prior))
pb.log.q.sd = (pb.log.q -log(pb.q.low))/4

###########################################################################################
###########################################################################################
###########################################################################################
###                          Assessments                                                ###
###########################################################################################
###########################################################################################
###########################################################################################


# CMSY --------------------------------------------------------------------
cod.cmsy = cmsy2(year = cod.data$year, catch = cod.data$catch, r.low = pb.r.low, r.hi = pb.r.hi)
plot_dlm (cod.cmsy)
cmsy_ref = cod.cmsy[["ref_pts"]]



# BSM ---------------------------------------------------------------------
### This will be using the r prior and q prior
### Not certain how it uses the k prior? Seems to be used as a start point?

cod.bsm = bsm(year = cod.data$year, catch = cod.data$catch,
              biomass=cod.data$index1, btype="CPUE",
              resilience = pb.resilience)

#plot_dlm (cod.bsm)
bsm_ref = cod.bsm[["ref_pts"]]


# SPICT -------------------------------------------------------------------
inp <- list()
inp$obsC <- cod.data$catch
inp$timeC <- cod.data$year
### Use both abundance index
# inp$obsI[[1]] <- cod.data$index1 [30:47]
# inp$obsI[[2]] <- cod.data$index2 [33:47]
# inp$timeI <- list(cod.data$year[30:47], cod.data$year[33:47])

### Use only one abundance index
inp$obsI <- cod.data$index1 [30:47]
inp$timeI <- cod.data$year[30:47]

### provide r prior
inp$priors$logr <- c(pb.log.r, pb.log.r.sd, 1)
### Does not seem to like having k prior
inp$priors$logK <- c(pb.log.k, pb.log.k.sd, 1)
### provide prior for q
inp$priors$logq <- c(pb.log.q, pb.log.q.sd, 1)


### fix n=2
inp$priors$logn <- c(log(2), 1e-3)

inp <- check.inp(inp)

#plots the Catch and index 
plotspict.data(inp)

# Run the model
res <- fit.spict(inp)

# What are the calculated values
summary(res)

# plot everything
# plot(res)
# # Calculate the residuals
# res <- calc.osa.resid(res)
# # are the residuals ok?
# plotspict.diagnostic(res)
# sumspict.diagnostics(res)

