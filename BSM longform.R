library (datalimited2)
library(spict)
options(scipen = 999)


cod.data = read.csv("Data//Cod_7_ek//Cod.data.csv")

# Priors ------------------------------------------------------------------
pb.resilience = "Medium"
pb.r.low = 0.2
pb.r.hi = 0.8


# Functions Needed --------------------------------------------------------

# Calculate moving average
ma <- function(x){
  x.1 <- stats::filter(x,rep(1/3,3),sides=1)
  x.1[1] <- x[1]
  x.1[2] <- (x[1]+x[2])/2
  return(x.1)
}

# Prior Setting

# Set R prior
r_prior <- function(r.low, r.hi, res){
  # initial range of r from input file
  if(is.na(r.low)==F & is.na(r.hi)==F) {
    start.r <- c(r.low,r.hi)
  }else{
    # initial range of r based on resilience
    if(res == "High") {
      start.r <- c(0.6,1.5)} else if(res == "Medium") {
        start.r <- c(0.2,0.8)}    else if(res == "Low") {
          start.r <- c(0.05,0.5)}  else { # i.e. res== "Very low"
            start.r <- c(0.015,0.1)}
  }
  return(start.r)
}

# Set start saturation prior
startbio_prior <- function(stb.low, stb.hi, start.yr){
  # use initial biomass range from input file if stated
  if(is.na(stb.low)==F & is.na(stb.hi)==F){
    startbio <- c(stb.low,stb.hi)
  }else{
    # if start year < 1960 assume high biomass
    if(start.yr < 1960){
      startbio <- c(0.5,0.9)
    }else{
      # else use medium prior biomass range
      startbio <- c(0.2,0.6)
    }
  }
  return(startbio)
}

# Set intermediate saturation prior
intbio_prior <- function(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct){
  # get index of years with lowest and highest catch between start+3 and end-3 years
  min.yr.i <- which.min(ct[4:(length(ct)-3)])+3
  max.yr.i <- which.max(ct[4:(length(ct)-3)])+3
  min.ct <- ct[min.yr.i]
  max.ct <- ct[max.yr.i]
  # use year and biomass range for intermediate biomass from input file
  if(is.na(intb.low)==F & is.na(intb.hi)==F){
    int.yr   <- int.yr
    intbio   <- c(intb.low,intb.hi)
    # if contrast in catch is low, use initial range again in mid-year
  }else if(min(ct)/max(ct) > 0.6) {
    int.yr    <- as.integer(mean(c(start.yr, end.yr)))
    intbio    <- startbio
    # else if year of minimum catch is after max catch then use min catch
  }else if(min.yr.i > max.yr.i) {
    int.yr    <- yr[min.yr.i-1]
    if(startbio[1]>=0.5 &  (int.yr-start.yr) < (end.yr-int.yr) &
       (min.ct/max.ct) > 0.3) intbio <- c(0.2,0.6) else intbio <- c(0.01,0.4)
       # else use max catch
  } else {
    # assume that biomass range in year before maximum catch was high or medium
    int.yr    <- yr[max.yr.i-1]
    intbio    <- if((startbio[1]>=0.5 & (int.yr-start.yr) < (end.yr-int.yr))| # if initial biomass is high, assume same for intermediate
                    # ((min.ct/max.ct < 0.3 & (max.yr.i - min.yr.i) < 25))) c(0.5,0.9) else c(0.2,0.6) }
                    (((max.ct-min.ct)/max.ct)/(max.yr.i-min.yr.i) > 0.04)) c(0.5,0.9) else c(0.2,0.6) } # if incease is steep, assume high, else medium
  out <- list(intbio, int.yr)
  return(out)
}

# Set end saturation prior
endbio_prior <- function(endb.low, endb.hi, nyr, ct.raw, ct){
  # final biomass range from input file
  if(is.na(endb.low)==F & is.na(endb.hi)==F){
    endbio   <- c(endb.low,endb.hi)
  }else{
    # else use mean final catch/max catch to estimate final biomass
    rawct.ratio=ct.raw[nyr]/max(ct)
    endbio  <- if(ct[nyr]/max(ct) > 0.8) {c(0.4,0.8)} else if(rawct.ratio < 0.5) {c(0.01,0.4)} else {c(0.2,0.6)}
    
    # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
    if(endbio[2]==0.4){
      if(rawct.ratio< 0.05) {endbio[2] <- 0.1} else
        if(rawct.ratio< 0.15) {endbio[2] <- 0.2} else
          if(rawct.ratio< 0.35) {endbio[2] <- 0.3} else {endbio[2] <- 0.4}
    }
  }
  return(endbio)
}

# Set K prior
k_prior <- function(endbio, start.r, ct){
  # initial prior range of k values, assuming min k will be larger than max catch / prior for r
  if(mean(endbio) <= 0.5){
    start.k <- c(max(ct)/start.r[2],4*max(ct)/start.r[1])
  }else{
    start.k <- c(2*max(ct)/start.r[2],12*max(ct)/start.r[1])
  }
  return(start.k)
}


# BSM Analysis ------------------------------------------------------------


year = cod.data$year
catch = cod.data$catch
biomass = cod.data$index2
btype = "CPUE"
resilience = pb.resilience
r.low=NA
r.hi=NA 
stb.low=NA 
stb.hi=NA 
int.yr=NA
intb.low=NA 
intb.hi=NA
endb.low=NA 
endb.hi=NA 
q.start=NA
q.end=NA
verbose=T

    
# Display 3 digits
options(digits=3)
  
  # Set model parameters
  FullSchaefer <- F # will automatically change to TRUE if enough abundance data available
  dataUncert <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
  sigmaR <- 0.1 # overall process error for CMSY; SD=0.1 is the default
  n <- 10000 # initial number of r-k pairs
  n.new <- n # initialize n.new
  ni <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
  nab <- 5 # default=5; minimum number of years with abundance data to run BSM
  duncert <- dataUncert # global defaults for uncertainty
  sigR <- sigmaR # global defaults for uncertainty
  
  # Setup data
  ##############################################################################
  
  # Build catch data (using original cMSY variable naming convention)
  catchData <- data.frame(yr=year, ct=catch, bt=biomass)
  
    # Transform catch data
  # 1. Convert to 1000s tons (or other units)
  # 2. Calculate 3-yr moving average (average of past 3 years)
  ct.raw <- catchData$ct / 1000
  ct <- ma(ct.raw)
  
  # Transform biomass data
  bt <- catchData$bt / 1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
  
  # Identify number of years and start/end years
  yr <- catchData$yr # functions use this quantity
  nyr <- length(yr)
  start.yr <- min(yr)
  end.yr <- max(yr)
  
  # Determine initial ranges for parameters and biomass
  ##############################################################################
  
  # Set priors
  res <- resilience # rename resilience
  start.r <- r_prior(r.low, r.hi, res)
  startbio <- startbio_prior(stb.low, stb.hi, start.yr)
  int_params <- intbio_prior(intb.low, intb.hi, int.yr, start.yr, end.yr, startbio, yr, ct)
  intbio <- int_params[[1]]
  int.yr <- int_params[[2]]
  endbio <- endbio_prior(endb.low, endb.hi, nyr, ct.raw, ct)
  start.k <- k_prior(endbio, start.r, ct)
  
  # Record priors into dataframe
  priors <- data.frame(cbind(c("r", "k", "startbio", "intbio", "endbio"), source="default",
                             rbind(start.r, start.k, startbio, intbio, endbio)), year=NA, stringsAsFactors=F)
  colnames(priors) <- c("param", "source", "lo", "hi", "year")
  rownames(priors) <- NULL
  priors$year[priors$param=="intbio"] <- int.yr
  priors$lo <- as.numeric(priors$lo)
  priors$hi <- as.numeric(priors$hi)
  if(!is.na(r.low)){priors$source[priors$param=="r"] <- "expert"}
  if(!is.na(stb.low)){priors$source[priors$param=="startbio"] <- "expert"}
  if(!is.na(intb.low)){priors$source[priors$param=="intbio"] <- "expert"}
  if(!is.na(endb.low)){priors$source[priors$param=="endbio"] <- "expert"}
  
  # Print priors (if desired)
  if(verbose==T){
    cat("startbio=",startbio,ifelse(is.na(stb.low)==T,"default","expert"),
        ", intbio=",int.yr,intbio,ifelse(is.na(intb.low)==T,"default","expert"),
        ", endbio=",endbio,ifelse(is.na(endb.low)==T,"default","expert"),"\n")
  }
  
  
  # Bayesian analysis of catch & biomass (or CPUE) with Schaefer model
  ##############################################################################
  
  # Indicate that BSM is being fit
  FullSchaefer <- T
  
  # Set inits for r-k in lower right corner of log r-k space
  init.r      <- start.r[1]+0.8*(start.r[2]-start.r[1])
  init.k      <- start.k[1]+0.1*(start.k[2]-start.k[1])
  
  # Vector with no penalty (=0) if predicted biomass is within viable range, else a penalty of 10 is set
  pen.bk = pen.F = rep(0,length(ct))
  
  # Add biomass priors
  b.yrs = c(1,length(start.yr:int.yr),length(start.yr:end.yr))
  b.prior = rbind(matrix(c(startbio[1],startbio[2],intbio[1],intbio[2],endbio[1],endbio[2]),2,3),rep(0,3)) # last row includes the 0 pen
  
  # Print
  if(verbose==T){cat("Running MCMC analysis....\n")}
  
  
  # Setup CPUE-based BSM
  ##########################################################
    
  # CPUE BSM

  # Catchability stuff
  ########################################
  
  # Expert-specified catchability (q)
  if(is.na(q.start)==F & is.na(q.end)==F){
    mean.last.ct      <-mean(ct[yr >= q.start & yr <= q.end], na.rm=T) # get mean catch of indicated years
    mean.last.cpue    <-mean(bt[yr >= q.start & yr <= q.end], na.rm=T) # get mean of CPUE of indicated years
    # Default catchability (q)
    # get prior range for q from mean catch and mean CPUE in recent years
  }else{
    lyr               <- ifelse(mean(start.r)>=0.5,5,10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish
    mean.last.ct      <-mean(ct[(nyr-lyr):nyr],na.rm=T) # get mean catch of last years
    mean.last.cpue    <-mean(bt[(nyr-lyr):nyr],na.rm=T) # get mean of CPUE of last years
  }
  gm.start.r      <- exp(mean(log(start.r))) # get geometric mean of prior r range
  if(mean(endbio) >= 0.5) {  # if biomass is high
    q.1           <- mean.last.cpue*0.25*gm.start.r/mean.last.ct
    q.2           <- mean.last.cpue*0.5*start.r[2]/mean.last.ct
  } else {
    q.1           <- mean.last.cpue*0.5*gm.start.r/mean.last.ct
    q.2           <- mean.last.cpue*start.r[2]/mean.last.ct
  }
  q.prior         <- c(q.1,q.2)
  init.q          <- mean(q.prior)
  
  # Setup JAGS model
  ########################################
  
  # Data to be passed on to JAGS
  jags.data        <- c('ct','bt','nyr', 'start.r', 'start.k', 'startbio', 'q.prior',
                        # 'init.q','init.r','init.k',
                        'pen.bk','pen.F','b.yrs','b.prior')
  # Parameters to be returned by JAGS
  jags.save.params <- c('r','k','q', 'P')
  
  # JAGS model
  Model = "model{
  # to reduce chance of non-convergence, Pmean[t] values are forced >= eps
  eps<-0.01
  penm[1] <- 0 # no penalty for first biomass
  Pmean[1] <- log(alpha)
  P[1] ~ dlnorm(Pmean[1],itau2)
  for (t in 2:nyr) {
  Pmean[t] <- ifelse(P[t-1] > 0.25,
  log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps)),  # Process equation
  log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
  P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
  penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*k*P[t])-log(q*k*(eps+0.001)),ifelse(P[t]>1,log(q*k*P[t])-log(q*k*(0.99)),0)) # penalty if Pmean is outside viable biomass
  }
  # ><> Biomass priors/penalties are enforced as follows
  for (i in 1:3) {
  penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[2,i]),0))
  b.prior[3,i] ~ dnorm(penb[i],100)
  }
  for (t in 1:nyr){
  Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) #><> Penalty term on F > 1, i.e. ct>B
  pen.F[t]  ~ dnorm(Fpen[t],1000)
  pen.bk[t] ~ dnorm(penm[t],10000)
  cpuem[t]  <- log(q*P[t]*k);
  bt[t]     ~ dlnorm(cpuem[t],isigma2);
  }
  # priors
  log.alpha               <- log((startbio[1]+startbio[2])/2) # needed for fit of first biomass
  sd.log.alpha            <- (log.alpha-log(startbio[1]))/4
  tau.log.alpha           <- pow(sd.log.alpha,-2)
  alpha                   ~  dlnorm(log.alpha,tau.log.alpha)
  # search in the k space starting from 20% of the range
  log.km              <- log(start.k[1]+0.2*(start.k[2]-start.k[1]))
  sd.log.k            <- (log.km-log(start.k[1]))/4
  tau.log.k           <- pow(sd.log.k,-2)
  k                   ~  dlnorm(log.km,tau.log.k)
  # set realistic prior for q
  log.qm              <- mean(log(q.prior))
  sd.log.q            <- (log.qm-log(q.prior[1]))/4
  tau.log.q           <- pow(sd.log.q,-2)
  q                   ~  dlnorm(log.qm,tau.log.q)
  # define process (tau) and observation (sigma) variances as inversegamma prios
  itau2 ~ dgamma(4,0.01)
  tau2  <- 1/itau2
  tau   <- pow(tau2,0.5)
  isigma2 ~ dgamma(2,0.01)
  sigma2 <- 1/isigma2
  sigma <- pow(sigma2,0.5)
  log.rm              <- mean(log(start.r))
  sigma.log.r         <- abs(log.rm - log(start.r[1]))/2
  tau.log.r           <- pow(sigma.log.r,-2)
  r                   ~  dlnorm(log.rm,tau.log.r)
}" # end of JAGS model for CPUE


  
  # Run JAGS model
  ##########################################################
  
  # Write JAGS model to temporary file
  wd <- getwd()
  jags_model <- paste(wd, "r2jags.bug", sep="/")
  cat(Model, file=jags_model)
  
  # Initialize JAGS model?
  if(btype=="biomass") {
    j.inits     <- function(){list("r"=stats::rnorm(1,mean=init.r,sd=0.2*init.r),
                                   "k"=stats::rnorm(1,mean=init.k,sd=0.1*init.k),
                                   "itau2"=1000,
                                   "isigma2"=1000)}} else {
                                     j.inits <- function(){list("r"=stats::rnorm(1,mean=init.r,sd=0.2*init.r),
                                                                "k"=stats::rnorm(1,mean=init.k,sd=0.1*init.k),
                                                                "q"=stats::rnorm(1,mean=init.q,sd=0.2*init.q),
                                                                "itau2"=1000,
                                                                "isigma2"=1000)}}
  # Run JAGS model
  jags_outputs <- R2jags::jags(data=jags.data,
                               working.directory=wd, inits=j.inits,
                               parameters.to.save=jags.save.params,
                               model.file="r2jags.bug", #n.chains = 2,
                               n.burnin = 30000, n.thin = 10,
                               n.iter = 60000)
  
  
  # Extract JAGS model results
  ##########################################################
  
  r_raw            <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$r))
  k_raw            <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$k))
  # Importance sampling: only accept r-k pairs where r is near the prior range
  r_out            <- r_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]
  k_out            <- k_raw[r_raw > 0.5*start.r[1] & r_raw < 1.5 * start.r[2]]
  
  mean.log.r.jags  <- mean(log(r_out))
  sd.log.r.jags    <- stats::sd(log(r_out))
  r.jags           <- exp(mean.log.r.jags)
  lcl.r.jags       <- exp(mean.log.r.jags - 1.96*sd.log.r.jags)
  ucl.r.jags       <- exp(mean.log.r.jags + 1.96*sd.log.r.jags)
  mean.log.k.jags  <- mean(log(k_out))
  sd.log.k.jags    <- stats::sd(log(k_out))
  k.jags           <- exp(mean.log.k.jags)
  lcl.k.jags       <- exp(mean.log.k.jags - 1.96*sd.log.k.jags)
  ucl.k.jags       <- exp(mean.log.k.jags + 1.96*sd.log.k.jags)
  MSY.posterior     <- r_out*k_out/4 # simpler
  mean.log.MSY.jags <- mean(log(MSY.posterior))
  sd.log.MSY.jags   <- stats::sd(log(MSY.posterior))
  MSY.jags          <- exp(mean.log.MSY.jags)
  lcl.MSY.jags      <- exp(mean.log.MSY.jags - 1.96*sd.log.MSY.jags)
  ucl.MSY.jags      <- exp(mean.log.MSY.jags + 1.96*sd.log.MSY.jags)
  
  # CPUE-based computations
  if(btype=="CPUE") {
    q_out           <- as.numeric(coda::mcmc(jags_outputs$BUGSoutput$sims.list$q))
    mean.log.q      <- mean(log(q_out))
    sd.log.q        <- stats::sd(log(q_out))
    mean.q          <- exp(mean.log.q)
    lcl.q           <- exp(mean.log.q-1.96*sd.log.q)
    ucl.q           <- exp(mean.log.q+1.96*sd.log.q)
    F.bt.cpue       <- mean.q*ct.raw/bt
    Fmsy.cpue       <- r.jags/2
  }
  
  # get F from observed biomass
  if(btype == "biomass"){
    F.bt       <- ct.raw/bt
    Fmsy.bt    <- r.jags/2
  }
  
  # get relative biomass P=B/k as predicted by BSM, including predictions for years with NA abundance
  all.P    <- jags_outputs$BUGSoutput$sims.list$P # matrix with P distribution by year
  quant.P  <- apply(all.P,2,stats::quantile,c(0.025,0.5,0.975),na.rm=T)
  
  # get k, r posterior ><>
  all.k  <- jags_outputs$BUGSoutput$sims.list$k # matrix with P distribution by year
  all.r  <- jags_outputs$BUGSoutput$sims.list$r # matrix with P distribution by year
  
  # get B/Bmys posterior
  all.b_bmsy=NULL
  for(t in 1:ncol(all.P)){
    all.b_bmsy  <- cbind(all.b_bmsy,all.P[,t]*2)}
  
  # get F/Fmys posterior ><>
  all.F_Fmsy=NULL
  for(t in 1:ncol(all.P)){
    all.F_Fmsy<- cbind(all.F_Fmsy,(ct.raw[t]/(all.P[,t]*all.k))/ifelse(all.P[,t]>0.25,all.r/2,all.r/2*4*all.P[,t]))}
  
  # Organize results
  ##############################################################################
  
  # Get management results
  MSY   <-MSY.jags; lcl.MSY<-lcl.MSY.jags; ucl.MSY<-ucl.MSY.jags
  Bmsy  <-k.jags/2; lcl.Bmsy<-lcl.k.jags/2; ucl.Bmsy<-ucl.k.jags/2
  Fmsy  <-r.jags/2; lcl.Fmsy<-lcl.r.jags/2; ucl.Fmsy<-ucl.r.jags/2
  B.Bmsy<-2*quant.P[2,];lcl.B.Bmsy<-2*quant.P[1,];ucl.B.Bmsy<-2*quant.P[3,]
  B          <-B.Bmsy*Bmsy;lcl.B<-lcl.B.Bmsy*Bmsy;ucl.B<-ucl.B.Bmsy*Bmsy
  Fm           <- ct.raw/B;lcl.F<-ct.raw/ucl.B;ucl.F<-ct.raw/lcl.B
  Fmsy.vec     <- ifelse(B.Bmsy>0.5,Fmsy,Fmsy*2*B.Bmsy)
  lcl.Fmsy.vec <- ifelse(B.Bmsy>0.5,lcl.Fmsy,lcl.Fmsy*2*B.Bmsy)
  ucl.Fmsy.vec <- ifelse(B.Bmsy>0.5,ucl.Fmsy,ucl.Fmsy*2*B.Bmsy)
  F.Fmsy       <- Fm/Fmsy.vec; lcl.F.Fmsy<-lcl.F/Fmsy.vec; ucl.F.Fmsy<-ucl.F/Fmsy.vec
  
  # Print results (if desired)
  if(verbose==T){
    # Print priors
    cat("---------------------------------------\n")
    cat("Catch data used from years", min(yr),"-", max(yr),", abundance =", btype, "\n")
    cat("Prior initial relative biomass =", startbio[1], "-", startbio[2],ifelse(is.na(stb.low)==T,"default","expert"), "\n")
    cat("Prior intermediate rel. biomass=", intbio[1], "-", intbio[2], "in year", int.yr,ifelse(is.na(intb.low)==T,"default","expert"), "\n")
    cat("Prior final relative biomass   =", endbio[1], "-", endbio[2],ifelse(is.na(endb.low)==T,"default","expert"), "\n")
    cat("Prior range for r =", format(start.r[1],digits=2), "-", format(start.r[2],digits=2),ifelse(is.na(r.low)==T,"default","expert,"),
        ", prior range for k =", start.k[1], "-", start.k[2],"\n")
    # Print BSM results
    cat("Results from Bayesian Schaefer model (BSM) using catch &",btype,"\n")
    cat("------------------------------------------------------------\n")
    if(btype == "CPUE") cat("q   =", mean.q,", lcl =", lcl.q, ", ucl =", ucl.q,"\n")
    cat("r   =", r.jags,", 95% CL =", lcl.r.jags, "-", ucl.r.jags,", k =", k.jags,", 95% CL =", lcl.k.jags, "-", ucl.k.jags,"\n")
    cat("MSY =", MSY.jags,", 95% CL =", lcl.MSY.jags, "-", ucl.MSY.jags,"\n")
    cat("Relative biomass in last year =", quant.P[2,][nyr], "k, 2.5th perc =",quant.P[1,][nyr],
        ", 97.5th perc =", quant.P[3,][nyr],"\n")
    cat("Exploitation F/(r/2) in last year =", (ct.raw[nyr]/(quant.P[2,][nyr]*k.jags))/(r.jags/2) ,"\n\n")
    # Print catchability prior (if CPUE-based BSM)
    if(btype=="CPUE") {
      cat("Prior range of q =",q.prior[1],"-",q.prior[2],"\n")
    }
  }
  
  # Format results for output
  ##############################################################################
  
  # Refence points dataframe
  ref_pts <- data.frame(rbind(c(est=r.jags, lo=lcl.r.jags, hi=ucl.r.jags),
                              c(k.jags*1000, lcl.k.jags*1000, ucl.k.jags*1000),
                              c(MSY.jags*1000, lcl.MSY.jags*1000, ucl.MSY.jags*1000),
                              c(Bmsy*1000, lcl.Bmsy*1000, ucl.Bmsy*1000),
                              c(Fmsy, lcl.Fmsy, ucl.Fmsy)))
  ref_pts$param <- c("r", "k", "msy", "bmsy", "fmsy")
  ref_pts <- subset(ref_pts, select=c(param, est, lo, hi))
  
  # # Format F time series data
  # if(btype=="CPUE"){f <- F.bt.cpue}else{f <- F.bt}
  # if(btype=="CPUE"){ffmsy <- F.bt.cpue/Fmsy.cpue}else{ffmsy <- F.bt/Fmsy.bt}
  #
  # # Format saturation data
  # if(btype=="CPUE"){s <-bt/(mean.q*k.jags)}else{s <-bt/k.jags}
  # if(btype=="CPUE"){s_lo <- bt/(mean.q*ucl.k.jags)}else{s_lo <- bt/ucl.k.jags}
  # if(btype=="CPUE"){s_hi <- bt/(mean.q*lcl.k.jags)}else{s_hi <- bt/lcl.k.jags}
  
  # Reference points time series
  ref_ts <- data.frame(year=yr, catch=ct.raw*1000, catch_ma=ct*1000,
                       b=B*1000, b_lo=lcl.B*1000, b_hi=ucl.B*1000,
                       bbmsy=B.Bmsy, bbmsy_lo=lcl.B.Bmsy, bbmsy_hi=ucl.B.Bmsy,
                       s=B.Bmsy/2, s_lo=lcl.B.Bmsy/2, s_hi=ucl.B.Bmsy/2,
                       f=Fm, f_lo=lcl.F, f_hi=ucl.F,
                       fmsy=Fmsy.vec, fmsy_lo=lcl.Fmsy.vec, fmsy_hi=ucl.Fmsy.vec,
                       ffmsy=F.Fmsy, ffmsy_lo=lcl.F.Fmsy, ffmsy_hi=ucl.F.Fmsy,
                       er=Fm/Fmsy)
  
  #stop parallel processing clusters
  # parallel::stopCluster(cl)
  # doParallel::stopImplicitCluster()
  
  # Assemble output
  output <- list(ref_pts=ref_pts, ref_ts=ref_ts, priors=priors,
                 r_viable=r_out, k_viable=k_out, method="BSM")
  
bsm_lonform = output[["ref_pts"]]
