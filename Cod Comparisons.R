library (datalimited2)
library(spict)
options(scipen = 999)


cod.data = read.csv("Data//Cod_7_ek//Cod.data.csv")



# Priors ------------------------------------------------------------------
pb.resilience = "Medium"
pb.r.low = 0.2
pb.r.hi = 0.8

pb.log.r = -0.3010
pb.log.r.sd = 0.04845

pb.log.k =  5.323
pb.log.k.sd = 0.232355763


# CMSY --------------------------------------------------------------------
cod.cmsy = cmsy2(year = cod.data$year, catch = cod.data$catch, r.low = pb.r.low, r.hi = pb.r.hi)
plot_dlm (cod.cmsy)
cmsy_ref = cod.cmsy[["ref_pts"]]



# BSM ---------------------------------------------------------------------
cod.bsm = bsm(year = cod.data$year, catch = cod.data$catch,
              biomass=cod.data$index2, btype="CPUE",
              resilience = pb.resilience)

#plot_dlm (cod.bsm)
bsm_ref = cod.bsm[["ref_pts"]]


# SPICT -------------------------------------------------------------------
inp <- list()
inp$obsC <- cod.data$catch
inp$timeC <- cod.data$year
### Use both abundance index
inp$obsI[[1]] <- cod.data$index1 [30:47]
inp$obsI[[2]] <- cod.data$index2 [33:47]
inp$timeI <- list(cod.data$year[30:47], cod.data$year[33:47])

### Use only one abundance index
# inp$obsI <- cod.data$index1 [30:47]
# inp$timeI <- cod.data$year[30:47]

### provide r prior
inp$priors$logr <- c(pb.log.r, pb.log.r.sd, 1)
### Does not seem to like having k prior
inp$priors$logK <- c(pb.log.k, pb.log.k.sd, 1)

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
plot(res)

# Calculate the residuals
res <- calc.osa.resid(res)

# are the residuals ok?
plotspict.diagnostic(res)

sumspict.diagnostics(res)

