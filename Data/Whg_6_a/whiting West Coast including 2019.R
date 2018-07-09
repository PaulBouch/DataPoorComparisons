
#################################################
# TSA run - final assessment 2 May 2017
#################################################

# rm(list=ls())

# source(file.path("..", "R functions", "tsa functions 2017 04 17.r"))

setwd("C:/My files/WGCSE/2018/TSA model/")
source("tsa functions 2017 04 17.r")

whiting_data <- tsa.input(
  landings = list(
    numbers = "whi.VIa.landings.numbers.dat", 
    weights = "whi.VIa.landings.weights.dat"
  ), 
  discards = list(
    numbers = "whi.VIa.discards.numbers.dat", 
    weights = "whi.VIa.discards.weights.dat"
  ), 
  survey = list(
    WCIBTS.Q1 = list(indices = "whi.VIa.IBTS.Q1.dat"),
    WCIBTS.Q4 = list(indices = "whi.VIa.IBTS.Q4.dat"),
    IRIBTS.Q4 = list(indices = "whi.VIa.IRGFS.Q4.dat"),
    SCO.Q1 = list(indices = "whi.VIa.IBTS.Q1new.dat"),
    SCO.Q4 = list(indices = "whi.VIa.IBTS.Q4new.dat")
  ),             
  auxiliary = list(
    natural.mortality = list(b = -0.29, a = 3.0 * (1000 ** -0.29)), 
    maturity = "whi.VIa.maturity.dat"
  )
)


whiting_setup <- tsa.setup(
  whiting_data,
  model = list(
    response = "discards", 
    model = "full", 
    ages = list(discards = 1:5), 
    scaling = -4, 
    years = c(1981:2018), 
    misrep = c(1995, 2005)
  ),
  surveys.control = list(
    WCIBTS.Q1 = list(scaling = -2),
    WCIBTS.Q4 = list(scaling = -2),
    IRIBTS.Q4 = list(ages = 1:4, scaling = -2, years = c(2003:2006, 2008:2017)),
    SCO.Q1 = list(scaling = -2),
    SCO.Q4 = list(scaling = -2, years = c(2011:2012, 2014:2017))
  ),
  ylim = c(1981, 2019), 
  f.est = c(1, 2, 4), 
  f.range = c(2, 4), 
  recruitment = "hockey"
)


whiting_setup <- tsa.update(
  whiting_setup, 
  gudmundssonH1 = c(2, 1, 1, 1, 1, 1, 1),
  survey = list(
    WCIBTS.Q1 = list(points = list(c(1992, 4, 5), c(1992, 5, 3), c(1993, 2, 3), c(2000, 1, 3), c(2000, 2, 3))),
    WCIBTS.Q4 = list(ages = c(1, 1, 1, 1, 1, 2), points = list(c(2007, 4, 3), c(2007, 5, 3)))
  ), 
  cvmult = list(
    landings=list(ages=c(2,1,1,1,1,1,2)), 
    discards = list(
      ages=c(1,1,1,1,2), 
      points = list(c(1987, 1, 3), c(1981, 1, 3), c(1991, 3, 3), c(2000, 1, 3), c(2013, 1, 3), c(2017, 5, 3)))
  )
)


whiting_setup$param <- within(whiting_setup$param, {

  fishing.selection["age 1", ]     <- c(0.10, 0.04, 0.2, 0.02, NA, 1)
  fishing.selection["age 2", ]     <- c(0.11, 0.05, 0.3, 0.01, NA, 1)
  fishing.selection["age 4", ]     <- c(0.33, 0.20, 0.6, 0.04, NA, 1)
  
  fishing.sd["F", ]                <- c(0.01, 0, 0.2, 0.05, NA, 1)  
  fishing.sd["U", ]                <- c(0.09, 0, 0.2, 0.02, NA, 1)
  fishing.sd["V", ]                <- c(0.01, 0, 0.2, 0.05, NA, 1)
  fishing.sd["Y", ]                <- c(0.28, 0, 0.4, 0.01, NA, 1)
  
  cv["landings", ]                 <- c(0.17, 0.05, 0.50, 0.02, NA, 1)
  cv["discards", ]                 <- c(0.53, 0.30, 0.90, 0.02, NA, 1)
  
  recruitment["max recruits", ]         <- c(29.0, 20.0, 40.00, 1.4, NA, 1)
  recruitment["ssb (change point)", ]   <- c(3.15, 2.0, 4.5, 0.29, NA, 1)
  recruitment["cv", ]                   <- c(0.32, 0.15, 0.40, 0.03, NA, 1)
  
  discard.rate.sd["transitory", ]  <- c(0.30, 0, 0.5, 0.05, NA, 1)
  discard.rate.sd["persistent", ]  <- c(0.19, 0, 0.3, 0.03, NA, 1)
  
  misrep["sd.transitory", ] <- c(0.01, 0, 0.4, 0.05, NA, 1)
  misrep["sd.persistent", ] <- c(0.18, 0, 0.4, 0.05, NA, 1)
  
  if (T)
    survey$WCIBTS.Q1 <- within(survey$WCIBTS.Q1, {
      selection["age 1", ] <- c(1.11, 0.7, 2.5, 0.14, 0.02, 1)
      selection["age 2", ] <- c(1.15, 0.6, 2.5, 0.14, 0.02, 1)        
      selection["age 3", ] <- c(0.98, 0.5, 2.5, 0.13, 0.02, 1)
      selection["age 4", ] <- c(0.82, 0.4, 2.5, 0.12, 0.02, 1)
      selection["age 5", ] <- c(0.68, 0.4, 1.5, 0.11, 0.02, 1)
      selection["age 6", ] <- c(0.60, 0.4, 1.5, 0.09, 0.02, 1)
      
      cv["sigma", ]       <- c(0.45, 0, 1.00, 0.03, 0.02, 1)
      cv["eta", ]         <- c(0.10, 0, 0.50, 0.03, 0.02, 1)
      cv["omega", ]       <- c(0.19, 0, 0.30, 0.05, 0.02, 1)
      cv["beta", ]        <- c(0.11, 0, 0.50, 0.03, 0.02, 1)
    })
  
  if (T)
    survey$WCIBTS.Q4 <- within(survey$WCIBTS.Q4, {
      
      selection["age 1", ] <- c(3.26, 2.5, 5.0, 0.28, NA, 1)
      selection["age 2", ] <- c(3.00, 2.5, 5.0, 0.22, NA, 1)
      selection["age 3", ] <- c(2.36, 2.0, 5.0, 0.19, NA, 1)
      selection["age 4", ] <- c(2.05, 1.5, 5.0, 0.19, NA, 1)
      selection["age 5", ] <- c(2.77, 2.0, 5.0, 0.46, NA, 1)
      selection["age 6", ] <- c(0.49, 0.1, 2.0, 0.27, NA, 1)
      
      cv["sigma", ]       <- c(0.20, 0, 0.50, 0.03, NA, 1)
      cv["eta", ]         <- c(0.19, 0, 0.50, 0.05, NA, 1)
      cv["omega", ]       <- c(0.01, 0, 0.30, 0.05, NA, 1)
      cv["beta", ]        <- c(0.14, 0, 0.50, 0.04, NA, 1)
    })
  
  if (T)
    survey$IRIBTS.Q4 <- within(survey$IRIBTS.Q4, {
      
      selection["age 1", ] <- c(11.7, 7, 25, 2.0, NA, 1)
      selection["age 2", ] <- c(12.0, 7, 25, 2.0, NA, 1)
      selection["age 3", ] <- c(13.1, 3, 20, 2.5, NA, 1)
      selection["age 4", ] <- c(10.9, 2, 20, 2.1, NA, 1)
    
      cv["sigma", ]       <- c(0.31, 0, 1.00, 0.05, NA, 1)
      cv["eta", ]         <- c(0.39, 0, 0.70, 0.05, NA, 1)
      cv["omega", ]       <- c(0.11, 0, 0.50, 0.05, NA, 1)
      cv["beta", ]        <- c(0.19, 0, 0.50, 0.05, NA, 1)
    })
  
  if (T)
    survey$SCO.Q1 <- within(survey$SCO.Q1, {
      
      selection["age 1", ] <- c(7.7, 1.0, 15.0, 2.2, NA, 1)
      selection["age 2", ] <- c(5.4, 1.0, 15.0, 0.9, NA, 1)        
      selection["age 3", ] <- c(9.7, 1.0, 15.0, 1.9, NA, 1)
      selection["age 4", ] <- c(7.5, 1.0, 15.0, 1.8, NA, 1)
      selection["age 5", ] <- c(7.3, 1.0, 15.0, 2.3, NA, 1)
      selection["age 6", ] <- c(6.3, 1.0, 15.0, 2.0, NA, 1)
      
      cv["sigma", ]       <- c(0.26, 0, 0.90, 0.05, NA, 1)
      cv["eta", ]         <- c(0.20, 0, 0.50, 0.05, NA, 1)
      cv["omega", ]       <- c(0.01, 0, 0.50, 0.05, NA, 1)
      cv["beta", ]        <- c(0.23, 0, 0.50, 0.05, NA, 1)
    })
  
  if (T)
    survey$SCO.Q4 <- within(survey$SCO.Q4, {
      
      selection["age 1", ] <- c( 2.9, 1.2, 10.0, 0.8, NA, 1)
      selection["age 2", ] <- c(10.8, 4.0, 15.0, 2.2, NA, 1)
      selection["age 3", ] <- c( 5.9, 2.0, 15.0, 1.3, NA, 1)
      selection["age 4", ] <- c( 7.5, 2.0, 15.0, 1.8, NA, 1)
      selection["age 5", ] <- c( 4.9, 1.0, 15.0, 1.6, NA, 1)
      selection["age 6", ] <- c( 7.1, 2.0, 15.0, 2.4, NA, 1)
      
      cv["sigma", ]       <- c(0.30, 0.2, 0.60, 0.05, NA, 1)
      cv["eta", ]         <- c(0.06, 0.0, 0.40, 0.05, NA, 1)
      cv["omega", ]       <- c(0.08, 0.0, 0.40, 0.05, NA, 1)
      cv["beta", ]        <- c(0.01, 0.0, 0.40, 0.05, NA, 1)
    })
  
})


# load("whiting fit mod.RData")

whiting_quick <- tsa.fit(whiting_setup, "TSA 2017 04 17.dll", filterType = "backward")
whiting_fit <- tsa.fit(whiting_setup, "TSA 2017 04 17.dll")

# save.image(file="whiting fit.RData")
# save.image(file="whiting fit mod.RData")


save(whiting_quick, file = file.path("backup fits", "backup May 2017", "whiting quick.RData"))
save(whiting_fit, file = file.path("backup fits", "backup May 2017", "whiting fit.RData"))
save.image(file="whiting fit including 2019.RData")
load("whiting fit including 2019.RData")

##### Retro analysis 

whiting_retro <- tsa.retro.fit(whiting_fit, 2008, "TSA 2017 04 17.dll", earlySurvey = c("WCIBTS.Q1", "SCO.Q1"))
# save(whiting_retro, file = file.path("backup fits", "backup May 2017", "whiting retro.RData"))
# load("whiting fit retro.RData")

tsa.retro.plot.stock.summary(whiting_retro)

# windows(9, 6)
# png(filename="figure.png", height=6, width=9, bg="white", pointsize=12, res=600, units="in")
# source("tsa retro functions HD.r")
# tsa.retro.plot.stock.summary(whiting_retro[paste(2008:2018)])
# dev.off()



##### TSA plots

tmp.wk <- whiting_fit
tsa.plot.stock.summary(tmp.wk)

tmp.wk$twice.loglik  
tmp.wk$summary

windows(5, 6)
# png(filename="figure.png", height=5, width=5, bg="white", pointsize=12, res=300, units="in") # report
# png(filename="figure.png", height=6, width=5, bg="white", pointsize=12, res=300, units="in") # presentation

tsa.plot.errors(tmp.wk$prediction.errors$landings)
tsa.plot.errors(tmp.wk$prediction.errors$discards)

tsa.plot.errors(tmp.wk$residuals$landings, "residuals")
tsa.plot.errors(tmp.wk$residuals$discards, "residuals")

windows(5, 6)
# png(filename="figure.png", height=5, width=5, bg="white", pointsize=12, res=300, units="in") # report
# png(filename="figure.png", height=6, width=5, bg="white", pointsize=12, res=300, units="in") # presentation

tsa.plot.errors(tmp.wk$prediction.errors$survey$WCIBTS.Q1)
tsa.plot.errors(tmp.wk$prediction.errors$survey$WCIBTS.Q4)
tsa.plot.errors(tmp.wk$prediction.errors$survey$IRIBTS.Q4)
tsa.plot.errors(tmp.wk$prediction.errors$survey$SCO.Q1)
tsa.plot.errors(tmp.wk$prediction.errors$survey$SCO.Q4)

tsa.plot.errors(tmp.wk$residuals$survey$WCIBTS.Q1, "residuals")
tsa.plot.errors(tmp.wk$residuals$survey$WCIBTS.Q4, "residuals")
tsa.plot.errors(tmp.wk$residuals$survey$IRIBTS.Q4, "residuals")
tsa.plot.errors(tmp.wk$residuals$survey$SCO.Q1, "residuals")
tsa.plot.errors(tmp.wk$residuals$survey$SCO.Q4, "residuals")

windows(9, 6)
# png(filename="figure.png", height=6, width=9, bg="white", pointsize=12, res=300, units="in")

tsa.plot.discards(tmp.wk)
tsa.plot.stock.recruit(tmp.wk)
tsa.plot.catchability(tmp.wk)
tsa.plot.stock.summary(tmp.wk)

# dev.off()











