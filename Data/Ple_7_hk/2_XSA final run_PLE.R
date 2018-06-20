###################################
## WGCSE PLE 27.7 h-k - XSA Assessment 
## Date Created:  "Thu May 04 15:25:49 2017"
## Date Edited: "Tue May 08 20:58:14 2018"
####################################

gc()
rm(list=ls())

# Using R Version 3.4.3
library(FLCore)         
library(FLEDA)          
library(FLAssess)       
library(FLXSA)          

setwd("\\\\GalwayFS03\\Fishdata\\Data for ICESWG\\2018\\WGCSE\\ple.27.7h-k\\XSA\\input data")

stock=readFLStock('PLE7jkID.txt')
tun=readFLIndices('PLE7jkTU.txt')

units(stock)[1:17]<- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

stock@catch.n = stock@landings.n
stock@catch.wt = stock@landings.wt
stock@catch = stock@landings

stock@stock.wt@.Data <- ifelse(stock@stock.wt@.Data==0,999,stock@stock.wt@.Data)

landings(stock) #this should be the 7jk landings only (from intercatch)
computeLandings(stock) #this should be the 7jk landings only (from intercatch)
#ok

#trim 
stock <- trim(stock,age=4:10)

#set plusgroup
stock=setPlusGroup(stock,plusgroup=8)

tun[[1]] <- trim(tun[[1]],age=4:8)


#run XSA
xsa.control <- FLXSA.control(tol = 1e-09, maxit = 99,  min.nse = 0.3,  fse  = 1.0,
                              rage = -1,   qage  = 6,   shk.n   = TRUE, shk.f = TRUE,
                              shk.yrs = 5, shk.ages= 3, window  = 100,  tsrange = 99,
                              tspower = 0)
xsa <- FLXSA(stock, tun, xsa.control)
xsa@desc <- 'Run1'
stock=stock+xsa
stock@desc=xsa@desc
run1=stock
xsa1 <- xsa

#bug fix
xsa@index.range[[1]][2] <- 6
diag <- capture.output(diagnostics(xsa))
write.csv(diag,'diag.csv',row.names=F)

write.csv(xsa@harvest@.Data,'F.csv')
write.csv(xsa@stock.n@.Data,'N.csv')
Ntab <- xsa@stock.n@.Data
Ftab <- xsa@harvest@.Data
Mtab <- stock@m@.Data
i <- dim(Ntab)[2]
x <- exp(-(Ftab+Mtab)[,20,,,,])*Ntab[,20,,,,]
j <- dim(Ntab)[1]
x1 <- c(0,x[1:j-1])
x1[j] <- x1[j]+x[j]
write.csv(x1,'N_int.csv',row.names=F)


xsa.control <- FLXSA.control(tol = 1e-09, maxit = 99,  min.nse = 0.5,  fse  = 2.0,
                              rage = -1,   qage  = 7,   shk.n   = TRUE, shk.f = TRUE,
                              shk.yrs = 5, shk.ages= 3, window  = 100,  tsrange = 99,
                              tspower = 0)
xsa <- FLXSA(stock, tun, xsa.control)
xsa@desc <- 'Run2'
stock=stock+xsa
stock@desc=xsa@desc
run2=stock

xsa.control <- FLXSA.control(tol = 1e-09, maxit = 99,  min.nse = 0.5,  fse  = 1.0,
                              rage = -1,   qage  = 7,   shk.n   = TRUE, shk.f = TRUE,
                              shk.yrs = 5, shk.ages= 4, window  = 100,  tsrange = 99,
                              tspower = 0)
xsa <- FLXSA(stock, tun, xsa.control)
xsa@desc <- 'Run3'
stock=stock+xsa
stock@desc=xsa@desc
run3=stock

# back to first run
stock <- run1
xsa <- xsa1

# for standard graphs
m(stock)
mat(stock)

stf0 <- stf(stock,nyears=1, wts.nyears=3, fbar.nyears=3)
srr <- FLSR(segreg)
R <- mean(stock.n(stock)[1,as.character(2003:2017)])
stock.n(stf0)[1,'2017'] <- R
ctrl <- projectControl(data.frame(year=2018,val=0,quantity='f'))
stf1  <- project(stf0, ctrl, srr)

stock.wt(stf1)[,'2018']
harvest(stf0)[,'2018']
catch.wt(stf1)[,'2018']

R
ssb(stf1)[,'2018']
#  142.217 kg

mean(ssb(stf1)[,c('2017','2018')])
mean(ssb(stf1)[,c('2014','2015','2016')])
mean(ssb(stf1)[,c('2017','2018')])/mean(ssb(stf1)[,c('2014','2015','2016')])
# raw values
#Index A (2016-2017) 	 76.7265
#Index B (2013-2015)	81.02705
#Index ratio (A/B)    0.9469245

mean(round(ssb(stf1)[,c('2017','2018')]))
mean(round(ssb(stf1)[,c('2014','2015','2016')]))
mean(round(ssb(stf1)[,c('2017','2018')]))/mean(round(ssb(stf1)[,c('2014','2015','2016')]))
# rounded values
#Index A (2016-2017) 	77
#Index B (2013-2015)	81
#Index ratio (A/B)    0.9506173


diagnostics(xsa)

# to get out values for standard graphs
ssb(run1)
rec(run1)
fbar(run1)
