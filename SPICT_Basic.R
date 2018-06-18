require (plyr)
require(ggplot2)
library(spict)
library(reshape2)
require (gridExtra)

inp <- list()
inp$obsC <- cod.data$catch
inp$timeC <- cod.data$year
inp$obsI[[1]] <- cod.data$index1 [30:47]
inp$obsI[[2]] <- cod.data$index2 [33:47]
inp$timeI <- list(cod.data$year[30:47], cod.data$year[33:47])

### Use only one abundance index
# inp$obsI <- cod.data$index1 [30:47]
# inp$timeI <- cod.data$year[30:47]

### provide r prior
inp$priors$logr <- c(-0.76, 0.36, 1)
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

# # Just plot the absolute biomass and relative F
# png (file=paste(stock, "relative.png"))
# par(mfrow = c(2, 1), mar = c(4, 4, 1, 1))
# plotspict.bbmsy(res, qlegend = FALSE, stamp = "")
# plotspict.ffmsy(res, qlegend = FALSE, stamp = "")
# dev.off()

# ####   Produce the reference value table
# tab1 <- sumspict.srefpoints(res);
# tab1 = as.data.frame(tab1)
# tab1$rel.diff.Drp = NULL
# tab2 <- sumspict.states(res);
#   
# tab = rbind(tab1, tab2)
# tab = round(tab, digits = 3)
#   
# png (file=paste(stock, "refs.png"), width = 400, height = 250)
# grid.table(tab)
# dev.off()
#   
# output = data.frame(stock,
#                       tab[2,1], tab[2,2], tab[2,3],
#                       tab[1,1], tab[1,2], tab[1,3],
#                       tab[6,1], tab[6,2], tab[6,3],
#                       tab[7,1], tab[7,2], tab[7,3])
#   
# write.table(output, file=outfile, append = T, sep = ",", 
#               dec = ".", row.names = FALSE, col.names = FALSE)            

  


