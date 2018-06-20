library(FLCore)
library (plyr)

### Load stock object
cod<-readFLStock("Data\\Cod_7_ek\\index2018disfinal.txt")
units(cod)[1:17] <- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

### Extract catch data from the stock object
cod.data = as.data.frame(landings(cod))
cod.data = cod.data[ , c(2,7)]
colnames(cod.data) = c("year", "catch")

#### Cod Survey indices import
cod.indices <- readFLIndices("Data\\Cod_7_ek\\SurveyIndex.txt")
names(cod.indices) <- c("FR-OTDEF","IR-FR COMBINED SURVEY")

### Indicator 1
ind1 = as.data.frame(cod.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
cod.data = merge(cod.data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(cod.indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
cod.data = merge(cod.data, ind2, all.x = T)

#write.csv(cod.data, "Data\\Cod_7_ek\\Cod.Data.csv")
