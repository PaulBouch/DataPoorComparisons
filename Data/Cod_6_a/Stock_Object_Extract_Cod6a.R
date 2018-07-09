library(FLCore)
library (plyr)

### Load stock object
discards <-read.csv("Data\\Cod_6_a\\cod via discards tonnes pb.csv")
discards$discards = rowSums(discards[, 2:8])

landings = read.csv("Data\\Cod_6_a\\cod via landings tonnes pb.csv")
landings$landings = rowSums(landings[, 2:8])

Cod.6a.Data = cbind (landings, discards)
Cod.6a.Data$catch = Cod.6a.Data$landings + Cod.6a.Data$discards

Cod.6a.Data = Cod.6a.Data [ ,c(1,19)]


cod.indices <- readFLIndices("Data\\Cod_6_a\\COD6A.surveys.dat")


### Indicator 1
ind1 = as.data.frame(cod.indices [[3]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(cod.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind2, all.x = T)

## Indicator 3
ind3 = as.data.frame(cod.indices [[2]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind3, all.x = T)

## Indicator 4
ind4 = as.data.frame(cod.indices [[4]]@index)
ind4 <- ddply(ind4, ~ year, summarize, index4 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind4, all.x = T)

## Indicator 5
ind5 = as.data.frame(cod.indices [[5]]@index)
ind5 <- ddply(ind5, ~ year, summarize, index5 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind5, all.x = T)

## Indicator 6
ind6 = as.data.frame(cod.indices [[6]]@index)
ind6 <- ddply(ind6, ~ year, summarize, index6 = sum(data))
Cod.6a.Data = merge(Cod.6a.Data, ind6, all.x = T)

write.csv(Cod.6a.Data, "Data\\Cod_6_a\\Cod.6a.Data.csv")
