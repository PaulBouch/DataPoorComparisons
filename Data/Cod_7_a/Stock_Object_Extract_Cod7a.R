library(FLCore)
library (plyr)


Cod.7a.Data <-read.csv("Data\\Cod_7_a\\COD7A Catch pb.csv")

cod.indices <- readFLIndices("Data\\Cod_7_a\\COD7ATUN.txt")


### Indicator 1
ind1 = as.data.frame(cod.indices [[2]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Cod.7a.Data = merge(Cod.7a.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(cod.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Cod.7a.Data = merge(Cod.7a.Data, ind2, all.x = T)

## Indicator 3
ind3 = as.data.frame(cod.indices [[3]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Cod.7a.Data = merge(Cod.7a.Data, ind3, all.x = T)

## Indicator 4
ind4 = as.data.frame(cod.indices [[4]]@index)
ind4 <- ddply(ind4, ~ year, summarize, index4 = sum(data))
Cod.7a.Data = merge(Cod.7a.Data, ind4, all.x = T)



write.csv(Cod.7a.Data, "Data\\Cod_7_a\\Cod.7a.Data.csv")
