library(FLCore)
library (plyr)


Had.NS.Data <-read.csv("Data\\Had_NS\\Had_NS Catch pb.csv")

indices <- readFLIndices("Data\\Had_NS\\survey.dat")


### Indicator 1
ind1 = as.data.frame(indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Had.NS.Data = merge(Had.NS.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Had.NS.Data = merge(Had.NS.Data, ind2, all.x = T)


write.csv(Had.NS.Data, "Data\\Had_NS\\Had.NS.Data.csv")
