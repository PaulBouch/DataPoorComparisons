library(FLCore)
library (plyr)


Whg.NS.Data <-read.csv("Data\\Whg_NS\\Whg_NS Catch pb.csv")

indices <- readFLIndices("Data\\Whg_NS\\survey.dat")


### Indicator 1
ind1 = as.data.frame(indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Whg.NS.Data = merge(Whg.NS.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Whg.NS.Data = merge(Whg.NS.Data, ind2, all.x = T)


write.csv(Whg.NS.Data, "Data\\Whg_NS\\Whg.NS.Data.csv")
