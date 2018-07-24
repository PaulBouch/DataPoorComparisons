library(FLCore)
library (plyr)


Mac.Data <-read.csv("Data\\Mac_Wide\\Mac Catch pb.csv")

mac.indices <- readFLIndices("Data\\Mac_Wide\\survey.dat")


### Indicator 1
ind1 = as.data.frame(mac.indices [[2]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Mac.Data = merge(Mac.Data, ind1, all.x = T)

### Indicator 2
ind2 = as.data.frame(mac.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Mac.Data = merge(Mac.Data, ind2, all.x = T)

### Indicator 3
ind3 = as.data.frame(mac.indices [[3]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Mac.Data = merge(Mac.Data, ind3, all.x = T)


write.csv(Mac.Data, "Data\\Mac_Wide\\Mac.Data.csv")
