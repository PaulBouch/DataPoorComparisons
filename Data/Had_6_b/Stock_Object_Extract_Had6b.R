library(FLCore)
library (plyr)


catch <-read.csv("Data\\Had_6_b\\Had6b catch pb.csv")
catch$catch = rowSums(catch[, 2:8])


Had.6b.Data = catch [ ,c(1,9)]


tun.indices <- readFLIndices("Data\\Had_6_b\\HAD6BTUN pb.txt")


### Indicator 1
ind1 = as.data.frame(tun.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Had.6b.Data = merge(Had.6b.Data, ind1, all.x = T)



write.csv(Had.6b.Data, "Data\\Had_6_b\\Had.6b.Data.csv")
