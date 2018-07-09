library(FLCore)
library (plyr)


Ple.7fg.Data = read.csv("Data\\Ple_7_fg\\ple7fg catch pb.csv")

tun.indices <- readFLIndices("Data\\Ple_7_fg\\P7FTUN3FULL.dat")

# tun.indices [[3]]

### Indicator 1
ind1 = as.data.frame(tun.indices [[3]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Ple.7fg.Data = merge(Ple.7fg.Data, ind1, all.x = T)


### Indicator 2
ind2 = as.data.frame(tun.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Ple.7fg.Data = merge(Ple.7fg.Data, ind2, all.x = T)


### Indicator 3
ind3 = as.data.frame(tun.indices [[2]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Ple.7fg.Data = merge(Ple.7fg.Data, ind3, all.x = T)


### Indicator 4
ind4 = as.data.frame(tun.indices [[4]]@index)
ind4 <- ddply(ind4, ~ year, summarize, index4 = sum(data))
Ple.7fg.Data = merge(Ple.7fg.Data, ind4, all.x = T)


write.csv(Ple.7fg.Data, "Data\\ple_7_fg\\Ple.7fg.Data.csv")
