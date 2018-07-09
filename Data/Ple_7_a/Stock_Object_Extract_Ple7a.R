library(FLCore)
library (plyr)


Ple.7a.Data = read.csv("Data\\Ple_7_a\\ple7a catch pb.csv")

tun.indices <- readFLIndices("Data\\Ple_7_a\\surveys pb.txt")


### Indicator 1
ind1 = as.data.frame(tun.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Ple.7a.Data = merge(Ple.7a.Data, ind1, all.x = T)



write.csv(Ple.7a.Data, "Data\\ple_7_a\\Ple.7a.Data.csv")
