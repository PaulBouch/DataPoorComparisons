library(FLCore)
library (plyr)


Sol.7a.Data = read.csv("Data\\Sol_7_a\\sol7a catch pb.csv")

tun.indices <- readFLIndices("Data\\Sol_7_a\\tun.txt")

#tun.indices [[1]]

### Indicator 1
ind1 = as.data.frame(tun.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Sol.7a.Data = merge(Sol.7a.Data, ind1, all.x = T)



write.csv(Sol.7a.Data, "Data\\Sol_7_a\\Sol.7a.Data.csv")
