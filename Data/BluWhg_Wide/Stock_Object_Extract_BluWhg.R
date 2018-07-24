library(FLCore)
library (plyr)


BluWhg.Data <-read.csv("Data\\BluWhg_Wide\\BluWhg Catch pb.csv")

indices <- readFLIndices("Data\\BluWhg_Wide\\survey.dat")


### Indicator 1
ind1 = as.data.frame(indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
BluWhg.Data = merge(BluWhg.Data, ind1, all.x = T)


write.csv(BluWhg.Data, "Data\\BluWhg_Wide\\BluWhg.Data.csv")
