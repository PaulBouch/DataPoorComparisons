library(FLCore)
library (plyr)


Sol.2024.Data <-read.csv("Data\\Sol_2024\\Sol2024 Catch pb.csv")

indices <- readFLIndices("Data\\Sol_2024\\survey pb.txt")


### Indicator 1
ind1 = as.data.frame(indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Sol.2024.Data = merge(Sol.2024.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Sol.2024.Data = merge(Sol.2024.Data, ind2, all.x = T)

## Indicator 3
ind3 = as.data.frame(indices [[3]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Sol.2024.Data = merge(Sol.2024.Data, ind3, all.x = T)


write.csv(Sol.2024.Data, "Data\\Sol_2024\\Sol.2024.Data.csv")
