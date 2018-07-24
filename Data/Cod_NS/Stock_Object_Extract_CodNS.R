library(FLCore)
library (plyr)


Cod.NS.Data <-read.csv("Data\\Cod_NS\\Cod NS Catch pb.csv")

cod.indices <- readFLIndices("Data\\Cod_NS\\survey.dat")

cod.indices[1]

### Indicator 1
ind1 = as.data.frame(cod.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Cod.NS.Data = merge(Cod.NS.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(cod.indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Cod.NS.Data = merge(Cod.NS.Data, ind2, all.x = T)



write.csv(Cod.NS.Data, "Data\\Cod_NS\\Cod.NS.Data.csv")
