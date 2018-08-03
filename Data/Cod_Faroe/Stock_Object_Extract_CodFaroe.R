library(FLCore)
library (plyr)


Cod.Faroe.Data <-read.csv("Data\\Cod_Faroe\\CODFaroe Catch pb.csv")

cod.indices <- readFLIndices("Data\\Cod_Faroe\\survey.dat")

cod.indices[1]

### Indicator 1
ind1 = as.data.frame(cod.indices [[2]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Cod.Faroe.Data = merge(Cod.Faroe.Data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(cod.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Cod.Faroe.Data = merge(Cod.Faroe.Data, ind2, all.x = T)



write.csv(Cod.Faroe.Data, "Data\\Cod_Faroe\\Cod.Faroe.Data.csv")
