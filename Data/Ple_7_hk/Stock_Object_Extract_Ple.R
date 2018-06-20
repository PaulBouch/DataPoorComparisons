library(FLCore)
library (plyr)

### Load stock object
ple=readFLStock('Data\\Ple_7_hk\\PLE7jkID.txt')
ple.indices=readFLIndices('Data\\Ple_7_hk\\PLE7jkTU.txt')


units(ple)[1:17] <- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

catch (ple)
landings(ple)
discards (ple)

### Extract catch data from the stock object
ple.data = as.data.frame(landings(ple))
ple.data = ple.data[ , c(2,7)]
colnames(ple.data) = c("year", "catch")



#### Cod Survey indices import
ple.indices=readFLIndices('Data\\Ple_7_hk\\PLE7jkTU.txt')

### Indicator 1
ind1 = as.data.frame(ple.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
ple.data = merge(ple.data, ind1, all.x = T)


write.csv(ple.data, "Data\\Ple_7_hk\\Ple.Data.csv")
