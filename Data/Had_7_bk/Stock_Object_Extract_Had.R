library(FLCore)
library (plyr)

### Load stock object
had <-stock
units(had)[1:17] <- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

catch (had)

### Extract catch data from the stock object
had.data = as.data.frame(catch(had))
had.data = had.data[ , c(2,7)]
colnames(had.data) = c("year", "catch")



#### Cod Survey indices import
had.indices <- tun

### Indicator 1
ind1 = as.data.frame(had.indices [[2]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
had.data = merge(had.data, ind1, all.x = T)

## Indicator 2
ind2 = as.data.frame(had.indices [[1]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
had.data = merge(had.data, ind2, all.x = T)

write.csv(had.data, "Data\\Had_7_bk\\Had.Data.csv")
