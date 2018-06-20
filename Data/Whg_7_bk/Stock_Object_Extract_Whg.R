library(FLCore)
library (plyr)

### Load stock object
whg <-stock
units(whg)[1:17] <- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

catch (whg)

### Extract catch data from the stock object
whg.data = as.data.frame(catch(whg))
whg.data = whg.data[ , c(2,7)]
colnames(whg.data) = c("year", "catch")



#### Cod Survey indices import
whg.indices <- tun

### Indicator 1
ind1 = as.data.frame(whg.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
whg.data = merge(whg.data, ind1, all.x = T)

# ## Indicator 2
# ind2 = as.data.frame(whg.indices [[2]]@index)
# ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
# whg.data = merge(whg.data, ind2, all.x = T)

write.csv(whg.data, "Data\\Whg_7_bk\\Whg.Data.csv")
