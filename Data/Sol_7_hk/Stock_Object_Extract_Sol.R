library(FLCore)
library (plyr)

### Load stock object
sol=readFLStock('Data\\Sol_7_hk\\sol7jkID.txt')

units(sol)[1:17] <- as.list(c(rep(c("tonnes","thousands","kg"),4), "NA", "NA", "f", "NA", "NA"))

catch (sol)
landings(sol)
discards (sol)

### Extract catch data from the stock object
sol.data = as.data.frame(landings(sol))
sol.data = sol.data[ , c(2,7)]
colnames(sol.data) = c("year", "catch")



#### Cod Survey indices import
sol.indices=readFLIndices('Data\\Sol_7_hk\\sol7jkTU.txt')

### Indicator 1
ind1 = as.data.frame(sol.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
sol.data = merge(sol.data, ind1, all.x = T)


write.csv(sol.data, "Data\\Sol_7_hk\\Sol.Data.csv")
