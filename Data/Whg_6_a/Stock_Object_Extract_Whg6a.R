library(FLCore)
library (plyr)


Whg.6a.Data = read.csv("Data\\Whg_6_a\\whg6a catch pb.csv")

tun.indices = readFLIndices("Data\\Whg_6_a\\tun combine pb.txt")

#tun.indices[[5]]


### Indicator 1
ind1 = as.data.frame(tun.indices [[1]]@index)
ind1 <- ddply(ind1, ~ year, summarize, index1 = sum(data))
Whg.6a.Data = merge(Whg.6a.Data, ind1, all.x = T)


### Indicator 2
ind2 = as.data.frame(tun.indices [[2]]@index)
ind2 <- ddply(ind2, ~ year, summarize, index2 = sum(data))
Whg.6a.Data = merge(Whg.6a.Data, ind2, all.x = T)


### Indicator 3
ind3 = as.data.frame(tun.indices [[3]]@index)
ind3 <- ddply(ind3, ~ year, summarize, index3 = sum(data))
Whg.6a.Data = merge(Whg.6a.Data, ind3, all.x = T)


# ### Indicator 4
# ind4 = as.data.frame(tun.indices [[4]]@index)
# ind4 <- ddply(ind4, ~ year, summarize, index4 = sum(data))
# Whg.6a.Data = merge(Whg.6a.Data, ind4, all.x = T)
# 
# ### Indicator 4
# ind5 = as.data.frame(tun.indices [[5]]@index)
# ind5 <- ddply(ind5, ~ year, summarize, index5 = sum(data))
# Whg.6a.Data = merge(Whg.6a.Data, ind5, all.x = T)


write.csv(Whg.6a.Data, "Data\\Whg_6_a\\Whg.6a.Data.Combine.csv")
