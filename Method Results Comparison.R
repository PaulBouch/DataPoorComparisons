require (plyr)
require(ggplot2)
require(reshape)
require (reshape2)
require (gridExtra)

options(scipen = 999)

### Read in data
Results <- read.csv("Results\\Results_All_030718.csv")

### Or read in all results folders at once
# file_names <- dir("Results") 
# setwd("Results")
# Results = do.call(rbind,lapply(file_names ,read.csv))
# setwd ("..")


ICES <- read.csv("Data\\ICES Refs 20_06_18.csv")

Results = Results[ grep("Convergence: 1", Results$Converge, invert = TRUE) , ]

Res2 = merge(Results, ICES, by = "Stock")

Res2$Name.N = ifelse(Res2$FixN2 == "Y", "N2", "")
Res2$Name.r = ifelse (Res2$Use.r == "Y", "r", "")
Res2$Name.k = ifelse (Res2$Use.k == "Y", "k", "")
Res2$Name.q = ifelse (Res2$Use.q == "Y", "q", "")
Res2$Name.Multi = ifelse (Res2$Index == "1", "Single", "Multi")

Res2$FullMethod = ifelse (Res2$Method == "SPICT", paste(Res2$Method, Res2$Name.N, Res2$Name.r, 
                                                        Res2$Name.k, Res2$Name.q, Res2$Name.Multi, sep = ""), 
                          ifelse (Res2$Method == "CMSY", "CMSY", "BSM"))



Res2$FMSY_ARE = 100* abs((Res2$fmsy - Res2$ICES_Fmsy)/Res2$ICES_Fmsy) 
Res2$F_ARE = 100* abs((Res2$f.end - Res2$ICES_F_end)/Res2$ICES_F_end) 
Res2$B_ARE = 100* abs((Res2$b.end - Res2$ICES_B_end)/Res2$ICES_B_end) 
Res2$F_FMSY_ARE = 100* abs(Res2$ffmsy - (Res2$ICES_F_end/Res2$ICES_Fmsy)/(Res2$ICES_F_end/Res2$ICES_Fmsy)) 

ARE_Summary = ddply(Res2, "FullMethod", summarise, 
                         Fmsy_mean = mean(FMSY_ARE) ,
                         Fmsy_median = median (FMSY_ARE),
                         F_Fmsy_mean = mean(F_FMSY_ARE, na.rm=TRUE),
                         F_Fmsy_median = median(F_FMSY_ARE, na.rm=TRUE),
                         B_mean = mean(B_ARE, na.rm=TRUE),
                         B_median = median(B_ARE, na.rm=TRUE))


##### Cod 7ek


ggplot(Res2[Res2$Stock == "Cod_7_ek", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
#  coord_cartesian(ylim=c(0, 1)) +
#  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Cod_7_ek", "ICES_Fmsy" ]), color = "red")

ggplot(Res2[Res2$Stock == "Had_7_bk", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
  coord_cartesian(ylim=c(0, 1)) +
  #  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Had_7_bk", "ICES_Fmsy" ]), color = "red")

ggplot(Res2[Res2$Stock == "Whg_7_bk", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
  #  coord_cartesian(ylim=c(0, 1)) +
  #  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Whg_7_bk", "ICES_Fmsy" ]), color = "red")

ggplot(Res2[Res2$Stock == "Whg_7_a", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
  coord_cartesian(ylim=c(0, 1)) +
  #  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Whg_7_a", "ICES_Fmsy" ]), color = "red")

ggplot(Res2[Res2$Stock == "Sol_7_hk", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
  #  coord_cartesian(ylim=c(0, 1)) +
  #  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Sol_7_hk", "ICES_Fmsy" ]), color = "red")

ggplot(Res2[Res2$Stock == "Ple_7_hk", ], 
       aes(FullMethod, fmsy))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  theme_bw() +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme(legend.title = element_blank(), legend.position = "bottom")+
  coord_cartesian(ylim=c(0, 1)) +
  #  scale_color_manual(values=c("black", "red", "cyan3", "orange")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_hline(yintercept =(Res2[Res2$Stock == "Ple_7_hk", "ICES_Fmsy" ]), color = "red")

