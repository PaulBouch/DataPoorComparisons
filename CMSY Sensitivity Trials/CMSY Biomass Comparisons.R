#### Plots to compare results
library (ggplot2)
library (gridExtra)
options(scipen = 999)

results_biomass_CMSY = read.csv("CMSY Sensitivity Trials\\Results\\results_biomass_CMSY.csv")
results_biomass_BSM = read.csv("CMSY Sensitivity Trials\\Results\\results_biomass_BSM.csv")

results_biomass_CMSY$Method = "CMSY"
results_biomass_BSM$Method = "BSM"

Biomass_test = rbind(results_biomass_CMSY, results_biomass_BSM)

fmsy.plot = ggplot (Biomass_test, aes(Method, fmsy, colour = Profile))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()

msy.plot = ggplot (Biomass_test, aes(Method, bmsy, colour = Profile))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = bmsy.low, ymax = bmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()

bmsy.plot = ggplot (Biomass_test, aes(Method, msy, colour = Profile))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = msy.low, ymax = msy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()


grid.arrange(fmsy.plot, msy.plot, bmsy.plot, ncol= 1)
