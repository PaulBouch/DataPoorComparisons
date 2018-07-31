#### Plots to compare results
library (ggplot2)
library (gridExtra)
options(scipen = 999)

results_res_CMSY = read.csv("CMSY Sensitivity Trials\\Results\\results_resilience_CMSY.csv")
results_res_BSM = read.csv("CMSY Sensitivity Trials\\Results\\results_resilience_BSM.csv")


results_res_CMSY$Method = "CMSY"
results_res_BSM$Method = "BSM"

Res_test = rbind(results_res_CMSY, results_res_BSM)

fmsy.plot = ggplot (Res_test, aes(Method, fmsy, colour = Resilience))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = fmsy.low, ymax = fmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()

msy.plot = ggplot (Res_test, aes(Method, bmsy, colour = Resilience))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = bmsy.low, ymax = bmsy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()

bmsy.plot = ggplot (Res_test, aes(Method, msy, colour = Resilience))+
  geom_point(pch=1, size = 3, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = msy.low, ymax = msy.high), position = position_dodge(width = 0.5), size =0.3)+
  theme_bw()


grid.arrange(fmsy.plot, msy.plot, bmsy.plot, ncol= 1)


