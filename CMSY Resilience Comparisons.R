#### Plots to compare results
library (ggplot2)
library (gridExtra)
options(scipen = 999)

results_resilience$Method = "CMSY"
results_resilience_BSM$Method = "BSM"

Res_test = rbind(results_resilience, results_resilience_BSM)

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


