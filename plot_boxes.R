setwd("C:/Users/test/OneDrive - Politecnico di Milano/Desktop/BHMixtures/BHMixtures")

load("simulation_A_results.RData")
load("simulation_B_results.RData")

library(ggplot2)

dfA <- rbind(
  data.frame(error = data_boxA[,1],
             par = data_boxA[,2],
             mix = "BH mixture"),
  data.frame(error = data_boxA_lm[,1],
             par = data_boxA_lm[,2],
             mix = "Linear mixture")
)

gA=ggplot(dfA, aes(x = factor(par), y = error, fill = mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  scale_fill_manual(values = c("lightblue","darkorchid4")) +
  ylim(0,1.4) +
  labs(x = "par_A",
       y = "L2 error - best vertex permutation",
       title = "Simulation A") +
  theme_minimal()

x11()
gA


library(ggplot2)

dfB <- rbind(
  data.frame(error = data_boxB[,1],
             par = data_boxB[,2],
             mix = "BH mixture"),
  data.frame(error = data_boxB_lm[,1],
             par = data_boxB_lm[,2],
             mix = "Linear mixture")
)

gB=ggplot(dfB, aes(x = factor(par), y = error, fill = mix)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  scale_fill_manual(values = c("lightblue","darkorchid4")) +
  ylim(0,1.4) +
  labs(x = "par_B",
       y = "L2 error - best vertex permutation",
       title = "Simulation B") +
  theme_minimal()

library(gridExtra)
x11()
grid.arrange(gA,gB,nrow=1)