rm(list = ls())

library(reshape)
library(tidyverse)

# relative efficiency vs. delta -------------------------------------------

## load simulation results
efficiency_table_delta = readRDS("simulations/results/efficiency_table_delta.RDS")
result_table_delta = readRDS("simulations/results/result_table_delta.RDS")

## construct date frame for making plot
efficiency_means <- rowMeans(efficiency_table_delta[,2:4])
efficiency_table_delta$efficiency_means <- efficiency_means
colnames(efficiency_table_delta) <- c("gamma", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")

var <- result_table_delta$sd^2 
result_table_delta$var <- var
mod_df <- result_table_delta[result_table_delta$est == "modified-EMEE",-2]
mod_avg = mod_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
mod_var <- mod_avg$var

ori_df <- result_table_delta[result_table_delta$est == "EMEE",-2]
ori_avg = ori_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
ori_var <- ori_avg$var

df_delta <- as.data.frame(rbind(cbind("Relative Efficiency", 1:10, efficiency_table_delta$efficiency_means),
                                cbind("Variance of pd-EMEE", 1:10, mod_var*200),
                                cbind("Variance of EMEE", 1:10, ori_var*200)))
colnames(df_delta) <- c("grps", "Delta","Relative Efficiency")
df_delta$`Delta` <- as.numeric(df_delta$`Delta`)
df_delta$`Relative Efficiency` <- as.numeric(df_delta$`Relative Efficiency`)

## make plot
p1 <- ggplot(data=df_delta, aes(x=`Delta`, y=`Relative Efficiency`, group=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_size_manual(values=c(1, 0, 0))+
  #labs(linetype="",x="Delta",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Delta")
p1

# relative efficiency vs. p (rand. prob.) ---------------------------------

efficiency_table_proba = readRDS("simulations/results/efficiency_table_proba.RDS")
result_table_proba = readRDS("simulations/results/result_table_proba.RDS")

## construct date frame for making plot
efficiency_means <- rowMeans(efficiency_table_proba[,2:4])
efficiency_table_proba$efficiency_means <- efficiency_means
colnames(efficiency_table_proba) <- c("proba", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")

var <- result_table_proba$sd^2 
result_table_proba$var <- var
mod_df <- result_table_proba[result_table_proba$est == "modified-EMEE",-2]
mod_avg = mod_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
mod_var <- mod_avg$var

ori_df <- result_table_proba[result_table_proba$est == "EMEE",-2]
ori_avg = ori_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
ori_var <- ori_avg$var

df_proba <- as.data.frame(rbind(cbind("Relative Efficiency", efficiency_table_proba$proba, efficiency_table_proba$efficiency_means),
                                cbind("Variance of pd-EMEE", efficiency_table_proba$proba, mod_var*50),
                                cbind("Variance of EMEE", efficiency_table_proba$proba, ori_var*50)))
colnames(df_proba) <- c("grps", "Randomization Probability of Treatment","Relative Efficiency")
df_proba$`Randomization Probability of Treatment` <- as.numeric(df_proba$`Randomization Probability of Treatment`)
df_proba$`Relative Efficiency` <- as.numeric(df_proba$`Relative Efficiency`)

## make plot
p2 <- ggplot(data=df_proba, aes(x=`Randomization Probability of Treatment`, y=`Relative Efficiency`, group=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_size_manual(values=c(1, 0, 0))+
  scale_y_continuous(sec.axis = sec_axis(~./50, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Randomization Probability of Treatment")
p2


# relative efficiency vs. K (reference regime) ----------------------------

## load simulation results
efficiency_table_K = readRDS("simulations/results/efficiency_table_K.RDS")
result_table_K = readRDS("simulations/results/result_table_K.RDS")

## construct date frame for making plot
efficiency_means <- rowMeans(efficiency_table_K[,2:4])
efficiency_table_K$efficiency_means <- efficiency_means
colnames(efficiency_table_K) <- c("K", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")

var <- result_table_K$sd^2 
result_table_K$var <- var
mod_df <- result_table_K[result_table_K$est == "modified-EMEE",-2]
mod_avg = mod_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
mod_var <- mod_avg$var

ori_df <- result_table_K[result_table_K$est == "EMEE",-2]
ori_avg = ori_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
ori_var <- ori_avg$var

df_K <- as.data.frame(rbind(cbind("Relative Efficiency", 1:10, efficiency_table_K$efficiency_means),
                                cbind("Variance of pd-EMEE", 1:10, mod_var*200),
                                cbind("Variance of EMEE", 1:10, ori_var*200)))
colnames(df_K) <- c("grps", "K","Relative Efficiency")
df_K$`K` <- as.numeric(df_K$`K`)
df_K$`Relative Efficiency` <- as.numeric(df_K$`Relative Efficiency`)

## make plot
p3 <- ggplot(data=df_K, aes(x=`K`, y=`Relative Efficiency`, group=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_size_manual(values=c(1, 0, 0))+
  #labs(linetype="",x="K",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of K")
p3


