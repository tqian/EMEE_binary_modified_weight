rm(list = ls())

library(reshape)
library(tidyverse)
library(latex2exp)
library(gridExtra)

# relative efficiency vs. delta -------------------------------------------

## load simulation results
efficiency_table_delta = as.data.frame(readRDS("pdEMEE_efficiency_Delta.RDS"))
efficiency_table_delta2 = as.data.frame(readRDS("pdEMEE2_efficiency_Delta.RDS"))
result_table_delta = as.data.frame(readRDS("result_table_delta_wEMEE2.RDS"))

## construct date frame for making plot
efficiency_means <- rowMeans(efficiency_table_delta[,2:4])
efficiency_means2 <- rowMeans(efficiency_table_delta2[,2:4])
efficiency_table_delta$efficiency_means <- efficiency_means
efficiency_table_delta2$efficiency_means <- efficiency_means2
colnames(efficiency_table_delta) <- c("gamma", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")
colnames(efficiency_table_delta2) <- c("gamma", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")

var <- result_table_delta$sd^2 
result_table_delta$var <- var

improved_df <- result_table_delta[result_table_delta$est == "improved-EMEE",-2]
improved_avg = improved_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
improved_var <- improved_avg$var

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

df_delta <- as.data.frame(rbind(cbind("Relative Efficiency of pd-EMEE over EMEE", 1:10, efficiency_table_delta$efficiency_means),
                                cbind("Relative Efficiency of pd-EMEE2 over EMEE", 1:10, efficiency_table_delta2$efficiency_means),
                                cbind("Variance of pd-EMEE2", 1:10, improved_var*200),
                                cbind("Variance of pd-EMEE", 1:10, mod_var*200),
                                cbind("Variance of EMEE", 1:10, ori_var*200)))
colnames(df_delta) <- c("grps", "Delta","Relative Efficiency")
df_delta$`Delta` <- as.numeric(df_delta$`Delta`)
df_delta$`Relative Efficiency` <- as.numeric(df_delta$`Relative Efficiency`)

## make plot
p1 <- ggplot(data=df_delta, aes(x=`Delta`, y=`Relative Efficiency`, color=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_linetype_manual(values=c("solid", "solid", "dashed", "twodash", "dotted"))+
  scale_size_manual(values=c(1, 1, 0, 0, 0))+
  scale_color_manual(values = c("red", "blue", "black", "black", "black")) +
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  xlab(TeX("$Delta$ (proximal outcome window length)"))
  #ggtitle("Relative Efficiency of pd-EMEE and pd-EMEE2 to EMEE over Different Values of Delta")
p1

# relative efficiency vs. p (rand. prob.) ---------------------------------

efficiency_table_proba = as.data.frame(readRDS("pdEMEE_efficiency_proba.RDS"))
efficiency_table_proba2 = as.data.frame(readRDS("pdEMEE2_efficiency_proba.RDS"))
result_table_proba = as.data.frame(readRDS("result_table_proba_wEMEE2.RDS"))

## construct date frame for making plot
efficiency_means <- rowMeans(efficiency_table_proba[,2:4])
efficiency_means2 <- rowMeans(efficiency_table_proba2[,2:4])
efficiency_table_proba$efficiency_means <- efficiency_means
efficiency_table_proba2$efficiency_means <- efficiency_means2
colnames(efficiency_table_proba) <- c("gamma", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")
colnames(efficiency_table_proba2) <- c("gamma", "SS = 30", "SS = 50", "SS = 100", "efficiency_means")

var <- result_table_proba$sd^2 
result_table_proba$var <- var

improved_df <- result_table_proba[result_table_proba$est == "wcls_improved",-2]
improved_avg = improved_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
improved_var <- improved_avg$var

mod_df <- result_table_proba[result_table_proba$est == "wcls_new",-2]
mod_avg = mod_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
mod_var <- mod_avg$var

ori_df <- result_table_proba[result_table_proba$est == "wcls",-2]
ori_avg = ori_df %>%
  group_by(group = gl(n()/3, 3)) %>%
  summarise_at(-1, mean, na.rm = TRUE)
ori_var <- ori_avg$var

df_proba <- as.data.frame(rbind(cbind("Relative Efficiency of pd-EMEE over EMEE", efficiency_table_proba$gamma, efficiency_table_proba$efficiency_means),
                                cbind("Relative Efficiency of pd-EMEE2 over EMEE", efficiency_table_proba2$gamma, efficiency_table_proba2$efficiency_means),
                                cbind("Variance of pd-EMEE2", efficiency_table_proba$gamma, improved_var*50),
                                cbind("Variance of pd-EMEE", efficiency_table_proba$gamma, mod_var*50),
                                cbind("Variance of EMEE", efficiency_table_proba$gamma, ori_var*50)))
colnames(df_proba) <- c("grps", "Randomization Probability of Treatment","Relative Efficiency")
df_proba$`Randomization Probability of Treatment` <- as.numeric(df_proba$`Randomization Probability of Treatment`)
df_proba$`Relative Efficiency` <- as.numeric(df_proba$`Relative Efficiency`)

## make plot
p2 <- ggplot(data=df_proba, aes(x=`Randomization Probability of Treatment`, y=`Relative Efficiency`, color=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_linetype_manual(values=c("solid", "solid", "dashed", "twodash", "dotted"))+
  scale_size_manual(values=c(1, 1, 0, 0, 0))+
  scale_color_manual(values = c("red", "blue", "black", "black", "black")) +
  scale_y_continuous(sec.axis = sec_axis(~./50, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  xlab(TeX("$p_a$ (randomization prob.)")) 
  #ggtitle("Relative Efficiency of pd-EMEE and pd-EMEE2 to EMEE over Different Values of Randomization Probability of Treatment")
p2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1)
p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
grid.arrange(p1, p2, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.2))

