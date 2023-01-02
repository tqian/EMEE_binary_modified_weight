library(reshape)
library(kableExtra)
library(knitr)
library(ggplot2)

#######
# 1. different deltas
result_df_delta <- readRDS("result_table_delta(with sd).RDS")
efficiency_df_delta <- readRDS("efficiency_table_delta(with sd).RDS")

result_df_delta$var <- result_df_delta$sd^2
#*(result_df_delta$ss - 1)
result_df_delta_mod <- result_df_delta[result_df_delta$est == 'modified-EMEE',]
mod_var <- sapply(seq(1, nrow(result_df_delta_mod), 3), function(j) mean(result_df_delta_mod[j+(0:2), 'var']))
result_df_delta_ori <- result_df_delta[result_df_delta$est == 'EMEE',]
ori_var <- sapply(seq(1, nrow(result_df_delta_ori), 3), function(j) mean(result_df_delta_ori[j+(0:2), 'var']))

efficiency_means <- rowMeans(efficiency_df_delta[,2:4])
#efficiency_sds <- apply(result_df_delta[,2:4], 1, sd)

efficiency_df_delta$efficiency_means <- efficiency_means
efficiency_df_delta$mod_var <- mod_var
efficiency_df_delta$ori_var <- ori_var

colnames(efficiency_df_delta) <- c("Delta", "SS = 30", "SS = 50", "SS = 100", 
                                   "efficiency_means","mod_var", "ori_var")

#geom_point(aes(x = Delta, y = efficiency_means)) +
p <- ggplot(efficiency_df_delta) +
  geom_line(aes(x = Delta, y = efficiency_means, linetype = "Relative Efficiency")) +
  geom_point(aes(x = Delta, y = efficiency_means)) +
  geom_line(aes(x = Delta, y = 200 * mod_var, linetype = "Variance of pd-EMEE")) +
  geom_line(aes(x = Delta, y = 200 * ori_var, linetype = "Variance of EMEE")) +
  theme_bw() +
  labs(linetype="",x="Delta",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./200, name = "Variance of Estimators"))+ 
  theme(legend.position="bottom",legend.text = element_text(size = 6),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8)) 
  #scale_linetype_manual(values=c("solid", "dotted", "dashed"))
  #guides(color = guide_legend(override.aes = list(linetype = c(1, 0, 0),shape = c(16, NA,NA)) ) )

df_delta <- as.data.frame(rbind(cbind("Relative Efficiency", efficiency_df_delta$Delta, efficiency_df_delta$efficiency_means),
                  cbind("Variance of pd-EMEE", efficiency_df_delta$Delta, efficiency_df_delta$mod_var*200),
                  cbind("Variance of EMEE", efficiency_df_delta$Delta, efficiency_df_delta$ori_var*200)))
colnames(df_delta) <- c("grps", "Delta","Relative Efficiency")
df_delta$Delta <- as.numeric(df_delta$Delta)
df_delta$`Relative Efficiency` <- as.numeric(df_delta$`Relative Efficiency`)

p1 <- ggplot(data=df_delta, aes(x=Delta, y=`Relative Efficiency`, group=grps)) +
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
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Randomization Probability of Delta")


######
# 2. different probability of treatments
result_df_proba <- readRDS("result_table_proba(with sd).RDS")
efficiency_df_proba <- readRDS("efficiency_table_proba(with sd).RDS")

result_df_proba$var <- result_df_proba$sd^2
#*(result_df_delta$ss - 1)
result_df_proba_mod <- result_df_proba[result_df_proba$est == 'wcls_new',]
mod_var <- sapply(seq(1, nrow(result_df_proba_mod), 3), function(j) mean(result_df_proba_mod[j+(0:2), 'var']))
result_df_proba_ori <- result_df_proba[result_df_proba$est == 'wcls',]
ori_var <- sapply(seq(1, nrow(result_df_proba_ori), 3), function(j) mean(result_df_proba_ori[j+(0:2), 'var']))

efficiency_means <- rowMeans(efficiency_df_proba[,2:4])
#efficiency_sds <- apply(result_df_proba[,2:4], 1, sd)

efficiency_df_proba$efficiency_means <- efficiency_means
efficiency_df_proba$mod_var <- mod_var
efficiency_df_proba$ori_var <- ori_var

colnames(efficiency_df_proba) <- c("Prob_a(1)", "SS = 30", "SS = 50", "SS = 100", 
                                   "efficiency_means","mod_var", "ori_var")


q <- ggplot(efficiency_df_proba) +
  geom_line(aes(x = `Prob_a(1)`, y = efficiency_means, linetype = "Relative Efficiency")) +
  geom_line(aes(x = `Prob_a(1)`, y = 50 * mod_var, linetype = "Variance of pd-EMEE")) +
  geom_line(aes(x = `Prob_a(1)`, y = 50 * ori_var, linetype = "Variance of EMEE")) +
  geom_point(aes(x = `Prob_a(1)`, y = efficiency_means)) +
  theme_bw() +
  labs(linetype="",x="Randomization Probability of Treatment",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./50, name = "Variance of Estimators"))+ 
  theme(legend.position="bottom",legend.text = element_text(size = 6),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8)) + 
  scale_linetype_manual(values=c("solid", "dotted", "dashed"))

df_proba <- as.data.frame(rbind(cbind("Relative Efficiency", efficiency_df_proba$`Prob_a(1)`, efficiency_df_proba$efficiency_means),
                                cbind("Variance of pd-EMEE", efficiency_df_proba$`Prob_a(1)`, efficiency_df_proba$mod_var*50),
                                cbind("Variance of EMEE", efficiency_df_proba$`Prob_a(1)`, efficiency_df_proba$ori_var*50)))
colnames(df_proba) <- c("grps", "Randomization Probability of Treatment", "Relative Efficiency")
df_proba$`Randomization Probability of Treatment` <- as.numeric(df_proba$`Randomization Probability of Treatment`)
df_proba$`Relative Efficiency` <- as.numeric(df_proba$`Relative Efficiency`)

p2 <- ggplot(data=df_proba, aes(x=`Randomization Probability of Treatment`, y=`Relative Efficiency`, group=grps)) +
  geom_line(aes(linetype=grps))+
  geom_point(aes(size=grps))+
  scale_size_manual(values=c(1, 0, 0))+
  #labs(linetype="",x="Delta",y="Relative Efficiency") + 
  scale_y_continuous(sec.axis = sec_axis(~./50, name = "Variance of Estimators"))+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        legend.title = element_blank(),
        legend.position="bottom",
        legend.text = element_text(size = 6), 
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))+ 
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Randomization Probability of Treatment")

require(gridExtra)
pdf("figure2(combined).pdf")
grid.arrange(p1, p2, respect=TRUE)
dev.off()
