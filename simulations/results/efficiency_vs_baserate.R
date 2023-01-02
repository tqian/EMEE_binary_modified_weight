library(reshape)
library(kableExtra)
library(knitr)
library(ggplot2)

#######
result_df_delta <- readRDS("result_table_delta(with sd).RDS")
efficiency_df_delta <- readRDS("efficiency_table_delta(with sd).RDS")

# base rate
base_rate_perDelta <- function(Delta) {
  C = (0.5^(0.5/Delta) + 1 + 0.5^(-0.5/Delta))
  prob_S_weight = c(0.5^(-0.5/Delta)/C, 1/C, 0.5^(0.5/Delta)/C)
  # E_S = 3*0.5^(1/Delta)/C
  
  base_rate_vals = c(1 - 0.5^(1 +0.5/Delta)* (3/C)^(Delta - 1),
                     1 - 0.5^(1 +0  /Delta)* (3/C)^(Delta - 1),
                     1 - 0.5^(1 -0.5/Delta)* (3/C)^(Delta - 1))
  e_base_rate = sum(prob_S_weight * base_rate_vals)
  return(e_base_rate)
}

base_rates = c()
for (Delta in 1:10) {
  base_rates = c(base_rates, base_rate_perDelta(Delta))
}

df_delta <- as.data.frame(rbind(cbind("Relative Efficiency", base_rates, efficiency_df_delta$efficiency_means),
                                cbind("Variance of pd-EMEE", base_rates, efficiency_df_delta$mod_var*200),
                                cbind("Variance of EMEE", base_rates, efficiency_df_delta$ori_var*200)))
colnames(df_delta) <- c("grps", "Base Rate","Relative Efficiency")
df_delta$`Base Rate` <- as.numeric(df_delta$`Base Rate`)
df_delta$`Relative Efficiency` <- as.numeric(df_delta$`Relative Efficiency`)

p1 <- ggplot(data=df_delta, aes(x=`Base Rate`, y=`Relative Efficiency`, group=grps)) +
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
  ggtitle("Relative Efficiency of pd-EMEE to EMEE over Different Values of Base Rates")

