# Plot --------------------------------------------------------------------

library(tidyverse)
library(patchwork) # to arrange plots
library(latex2exp)
library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

###########
# data engineering
dta <- readRDS("FINAL Dataset_A.rds")
dta <- as.data.frame(dta)
dta$days_since_download <- dta$days_since_download - 1
dta$after_9pm_day_before_before_8pm <- dta$after_9pm_day_before* dta$before_8pm

head(dta)

moderators <- c("days_since_download", "treatment_day_before", 
                    "after_9pm_day_before_before_8pm", "NULL")

plot_est_collected <- c()
plot_re_collected <- c()

for (moderator in moderators) {
  if (moderator == "NULL") {
    plot_title <- "Moderator: None (Marginal)"
  } else if (moderator == "days_since_download") {
    plot_title <- "Moderator: Days since download"
  } else if (moderator == "treatment_day_before") {
    plot_title <- "Moderator: Whether user receive treatment \
    the day before"
  } else if (moderator == "after_9pm_day_before_before_8pm") {
    plot_title <- "Moderator: whether the user opened the app \
    between 8 PM and 9 PM the day before"
  } 
  xlabel <- TeX('$\\Delta$ (proximal outcome window length)')
  
  if (moderator == "NULL"){
    ylab_label_p1 <- TeX(r'($\hat{\beta}_0$ ($\pm$1.96 SE))')
    ylab_label_p2 <- TeX(r'(RE($\hat{\beta}_0$))')
    
    beta_collected <- data.frame(readRDS(paste0("analysis_result/moderator=", moderator, ".RDS")))
    colnames(beta_collected) <- c("Delta", "est", "se", "rci", "lci", "method", "estimand", "re")
    rownames(beta_collected) <- NULL
    beta_collected[, c(1:5, 8)] <- round(apply(beta_collected[, c(1:5, 8)], 2, as.numeric), 3)
    
    beta_collected$method <- ifelse(beta_collected$method == "IPW", "EMEE", "pd-EMEE")
    
    p1 <- beta_collected %>%
      ggplot(aes(x = factor(Delta), y = est, color = method)) + 
      geom_point(position = position_dodge(width = 0.6),
                 size = 4) + 
      geom_linerange(aes(ymin = lci, ymax = rci),
                     position = position_dodge2(width = 0.6),
                     linewidth = 1.2) +
      theme_bw() +
      ggtitle(plot_title) +
      scale_color_manual(values = mypalette) +
      xlab(xlabel) +
      ylab(ylab_label_p1) +
      theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 18),
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))
    
    p2 <- beta_collected %>%
      filter(method == "pd-EMEE") %>%
      ggplot(aes(x = factor(Delta), y = re, group = 1)) +
      geom_point(size = 4, color = "black") +
      geom_hline(yintercept = 1, linetype = 2) +
      # coord_cartesian(ylim = c(0.95, 1.5)) +
      theme_bw() +
      xlab(NULL) +
      ylab(ylab_label_p2) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 18),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 18))
    
    plot_est_collected <- c(plot_est_collected, list(p1))
    plot_re_collected <- c(plot_re_collected, list(p2))
    
  } else {
    print(plot_title)
    
    for (this_estimand in c("beta0", "beta1")) {
      
      beta_collected <- data.frame(readRDS(paste0("analysis_result/moderator=", moderator, ".RDS")))
      colnames(beta_collected) <- c("Delta", "est", "se", "rci", "lci", "method", "estimand", "re")
      rownames(beta_collected) <- NULL
      beta_collected[, c(1:5, 8)] <- round(apply(beta_collected[, c(1:5, 8)], 2, as.numeric), 3)
      
      beta_collected$method <- ifelse(beta_collected$method == "IPW", "EMEE", "pd-EMEE")
      beta_collected <- filter(beta_collected, estimand == this_estimand)
      
      if (this_estimand == "beta0") {
        ylab_label_p1 <- TeX(r'($\hat{\beta}_0$ ($\pm$1.96 SE))')
        ylab_label_p2 <- TeX(r'(RE($\hat{\beta}_0$))')
      } else {
        ylab_label_p1 <- TeX(r'($\hat{\beta}_1$ ($\pm$1.96 SE))')
        ylab_label_p2 <- TeX(r'(RE($\hat{\beta}_1$))')
      }
      
      
      p1 <- beta_collected %>%
        ggplot(aes(x = factor(Delta), y = est, color = method)) + 
        geom_point(position = position_dodge(width = 0.6),
                   size = 4) + 
        geom_linerange(aes(ymin = lci, ymax = rci),
                       position = position_dodge2(width = 0.6),
                       linewidth = 1.2) +
        theme_bw() +
        scale_color_manual(values = mypalette) +
        xlab(xlabel) +
        ylab(ylab_label_p1) +
        theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 18),
              plot.title = element_text(size = 20, hjust = 0.5),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 18))
      if (this_estimand == "beta0") {
        p1 <- p1 + ggtitle(plot_title)
      }
      
      p2 <- beta_collected %>%
        filter(method == "pd-EMEE") %>%
        ggplot(aes(x = factor(Delta), y = re, group = 1)) +
        geom_point(size = 4, color = "black") +
        geom_hline(yintercept = 1, linetype = 2) +
        # coord_cartesian(ylim = c(0.95, 1.5)) +
        theme_bw() +
        xlab(NULL) +
        ylab(ylab_label_p2) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 18),
              legend.title = element_text(size = 20),
              legend.text = element_text(size = 18))
      
      plot_est_collected <- c(plot_est_collected, list(p1))
      plot_re_collected <- c(plot_re_collected, list(p2))
    }
  }
}

pdf(file = "Drinkless_analysis_plot.pdf", width = 14, height = 16)
plot_est_collected[[1]] + plot_re_collected[[1]] + plot_est_collected[[2]] + plot_re_collected[[2]] +
  plot_est_collected[[3]] + plot_re_collected[[3]] + plot_est_collected[[4]] + plot_re_collected[[4]] +
  plot_est_collected[[5]] + plot_re_collected[[5]] + plot_est_collected[[6]] + plot_re_collected[[6]] +
  plot_est_collected[[7]] + plot_re_collected[[7]] +  
  guide_area() +
  plot_layout(nrow = 8, ncol = 2, byrow = FALSE, height = rep(c(1, 2), 9), guides = "collect")
dev.off()
