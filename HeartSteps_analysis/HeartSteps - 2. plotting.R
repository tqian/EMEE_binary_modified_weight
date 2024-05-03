# Tianchen Qian
# 2023.12.01

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Plot --------------------------------------------------------------------

library(tidyverse)
library(patchwork) # to arrange plots
library(latex2exp)
library(RColorBrewer)
mypalette <- brewer.pal(n = 5, "Set1")

plot_est_collected <- c()
plot_re_collected <- c()

moderators <- c("decision.index.nogap.new", "at_home_or_work", "is_weekday", "NULL")

for (moderator in moderators) {
    
    if (moderator == "NULL") {
        plot_title <- "Moderator: None"
    } else if (moderator == "decision.index.nogap.new") {
        plot_title <- "Moderator: Decision point index"
    } else if (moderator == "at_home_or_work") {
        plot_title <- "Moderator: At home/work "
    } else if (moderator == "is_weekday") {
        plot_title <- "Moderator: Is weekday"
    }
    
    if (moderator == "NULL") {
        ylab_label_p1 <- TeX(r'($\hat{\beta}_0$ ($\pm$1.96 SE))')
        ylab_label_p2 <- TeX(r'(RE($\hat{\beta}_0$))')
        
        beta_collected <- readRDS(paste0("analysis_result/moderator=", moderator, ".RDS"))
        
        beta_collected$method <- ifelse(beta_collected$method == "IPW", "EMEE", "pd-EMEE")
        
        p1 <- beta_collected %>%
            ggplot(aes(x = factor(threshold/1000), y = est, color = method)) + 
            geom_point(position = position_dodge(width = 0.6),
                       size = 4) + 
            geom_linerange(aes(ymin = lci, ymax = rci),
                           position = position_dodge2(width = 0.6),
                           linewidth = 1.2) +
            theme_bw() +
            ggtitle(plot_title) +
            scale_color_manual(values = mypalette) +
            xlab("Outcome threshold (x1000 steps)") +
            ylab(ylab_label_p1) +
            theme(axis.text.x = element_text(size = 14, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 12),
                  axis.title = element_text(size = 18),
                  plot.title = element_text(size = 20, hjust = 0.5),
                  legend.title = element_text(size = 20),
                  legend.text = element_text(size = 18))
        
        p2 <- beta_collected %>%
            filter(method == "pd-EMEE") %>%
            ggplot(aes(x = factor(threshold/1000), y = re, group = 1)) +
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
        for (this_estimand in c("beta0", "beta1")) {
            
            beta_collected <- readRDS(paste0("analysis_result/moderator=", moderator, ".RDS"))
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
                ggplot(aes(x = factor(threshold/1000), y = est, color = method)) + 
                geom_point(position = position_dodge(width = 0.6),
                           size = 4) + 
                geom_linerange(aes(ymin = lci, ymax = rci),
                               position = position_dodge2(width = 0.6),
                               linewidth = 1.2) +
                theme_bw() +
                scale_color_manual(values = mypalette) +
                xlab("Outcome threshold (x1000 steps)") +
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
                ggplot(aes(x = factor(threshold/1000), y = re, group = 1)) +
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

pdf(file = "HeartSteps_analysis_plot.pdf", width = 21, height = 8)
plot_est_collected[[1]] + plot_re_collected[[1]] + plot_est_collected[[2]] + plot_re_collected[[2]] +
    plot_est_collected[[3]] + plot_re_collected[[3]] + plot_est_collected[[4]] + plot_re_collected[[4]] +
    plot_est_collected[[5]] + plot_re_collected[[5]] + plot_est_collected[[6]] + plot_re_collected[[6]] +
    plot_est_collected[[7]] + plot_re_collected[[7]] + guide_area() +
    plot_layout(nrow = 4, ncol = 4, byrow = FALSE, height = rep(c(2, 1), 7), guides = "collect")
dev.off()
