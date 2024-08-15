# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggpubr))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <-
  glue("{path_data}/../output/paper_output/metrics")

input_metrics_table <-
  fread(glue("{path_output}/fig6_fig7_accuracy_f1_score_bray_curtis_contamination_ratio_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

# 2. figure --------------------------------------------------------------------

data_boxplot_input <-
  input_metrics_table %>%
  mutate(contamination_90 = ifelse(Percent_Contaminants >= 90, "90 UP", "90 DOWN")) %>%
  mutate(contamination_90 = factor(contamination_90, levels = c("90 DOWN", "90 UP")))

value_comparisons <- list(c("SCRuB", "Green Cleaner"))

fig_f1_score <- 
  ggboxplot(data_boxplot_input %>%
              rename("Type" = "type") %>%
              filter(metrics == "f1_score"), x = "Type", y = "value",
            color = "Type", palette = "npg",
            add = "jitter",
            facet.by = "contamination_90", short.panel.labs = T) +
  stat_compare_means(comparisons = value_comparisons, label.y = 1, method = "wilcox", paired = T, label = "p.signif") +
  theme_bw() + 
  scale_color_manual(values = c("SCRuB" = "blue3", "Green Cleaner" = "green3")) +
  labs(x = "Type", y = "F1-score")

fig_diss <- 
  ggboxplot(data_boxplot_input %>%
              rename("Type" = "type") %>%
              filter(metrics == "bray_curtis"), x = "Type", y = "value",
            color = "Type", palette = "npg",
            add = "jitter",
            facet.by = "contamination_90", short.panel.labs = T) +
  stat_compare_means(comparisons = value_comparisons, label.y = 1, method = "wilcox", paired = T, label = "p.signif") +
  theme_bw() + 
  scale_color_manual(values = c("SCRuB" = "blue3", "Green Cleaner" = "green3")) +
  labs(x = "Type", y = "Bray-Curtis")

fig_merged <- ggarrange(fig_f1_score, fig_diss,
                        nrow = 1, common.legend = T, legend = "right", labels = c("A", "B"))

# save -------------------------------------------------------------------------

png(glue("{path_output}/fig8_f1_score_bray_curtis_90_up_down_boxplot.png"), width = 1000, height = 400)

print(fig_merged)

dev.off()