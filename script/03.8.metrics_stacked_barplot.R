# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

# plot
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

rm(list = ls())

options(scipen = 999)

# 1. Variables and Function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <-
  glue("{path_data}/../output/paper_output/metrics")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

input_metrics <-
  fread(glue("{path_output}/metrics_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

# 2. pre processing ----------------------------------------------------

data_metrics <- 
  input_metrics %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}"))) %>%
  mutate(d_stage = glue("D{d_stage}")) %>%
  mutate(d_stage = factor(d_stage, levels = glue("D{seq(7)}"))) %>%
  filter(type %in% c("TP", "TN", "FP", "FN")) %>%
  mutate(type = factor(type, levels = c("TP", "TN", "FP", "FN"))) %>%
  mutate(scrub_relab = 0) %>%
  mutate(gc_relab = 0)

value_dataset_list <- glue("Dataset {seq(10)}")
value_d_stage_list <- glue("D{seq(7)}")

for (value_dataset in value_dataset_list) {
  
  # value_dataset <- value_dataset_list[1]
  
  data_tmp_dataset_info <- 
    data_metrics %>%
    filter(Dataset == value_dataset)
  
  for (value_d_stage in value_d_stage_list) {
    
    # value_d_stage <- value_d_stage_list[1]
    
    data_tmp_d_stage_info <-
      data_tmp_dataset_info %>%
      filter(d_stage == value_d_stage)
    
    value_sum_scrub <- 
      sum(data_tmp_d_stage_info$SCRuB)
    
    value_sum_gc <- 
      sum(data_tmp_d_stage_info$`Green Cleaner`)
    
    value_relab_scrub <- 
      data_tmp_d_stage_info$SCRuB / value_sum_scrub
    
    value_relab_gc <- 
      data_tmp_d_stage_info$`Green Cleaner` / value_sum_gc
    
    logic_row <-
      (data_metrics$Dataset == value_dataset) & (data_metrics$d_stage == value_d_stage)
    
    data_metrics$scrub_relab[logic_row] <- value_relab_scrub
    data_metrics$gc_relab[logic_row] <- value_relab_gc
    
  }
  
}

# 3. plot ----------------------------------------------------

data_scrub_plot <- 
  data_metrics %>%
  select(Dataset, d_stage, type, scrub_relab) %>%
  rename("Proportion" = "scrub_relab")

data_gc_plot <- 
  data_metrics %>%
  select(Dataset, d_stage, type, gc_relab) %>%
  rename("Proportion" = "gc_relab")
  
fig_scrub <- 
  ggplot(data_scrub_plot, aes(x = d_stage, y = Proportion, fill = type)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_wrap(~ Dataset, nrow = 2) +
  scale_fill_manual("Class", values = c("#F8766D", "#619CFF", "#BDBDBD", "#6C6C6C")) +
  labs(x = "Dilution series")

fig_gc <-
  ggplot(data_gc_plot, aes(x = d_stage, y = Proportion, fill = type)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  facet_wrap(~ Dataset, nrow = 2) +
  scale_fill_manual("Class", values = c("#F8766D", "#619CFF", "#BDBDBD", "#6C6C6C")) +
  labs(x = "Dilution series")

fig_merged <- ggarrange(fig_scrub, fig_gc, 
                        nrow = 2, common.legend = T, legend = "right", labels = c("A", "B"))

png(glue("{path_output}/fig5_metrics_stacked_barplot_figure.png"), width = 600, height = 600, pointsize = 12)

print(fig_merged)

dev.off()

write.table(data_metrics, glue("{path_output}/fig5_metrics_stacked_barplot_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)