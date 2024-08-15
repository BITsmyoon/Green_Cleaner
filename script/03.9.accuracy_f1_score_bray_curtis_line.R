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
path_taxa_relab <-
  glue("{path_data}/../output/dataset_output/taxa_relab")

input_all_metrics <-
  fread(glue("{path_output}/metrics_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_cotamination_ratio <-
  fread(glue("{path_output}/../contamination_ratio/fig2_contamination_ratio_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

func.dissim.cal <- 
  function(value_tmp_dataset) {
    
    # value_tmp_dataset <- "dataset_1"
    
    data_all_relab <- 
      fread(glue("{path_taxa_relab}/{value_tmp_dataset}_merged_taxa_relab.txt"), 
            sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
    
    data_sim_input <- 
      data_all_relab %>%
      column_to_rownames("taxa") %>%
      t() %>%
      as.data.frame()
    
    value_bc_dissimilarity <- 
      vegdist(data_sim_input, method = "bray") %>%
      as.matrix() %>%
      as.data.frame() %>%
      select(D0) %>%
      .[-1, ]
    
    data_bc_dissimilarity <-
      data.frame(Dataset = gsub(pattern = "dataset_", replacement = "Dataset ", x = value_tmp_dataset),
                 d_stage = rep(x = seq(7), 2),
                 type = c(rep("SCRuB", 7), rep("Green Cleaner", 7)),
                 value = value_bc_dissimilarity)
    
    return(data_bc_dissimilarity)
    
  }

# 2. pre-processing ----------------------------------------------------

data_bc_batch_1 <- 
  func.dissim.cal("dataset_1")
data_bc_batch_2 <- 
  func.dissim.cal("dataset_2")
data_bc_batch_3 <- 
  func.dissim.cal("dataset_3")
data_bc_batch_4 <- 
  func.dissim.cal("dataset_4")
data_bc_batch_5 <- 
  func.dissim.cal("dataset_5")
data_bc_batch_6 <- 
  func.dissim.cal("dataset_6")
data_bc_batch_7 <- 
  func.dissim.cal("dataset_7")
data_bc_batch_8 <- 
  func.dissim.cal("dataset_8")
data_bc_batch_9 <- 
  func.dissim.cal("dataset_9")
data_bc_batch_10 <- 
  func.dissim.cal("dataset_10")

data_accuracy_plot_input <- 
  input_all_metrics %>%
  filter(type == "Accuracy") %>%
  select(-type) %>%
  pivot_longer(cols = c("SCRuB", "Green Cleaner"), names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("SCRuB", "Green Cleaner"))) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}")))

data_f1_score_plot_input <- 
  input_all_metrics %>%
  filter(type == "F1-score") %>%
  select(-type) %>%
  pivot_longer(cols = c("SCRuB", "Green Cleaner"), names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("SCRuB", "Green Cleaner"))) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}")))

data_diss_plot_input <- 
  rbind(data_bc_batch_1, data_bc_batch_2) %>%
  rbind(., data_bc_batch_3) %>%
  rbind(., data_bc_batch_4) %>%
  rbind(., data_bc_batch_5) %>%
  rbind(., data_bc_batch_6) %>%
  rbind(., data_bc_batch_7) %>%
  rbind(., data_bc_batch_8) %>%
  rbind(., data_bc_batch_9) %>%
  rbind(., data_bc_batch_10) %>%
  mutate(type = factor(type, levels = c("SCRuB", "Green Cleaner"))) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}")))

data_tmp_merged_table <- 
  rbind(data_accuracy_plot_input %>%
          mutate(metrics = "accuracy"),
        data_f1_score_plot_input %>%
          mutate(metrics = "f1_score")) %>%
  rbind(., data_diss_plot_input %>%
          mutate(metrics = "bray_curtis")) %>%
  mutate(key_value = glue('{gsub(pattern = "Dataset ", replacement = "dataset_", Dataset)}_D{d_stage}'))

data_merged_table <-
  merge(data_tmp_merged_table, 
        input_cotamination_ratio %>%
          mutate(key_value = glue("{dataset}_{Dilution_Series}")) %>%
          select(key_value, Percent_Contaminants),
        by = "key_value") %>%
  select(-key_value)

# 3. figure --------------------------------------------------------------------

fig_accuracy <- 
  ggplot(data_accuracy_plot_input, aes(x = d_stage, y = value, colour = type, group = type)) +
  geom_line() +  
  geom_point() +  
  labs(x = "Dilution series", y = "Accuracy") +
  theme_bw() +
  facet_wrap(~ Dataset, nrow = 5) +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_f1_score <- 
  ggplot(data_f1_score_plot_input, aes(x = d_stage, y = value, colour = type, group = type)) +
  geom_line() +  
  geom_point() +  
  labs(x = "Dilution series", y = "F1-score") +
  theme_bw() +
  facet_wrap(~ Dataset, nrow = 5) +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_diss <- 
  ggplot(data_diss_plot_input, aes(x = d_stage, y = value, colour = type, group = type)) +
  geom_line() +  
  geom_point() +  
  labs(x = "Dilution series", y = "Bray-Curtis") +
  theme_bw() +
  facet_wrap(~ Dataset, nrow = 5) +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_6_merged <- 
  ggarrange(fig_accuracy, fig_f1_score, fig_diss,
            nrow = 1, common.legend = T, legend = "right", labels = c("A", "B", "C"))

fig_line_accuracy <- 
  ggplot(data_merged_table %>%
           filter(metrics == "accuracy"), 
         aes(x = Percent_Contaminants, y = value, color = type)) +
  geom_point() + 
  geom_smooth(method = "loess", se = F) + 
  labs(x = "Proportion of total contaminants", y = "Accuracy") +
  theme_bw() +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_line_f1_score <- 
  ggplot(data_merged_table %>%
           filter(metrics == "f1_score"), 
         aes(x = Percent_Contaminants, y = value, color = type)) +
  geom_point() + 
  geom_smooth(method = "loess", se = F) + 
  labs(x = "Proportion of total contaminants", y = "F1-score") +
  theme_bw() +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_line_diss <- 
  ggplot(data_merged_table %>%
           filter(metrics == "bray_curtis"), 
         aes(x = Percent_Contaminants, y = value, color = type)) +
  geom_point() + 
  geom_smooth(method = "loess", se = F) + 
  labs(x = "Proportion of total contaminants", y = "Bray-Curtis") +
  theme_bw() +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_7_merged <- 
  ggarrange(fig_line_accuracy, fig_line_f1_score, fig_line_diss,
            nrow = 1, common.legend = T, legend = "right", labels = c("A", "B", "C"))

# save -------------------------------------------------------------------------

png(glue("{path_output}/fig6_accuracy_f1_score_bray_curtis_dot_figure.png"), width = 1200, height = 800)

print(fig_6_merged)

dev.off()

png(glue("{path_output}/fig7_accuracy_f1_score_bray_curtis_line_figure.png"), width = 1000, height = 400)

print(fig_7_merged)

dev.off()

write.table(data_merged_table, glue("{path_output}/fig6_fig7_accuracy_f1_score_bray_curtis_contamination_ratio_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)