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

input_all_metrics <-
  fread(glue("{path_output}/metrics_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_cotamination_ratio <-
  fread(glue("{path_output}/../contamination_ratio/fig2_contamination_ratio_table.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

# 2. pre-processing ----------------------------------------------------

data_specificity_plot_input <- 
  input_all_metrics %>%
  filter(type == "specificity") %>%
  select(-type) %>%
  pivot_longer(cols = c("SCRuB", "Green Cleaner"), names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("SCRuB", "Green Cleaner"))) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}")))

data_ppv_plot_input <- 
  input_all_metrics %>%
  filter(type == "PPV") %>%
  select(-type) %>%
  pivot_longer(cols = c("SCRuB", "Green Cleaner"), names_to = "type", values_to = "value") %>%
  mutate(type = factor(type, levels = c("SCRuB", "Green Cleaner"))) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}")))

data_tmp_merged_table <- 
  rbind(data_specificity_plot_input %>%
          mutate(metrics = "specificity"),
        data_ppv_plot_input %>%
          mutate(metrics = "PPV")) %>%
  mutate(key_value = glue('{gsub(pattern = "Dataset ", replacement = "dataset_", Dataset)}_D{d_stage}'))

data_merged_table <-
  merge(data_tmp_merged_table, 
        input_cotamination_ratio %>%
          mutate(key_value = glue("{dataset}_{Dilution_Series}")) %>%
          select(key_value, Percent_Contaminants),
        by = "key_value") %>%
  select(-key_value)

# 3. figure --------------------------------------------------------------------

fig_line_specificity <- 
  ggplot(data_merged_table %>%
           filter(metrics == "specificity"), 
         aes(x = Percent_Contaminants, y = value, color = type)) +
  geom_point() + 
  geom_smooth(method = "loess", se = F) + 
  labs(x = "Proportion of total contaminants", y = "Specificity") +
  theme_bw() +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_line_ppv <- 
  ggplot(data_merged_table %>%
           filter(metrics == "PPV"), 
         aes(x = Percent_Contaminants, y = value, color = type)) +
  geom_point() + 
  geom_smooth(method = "loess", se = F) + 
  labs(x = "Proportion of total contaminants", y = "PPV") +
  theme_bw() +
  scale_color_manual("Type", values = c("SCRuB" = "blue3", "Green Cleaner" = "green3"))

fig_merged <- 
  ggarrange(fig_line_specificity, fig_line_ppv,
            nrow = 1, common.legend = T, legend = "right", labels = c("A", "B", "C"))

# save -------------------------------------------------------------------------

png(glue("{path_output}/additional_file4_specificity_ppv_line_figure.png"), width = 900, height = 400)

print(fig_merged)

dev.off()

write.table(data_merged_table, glue("{path_output}/additional_file4_specificity_ppv_line_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)