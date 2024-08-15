# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))

# Plotting
suppressPackageStartupMessages(library(ggplot2))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_raw_asv_count <-
  fread(glue("{path_data}/input_raw_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

path_output <- 
  glue("{path_data}/../output/paper_output/contamination_ratio")
path_asv_count <-
  glue("{path_data}/../output/dataset_output")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

input_dataset_list <- glue("dataset_{seq(10)}")

# 2. pre processing ------------------------------------------------------------

for (value_tmp_dataset in input_dataset_list) {
  
  # value_tmp_dataset <- input_dataset_list[1]
  
  data_selected_dataset_meta <-
    input_meta %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series != "-")
  
  value_d0_sample_name <- 
    data_selected_dataset_meta %>%
    filter(dilution_series == "d0") %>%
    pull("sample_id")
  
  data_selected_asv_count <- 
    input_raw_asv_count %>%
    select("asv_id", data_selected_dataset_meta$sample_id) %>%
    column_to_rownames("asv_id")
  
  logic_contam <- 
    data_selected_asv_count[, value_d0_sample_name] == 0
  
  value_tmp_total_sum <- 
    colSums(data_selected_asv_count) %>%
    as.numeric()
  
  value_tmp_contam_sum <- 
    data_selected_asv_count[logic_contam, ] %>%
    colSums() %>%
    as.numeric()
  
  data_tmp_merged_info <- 
    data.frame(Dilution_Series = glue("D{0:7}"),
               Percent_Contaminants = as.numeric((value_tmp_contam_sum / value_tmp_total_sum) * 100),
               dataset = rep(value_tmp_dataset, 8))
  
  if (value_tmp_dataset == input_dataset_list[1]) {
    
    data_merged_info <- data_tmp_merged_info
    
  } else {
    
    data_merged_info <- rbind(data_merged_info, data_tmp_merged_info)
    
  }
  
}

data_modi_merged_info <-
  data_merged_info %>%
  mutate(Dilution_Series = factor(Dilution_Series)) %>%
  mutate(dataset = gsub("_", " ", dataset)) %>%
  mutate(dataset = gsub("data", "Data", dataset)) %>%
  mutate(dataset = factor(dataset, levels = glue("Dataset {seq(10)}")))

# 3. plot ------------------------------------------------------------

png(glue("{path_output}/fig2_contamination_ratio_figure.png"), width = 700, height = 400)

ggplot(data_modi_merged_info, aes(x = Dilution_Series, y = Percent_Contaminants)) + 
  geom_point() +
  facet_wrap(vars(dataset), ncol = 5) +
  theme_bw() + 
  ylab("Proportion of total contaminants") +
  xlab("Dilution series")

dev.off()

write.table(data_merged_info, glue("{path_output}/fig2_contamination_ratio_table.txt"), quote = F, sep = "\t", row.names = F, col.names = T)