# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

rm(list = ls())

options(scipen = 999)

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_dataset_decontamed <-
  glue("{path_data}/../output/dataset_output/asv_count")
path_output <-
  glue("{path_data}/../output/paper_output/removed_asv_ratio")

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_sample_asv_count <-
  fread(glue("{path_data}/input_raw_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_scrub_asv_count <-
  fread(glue("{path_data}/input_scrub_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

# 2. pre processing ------------------------------------------------------------

data_selected_meta <- 
  input_meta %>%
  filter(dataset != "-") %>%
  select(sample_id, dilution_series, dataset) %>%
  mutate(removed_scrub = 0) %>%
  mutate(removed_gc = 0)

for (value_num in seq(nrow(data_selected_meta))) {
  
  value_tmp_dataset <- 
    data_selected_meta$dataset[value_num]
  value_tmp_dilution <- 
    data_selected_meta$dilution_series[value_num]
  value_tmp_sample <- 
    data_selected_meta$sample_id[value_num]
  
  data_tmp_d0_info <-
    data_selected_meta %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series == "d0")
  
  if ((value_tmp_dilution != "d0") & (value_tmp_dilution != "-")) {
    
    value_tmp_scrub_total_asv_count <- 
      input_scrub_asv_count %>%
      select(all_of(value_tmp_sample)) %>%
      colSums() %>%
      as.numeric()
    
    value_tmp_gc_total_asv_count <- 
      fread(glue("{path_dataset_decontamed}/{value_tmp_dataset}_merged_asv_count.txt"),
            sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F) %>%
      select(all_of(value_tmp_sample)) %>%
      colSums() %>%
      as.numeric()
    
    value_tmp_raw_total_asv_count <- 
      input_sample_asv_count %>%
      select(all_of(value_tmp_sample)) %>%
      colSums() %>%
      as.numeric()
    
    data_selected_meta[data_selected_meta$sample_id == value_tmp_sample, ]$removed_scrub <- 
      round((value_tmp_raw_total_asv_count - value_tmp_scrub_total_asv_count) / value_tmp_raw_total_asv_count, 4) * 100
    data_selected_meta[data_selected_meta$sample_id == value_tmp_sample, ]$removed_gc <- 
      round((value_tmp_raw_total_asv_count - value_tmp_gc_total_asv_count) / value_tmp_raw_total_asv_count, 4) * 100
    
  }
  
}

data_removed_asv_ratio <- 
  data_selected_meta %>%
  filter(dilution_series != "d0") %>%
  filter(dilution_series != "-")

write.table(data_removed_asv_ratio, glue("{path_output}/removed_asv_ratio_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)