# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

input_meta_data <- 
  fread(glue("{path_data}/input_sample_meta_data.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_asv_count <-
  fread(glue("{path_data}/input_raw_asv_count.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

input_dataset_list <- 
  glue("dataset_{seq(10)}")

# 2. processing ----------------------------------------------------------------

for (value_tmp_dataset in input_dataset_list) {
  
  # value_tmp_dataset <- input_dataset_list[1]
  
  path_output <- 
    glue("{path_data}/../output/dataset_output/asv_count")
  
  if (!dir.exists(path_output)) {
    
    dir.create(path = path_output, recursive = T)
    dir.create(path = glue("{path_output}/../meta_data_for_qiime"), recursive = T)
    
  }
  
  data_ntc_meta <-
    input_meta_data %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series == "-")
  
  data_d0_meta <-
    input_meta_data %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series == "d0")
  
  data_sample_meta <-
    input_meta_data %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series != "d0") %>%
    filter(dilution_series != "-")

  value_batch_name <-
    unique(data_sample_meta$decontamination_batch)
  
  data_ntc_asv_count <- 
    input_asv_count %>%
    select(asv_id, data_ntc_meta$sample_id)
  
  data_d0_asv_count <- 
    input_asv_count %>%
    select(asv_id, data_d0_meta$sample_id)
  
  data_decontamed_asv_count <- 
    fread(glue("{path_data}/../output/batch_decontaminated/{value_batch_name}_Green_Cleaner_decontaminated_asv_count.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
  
  data_sample_asv_count <- 
    data_decontamed_asv_count %>%
    select(asv_id, data_sample_meta$sample_id)
  
  data_inter_merged_asv_count <- 
    full_join(data_d0_asv_count, data_sample_asv_count, by = "asv_id") %>%
    full_join(., data_ntc_asv_count, by = "asv_id")
  
  data_inter_merged_asv_count[is.na(data_inter_merged_asv_count)] <- 0
  
  logic_rowsum_not_zero <- 
    data_inter_merged_asv_count %>%
    column_to_rownames("asv_id") %>%
    rowSums() != 0
  
  data_final_asv_count <- 
    data_inter_merged_asv_count[logic_rowsum_not_zero, ]
  
  write.table(data_final_asv_count, glue("{path_output}/{value_tmp_dataset}_merged_asv_count.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = T)
  
  value_sample_name <- 
    data_final_asv_count %>%
    select(-asv_id) %>%
    colnames()
  
  data_meta <- 
    data.frame(sample = value_sample_name,
             dataset = value_tmp_dataset) %>%
    rename("#SampleID" = "sample")
  
  write.table(data_meta, glue("{path_output}/../meta_data_for_qiime/{value_tmp_dataset}_metadata.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = T)
  
}
