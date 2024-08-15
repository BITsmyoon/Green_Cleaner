# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))

rm(list = ls())

options(scipen = 999)

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <-
  glue("{path_data}/../output/paper_output/metrics")
path_taxa_relab_output <-
  glue("{path_data}/../output/dataset_output/taxa_relab")
path_dataset_asv_count <-
  glue("{path_data}/../output/dataset_output/asv_count")
path_qiime <-
  glue("{path_data}/../output/qiime")

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_sample_asv_count <-
  fread(glue("{path_data}/input_raw_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_scrub_asv_count <-
  fread(glue("{path_data}/input_scrub_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

input_raw_species_count <-
  read.csv(glue("{path_data}/../output/qiime/qiime_barplot/level-7.csv"), 
           sep = ",", header = T, check.names = F, stringsAsFactors = F)
input_scrub_species_count <-
  read.csv(glue("{path_qiime}/scrub_qiime_barplot/level-7.csv"), 
           sep = ",", header = T, check.names = F, stringsAsFactors = F)

input_dataset_list <- glue("dataset_{seq(10)}")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  dir.create(path = path_taxa_relab_output, recursive = T)
  
}

func.confusion.check <- function(value){
  
  value_tmp_true_signal <- value[2]
  value_tmp_contam_check <- value[3]
  
  if ((value_tmp_true_signal == T) & (value_tmp_contam_check == T)) {
    
    value_tmp_check <- "TP"
    
  } else if ((value_tmp_true_signal == F) & (value_tmp_contam_check == F)) {
    
    value_tmp_check <- "TN"
    
  } else if ((value_tmp_true_signal == F) & (value_tmp_contam_check == T)) {
    
    value_tmp_check <- "FP"
    
  } else if ((value_tmp_true_signal == T) & (value_tmp_contam_check == F)) {
    
    value_tmp_check <- "FN"
    
  }
  
  return(value_tmp_check)
  
}

# 2. metrics -------------------------------------------------------------------

for (value_tmp_batch in input_dataset_list) {
  
  # value_tmp_batch <- input_dataset_list[1]
  
  data_tmp_meta <-
    input_meta %>%
    filter(dataset == value_tmp_batch) %>%
    filter(ntc_check != "ntc")
  
  value_not_d0_sample_name <- 
    data_tmp_meta %>%
    filter(dilution_series != "d0") %>%
    pull("sample_id")
  value_d0_sample_name <- 
    data_tmp_meta %>%
    filter(dilution_series == "d0") %>%
    pull("sample_id")
  
  data_d0_asv_count <- 
    input_sample_asv_count %>%
    select(asv_id, all_of(value_d0_sample_name)) %>%
    rename("D0" = 2)
  data_scrub_asv_count <- 
    input_scrub_asv_count %>%
    select(asv_id, all_of(value_not_d0_sample_name))
  colnames(data_scrub_asv_count) <- 
    c("asv_id", glue("D{seq(7)}_SCRuB"))
  
  data_gc_asv_count <-
    fread(glue("{path_dataset_asv_count}/{value_tmp_batch}_merged_asv_count.txt"), 
          sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F) %>%
    select(asv_id, all_of(value_not_d0_sample_name))
  colnames(data_gc_asv_count) <- 
    c("asv_id", glue("D{seq(7)}_Green_Cleaner"))
  
  data_tmp_info <- 
    data.frame(sample_name = c("D0", glue("D{seq(7)}_SCRuB"), glue("D{seq(7)}_Green_Cleaner")),
               type = c("D0", rep("SCRuB", 7), rep("Green_Cleaner", 7)),
               d_stage = c("D0", rep(seq(7), 2)))
  
  data_tmp_simulated_info <- 
    data.frame(sample_name = glue("D{seq(7)}_contamed"),
               type = rep("Contamed", 7),
               d_stage = seq(7))
  
  data_selected_tmp_simulated_asv_count <- 
    input_sample_asv_count %>%
    select(c("asv_id", all_of(value_not_d0_sample_name))) %>%
    column_to_rownames("asv_id")
  colnames(data_selected_tmp_simulated_asv_count) <- 
    glue("D{seq(7)}_contamed")
  
  data_selected_tmp_asv_count <- 
    full_join(data_d0_asv_count, data_scrub_asv_count, by = "asv_id") %>%
    full_join(., data_gc_asv_count, by = "asv_id")
  
  data_selected_tmp_asv_count[is.na(data_selected_tmp_asv_count)] <- 0
  
  data_tmp_true_signal <- 
    data_selected_tmp_asv_count %>%
    select(asv_id, 2) %>%
    rename("D0" = 2) %>%
    mutate(true_signal = ifelse(D0 != 0, T, F)) %>%
    select(asv_id, true_signal)
  
  data_all_confusion_scrub <- data.frame()
  data_all_confusion_gc <- data.frame()
  
  for (value_tmp_d_stage in seq(7)) {
    
    # value_tmp_d_stage <- 1
    
    data_terget_simulated_count <- 
      data_selected_tmp_simulated_asv_count %>%
      select(all_of(value_tmp_d_stage))
    
    value_tmp_scrub_sample_name <- 
      data_tmp_info %>%
      filter(d_stage == value_tmp_d_stage) %>%
      filter(type == "SCRuB") %>%
      pull("sample_name")
    value_tmp_gc_sample_name <- 
      data_tmp_info %>%
      filter(d_stage == value_tmp_d_stage) %>%
      filter(type == "Green_Cleaner") %>%
      pull("sample_name")
    
    data_tmp_scrub_count <- 
      data_selected_tmp_asv_count[, c("asv_id", value_tmp_scrub_sample_name)] %>%
      rename("scrub" = 2)
    
    data_tmp_scrub_add_contam_check <- 
      data_tmp_scrub_count %>%
      mutate(contam_check = ifelse(scrub != 0, T, F)) %>%
      select(asv_id, contam_check)
    
    data_tmp_merged_scrub <- 
      cbind(data_tmp_true_signal, data_tmp_scrub_add_contam_check$contam_check)
    
    data_tmp_merged_scrub$confusion_check <-
      apply(data_tmp_merged_scrub, 1, func.confusion.check)
    
    data_tmp_gc_count <- 
      data_selected_tmp_asv_count[, c("asv_id", value_tmp_gc_sample_name)] %>%
      rename("gc" = 2)
    
    data_tmp_gc_add_contam_check <- 
      data_tmp_gc_count %>%
      mutate(contam_check = ifelse(gc != 0, T, F)) %>%
      select(asv_id, contam_check)
    
    data_tmp_merged_gc <- 
      cbind(data_tmp_true_signal, data_tmp_gc_add_contam_check$contam_check)
    
    data_tmp_merged_gc$confusion_check <-
      apply(data_tmp_merged_gc, 1, func.confusion.check)
    
    for (value_tmp_target in c("TP", "TN", "FP", "FN", "accuracy", "sensitivity", "specificity", "ppv", "f1_score")) {

      # value_tmp_target <- "FP"

      if (value_tmp_target %in% c("TP", "TN", "FP", "FN")) {

        value_tmp_target_asv <- 
          data_tmp_merged_scrub %>%
          filter(confusion_check == value_tmp_target) %>%
          pull("asv_id")
        
        # length(data_terget_simulated_count[value_tmp_target_asv, ])
        
        value_target_sum <-
          sum(data_terget_simulated_count[value_tmp_target_asv, ])

        data_tmp_confusion_scrub <- data.frame(batch = value_tmp_batch,
                                               d_stage = value_tmp_d_stage,
                                               type = value_tmp_target,
                                               value = value_target_sum)

        value_tmp_target_asv <- 
          data_tmp_merged_gc %>%
          filter(confusion_check == value_tmp_target) %>%
          pull("asv_id")
        
        value_target_sum <-
          sum(data_terget_simulated_count[value_tmp_target_asv, ])

        data_tmp_confusion_gc <- data.frame(batch = value_tmp_batch,
                                            d_stage = value_tmp_d_stage,
                                            type = value_tmp_target,
                                            value = value_target_sum)

        if ((value_tmp_d_stage == 1) & (value_tmp_target == "TP")) {

          data_all_confusion_scrub <- data_tmp_confusion_scrub
          data_all_confusion_gc <- data_tmp_confusion_gc

        } else {

          data_all_confusion_scrub <- rbind(data_all_confusion_scrub, data_tmp_confusion_scrub)
          data_all_confusion_gc <- rbind(data_all_confusion_gc, data_tmp_confusion_gc)

        }

      } else if (value_tmp_target == "accuracy") {

        value_tmp_scrub_under <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          pull("value") %>%
          sum()

        value_tmp_scrub_upper <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "TN")) %>%
          pull("value") %>%
          sum()

        data_all_confusion_scrub <-
          rbind(data_all_confusion_scrub, data.frame(batch = value_tmp_batch,
                                                     d_stage = value_tmp_d_stage,
                                                     type = "accuracy",
                                                     value = as.numeric(value_tmp_scrub_upper / value_tmp_scrub_under)))

        value_tmp_gc_under <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          pull("value") %>%
          sum()

        value_tmp_gc_upper <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "TN")) %>%
          pull("value") %>%
          sum()

        data_all_confusion_gc <-
          rbind(data_all_confusion_gc, data.frame(batch = value_tmp_batch,
                                                  d_stage = value_tmp_d_stage,
                                                  type = "accuracy",
                                                  value = as.numeric(value_tmp_gc_upper / value_tmp_gc_under)))

      } else if (value_tmp_target == "sensitivity") {

        value_tmp_scrub_under <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "FN")) %>%
          pull("value") %>%
          sum()

        value_tmp_scrub_upper <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP")) %>%
          pull("value") %>%
          sum()

        data_all_confusion_scrub <-
          rbind(data_all_confusion_scrub, data.frame(batch = value_tmp_batch,
                                                     d_stage = value_tmp_d_stage,
                                                     type = "sensitivity",
                                                     value = as.numeric(value_tmp_scrub_upper / value_tmp_scrub_under)))

        value_tmp_gc_under <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "FN")) %>%
          pull("value") %>%
          sum()

        value_tmp_gc_upper <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP")) %>%
          pull("value") %>%
          sum()

        data_all_confusion_gc <-
          rbind(data_all_confusion_gc, data.frame(batch = value_tmp_batch,
                                                  d_stage = value_tmp_d_stage,
                                                  type = "sensitivity",
                                                  value = as.numeric(value_tmp_gc_upper / value_tmp_gc_under)))

      } else if (value_tmp_target == "specificity") {
        
        value_tmp_scrub_under <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("FP", "TN")) %>%
          pull("value") %>%
          sum()
        
        value_tmp_scrub_upper <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TN")) %>%
          pull("value") %>%
          sum()
        
        data_all_confusion_scrub <-
          rbind(data_all_confusion_scrub, data.frame(batch = value_tmp_batch,
                                                     d_stage = value_tmp_d_stage,
                                                     type = "specificity",
                                                     value = as.numeric(value_tmp_scrub_upper / value_tmp_scrub_under)))
        
        value_tmp_gc_under <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("FP", "TN")) %>%
          pull("value") %>%
          sum()
        
        value_tmp_gc_upper <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TN")) %>%
          pull("value") %>%
          sum()
        
        data_all_confusion_gc <-
          rbind(data_all_confusion_gc, data.frame(batch = value_tmp_batch,
                                                  d_stage = value_tmp_d_stage,
                                                  type = "specificity",
                                                  value = as.numeric(value_tmp_gc_upper / value_tmp_gc_under)))
        
      } else if (value_tmp_target == "ppv") {
        
        value_tmp_scrub_under <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "FP")) %>%
          pull("value") %>%
          sum()
        
        value_tmp_scrub_upper <-
          data_all_confusion_scrub %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP")) %>%
          pull("value") %>%
          sum()
        
        data_all_confusion_scrub <-
          rbind(data_all_confusion_scrub, data.frame(batch = value_tmp_batch,
                                                     d_stage = value_tmp_d_stage,
                                                     type = "ppv",
                                                     value = as.numeric(value_tmp_scrub_upper / value_tmp_scrub_under)))
        
        value_tmp_gc_under <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP", "FP")) %>%
          pull("value") %>%
          sum()
        
        value_tmp_gc_upper <-
          data_all_confusion_gc %>%
          filter(d_stage == value_tmp_d_stage) %>%
          filter(type %in% c("TP")) %>%
          pull("value") %>%
          sum()
        
        data_all_confusion_gc <-
          rbind(data_all_confusion_gc, data.frame(batch = value_tmp_batch,
                                                  d_stage = value_tmp_d_stage,
                                                  type = "ppv",
                                                  value = as.numeric(value_tmp_gc_upper / value_tmp_gc_under)))
        
      } else if (value_tmp_target == "f1_score") {

        value_tmp_scrub_precision <-
          (data_all_confusion_scrub %>%
             filter(d_stage == value_tmp_d_stage) %>%
             filter(type %in% c("TP")) %>%
             pull("value") %>%
             sum()) / (data_all_confusion_scrub %>%
                         filter(d_stage == value_tmp_d_stage) %>%
                         filter(type %in% c("TP", "FP")) %>%
                         pull("value") %>%
                         sum())

        value_tmp_scrub_recall <-
          (data_all_confusion_scrub %>%
             filter(d_stage == value_tmp_d_stage) %>%
             filter(type %in% c("TP")) %>%
             pull("value") %>%
             sum()) / (data_all_confusion_scrub %>%
                         filter(d_stage == value_tmp_d_stage) %>%
                         filter(type %in% c("TP", "FN")) %>%
                         pull("value") %>%
                         sum())

        data_all_confusion_scrub <-
          rbind(data_all_confusion_scrub, data.frame(batch = value_tmp_batch,
                                                     d_stage = value_tmp_d_stage,
                                                     type = "f1_score",
                                                     value = (2 * value_tmp_scrub_precision * value_tmp_scrub_recall) / (value_tmp_scrub_precision + value_tmp_scrub_recall)))

        value_tmp_gc_precision <-
          (data_all_confusion_gc %>%
             filter(d_stage == value_tmp_d_stage) %>%
             filter(type %in% c("TP")) %>%
             pull("value") %>%
             sum()) / (data_all_confusion_gc %>%
                         filter(d_stage == value_tmp_d_stage) %>%
                         filter(type %in% c("TP", "FP")) %>%
                         pull("value") %>%
                         sum())

        value_tmp_gc_recall <-
          (data_all_confusion_gc %>%
             filter(d_stage == value_tmp_d_stage) %>%
             filter(type %in% c("TP")) %>%
             pull("value") %>%
             sum()) / (data_all_confusion_gc %>%
                         filter(d_stage == value_tmp_d_stage) %>%
                         filter(type %in% c("TP", "FN")) %>%
                         pull("value") %>%
                         sum())

        data_all_confusion_gc <-
          rbind(data_all_confusion_gc, data.frame(batch = value_tmp_batch,
                                                  d_stage = value_tmp_d_stage,
                                                  type = "f1_score",
                                                  value = (2 * value_tmp_gc_precision * value_tmp_gc_recall) / (value_tmp_gc_precision + value_tmp_gc_recall)))

      }

    }
    
  }
  
  if (value_tmp_batch == input_dataset_list[1]) {
    
    data_all_confusion_scrub_all_batch <- data_all_confusion_scrub
    data_all_confusion_gc_all_batch <- data_all_confusion_gc
    
  } else {
    
    data_all_confusion_scrub_all_batch <- rbind(data_all_confusion_scrub_all_batch, data_all_confusion_scrub)
    data_all_confusion_gc_all_batch <- rbind(data_all_confusion_gc_all_batch, data_all_confusion_gc)
    
  }
  
}

data_all_confusion_scrub_all_batch$type <-
  gsub(pattern = "accuracy", replacement = "Accuracy", x = data_all_confusion_scrub_all_batch$type)
data_all_confusion_scrub_all_batch$type <-
  gsub(pattern = "f1_score", replacement = "F1-score", x = data_all_confusion_scrub_all_batch$type)
data_all_confusion_scrub_all_batch$type <-
  gsub(pattern = "ppv", replacement = "PPV", x = data_all_confusion_scrub_all_batch$type)

data_all_confusion_gc_all_batch$type <-
  gsub(pattern = "accuracy", replacement = "Accuracy", x = data_all_confusion_gc_all_batch$type)
data_all_confusion_gc_all_batch$type <-
  gsub(pattern = "f1_score", replacement = "F1-score", x = data_all_confusion_gc_all_batch$type)
data_all_confusion_gc_all_batch$type <-
  gsub(pattern = "ppv", replacement = "PPV", x = data_all_confusion_gc_all_batch$type)

data_all_merge <- 
  cbind(data_all_confusion_scrub_all_batch %>%
          rename("SCRuB" = "value"), 
        data_all_confusion_gc_all_batch %>%
          select(value) %>%
          rename("Green Cleaner" = "value")) %>%
  rename("Dataset" = "batch") %>%
  mutate(Dataset = gsub(pattern = "dataset_", replacement = "Dataset ", Dataset)) %>%
  mutate(Dataset = factor(Dataset, levels = glue("Dataset {seq(10)}"))) %>%
  mutate(d_stage = as.factor(d_stage))

write.table(data_all_merge, glue("{path_output}/metrics_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)

# 3. dissimilarity -------------------------------------------------------------

data_raw_species_relab <- 
  input_raw_species_count %>%
  select(-dataset) %>%
  column_to_rownames("index") %>%
  apply(., 1, function(value){
    (value / sum(value)) * 100
  }) %>%
  as.data.frame() %>%
  rownames_to_column("taxa")

data_scrub_species_relab <- 
  input_scrub_species_count %>%
  select(-dataset, -dilution_series) %>%
  column_to_rownames("index") %>%
  apply(., 1, function(value){
    (value / sum(value)) * 100
  }) %>%
  as.data.frame() %>%
  rownames_to_column("taxa")

for (value_tmp_dataset in input_dataset_list) {
  
  # value_tmp_dataset <- input_dataset_list[1]
  
  value_tmp_d0_name <- 
    input_meta %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series == "d0") %>%
    pull("sample_id")
  
  value_tmp_sample_meta <-
    input_meta %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series != "d0") %>%
    filter(dilution_series != "-")
  
  data_tmp_species_relab <- 
    read.csv(glue("{path_qiime}/dataset/{value_tmp_dataset}_qiime_barplot/level-7.csv"), 
             header = T, sep = ",", check.names = F, stringsAsFactors = F) %>%
    select(-dataset) %>%
    column_to_rownames("index") %>%
    apply(., 1, function(value){
      (value / sum(value)) * 100
    }) %>%
    as.data.frame() %>%
    rownames_to_column("taxa")
  
  data_d0_relab <- 
    data_raw_species_relab %>%
    select(taxa, all_of(value_tmp_d0_name)) %>%
    rename("D0" = 2)
  
  data_scrub_relab <- 
    data_scrub_species_relab %>%
    select(taxa, value_tmp_sample_meta$sample_id) %>%
    rename("SCRuB_D1" = 2) %>%
    rename("SCRuB_D2" = 3) %>%
    rename("SCRuB_D3" = 4) %>%
    rename("SCRuB_D4" = 5) %>%
    rename("SCRuB_D5" = 6) %>%
    rename("SCRuB_D6" = 7) %>%
    rename("SCRuB_D7" = 8)
  
  data_gc_relab <- 
    data_tmp_species_relab %>%
    select(taxa, value_tmp_sample_meta$sample_id) %>%
    rename("Green_Cleaner_D1" = 2) %>%
    rename("Green_Cleaner_D2" = 3) %>%
    rename("Green_Cleaner_D3" = 4) %>%
    rename("Green_Cleaner_D4" = 5) %>%
    rename("Green_Cleaner_D5" = 6) %>%
    rename("Green_Cleaner_D6" = 7) %>%
    rename("Green_Cleaner_D7" = 8)
  
  data_merged_species_relab <- 
    full_join(data_d0_relab,
              data_scrub_relab,
            by = "taxa") %>%
    full_join(., data_gc_relab,
              by = "taxa")
  
  data_merged_species_relab[is.na(data_merged_species_relab)] <- 0
  
  write.table(data_merged_species_relab, glue("{path_taxa_relab_output}/{value_tmp_dataset}_merged_taxa_relab.txt"), 
              quote = F, sep = "\t", row.names = F, col.names = T)
  
}