# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

# Plotting
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <-
  glue("{path_data}/../output/paper_output/category")
path_gc_dataset_asv <-
  glue("{path_data}/../output/dataset_output/asv_count")

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_raw_asv_count <-
  fread(glue("{path_data}/input_raw_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_scrub_asv_count <-
  fread(glue("{path_data}/input_scrub_asv_count.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_dataset_list <- glue("dataset_{seq(10)}")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

func.category.check <- 
  function(value) {
    
    category_1_check <- value[3]
    category_2_check <- value[4]
    category_3_check <- value[5]
    
    if ((category_1_check == T) & (category_2_check == F) & (category_3_check == F)) {
      value_category <- "category_1"
    } else if ((category_1_check == F) & (category_2_check == T) & (category_3_check == F)) {
      value_category <- "category_2"
    } else if ((category_1_check == F) & (category_2_check == F) & (category_3_check == T)) {
      value_category <- "category_3"
    } else {
      value_category <- "not_match"
    }
    
    return(value_category)
    
  }

# 2. pre processing ----------------------------------------------------

for (value_tmp_dataset in input_dataset_list) {
  
  # value_tmp_dataset <- input_dataset_list[2]
  
  data_tmp_meta <-
    input_meta %>%
    filter(dataset == value_tmp_dataset)
  
  value_ntc_sample_name <-
    input_meta %>%
    filter(ntc_check == "ntc") %>%
    filter(grepl(pattern = value_tmp_dataset, dataset)) %>%
    pull("sample_id")
  value_d0_sample_name <-
    data_tmp_meta %>%
    filter(dilution_series == "d0") %>%
    pull("sample_id")
  value_target_sample_name <-
    data_tmp_meta %>%
    filter(dilution_series != "d0") %>%
    filter(ntc_check != "ntc") %>%
    pull("sample_id")
  
  data_ntc_asv_count <-
    input_raw_asv_count %>%
    select(asv_id, all_of(value_ntc_sample_name)) %>%
    rename("ntc" = 2)
  data_d0_asv_count <-
    input_raw_asv_count %>%
    select(asv_id, all_of(value_d0_sample_name)) %>%
    rename("D0" = 2)
  data_scrub_asv_count <- 
    input_scrub_asv_count %>%
    select(asv_id, all_of(value_target_sample_name))
  data_gc_asv_count <- 
    fread(glue("{path_gc_dataset_asv}/{value_tmp_dataset}_merged_asv_count.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F) %>%
    select(asv_id, all_of(value_target_sample_name))
  
  colnames(data_scrub_asv_count) <- 
    c("asv_id", glue("D{seq(7)}_SCRuB"))
  colnames(data_gc_asv_count) <- 
    c("asv_id", glue("D{seq(7)}_Green_Cleaner"))
  
  data_merged_asv_count <- 
    full_join(data_ntc_asv_count, data_d0_asv_count, by = "asv_id") %>%
    full_join(., data_scrub_asv_count, by = "asv_id") %>%
    full_join(., data_gc_asv_count, by = "asv_id")
  
  data_merged_asv_count[is.na(data_merged_asv_count)] <- 0
  
  data_ntc_asv_relab <-
    data_merged_asv_count %>%
    column_to_rownames("asv_id") %>%
    select(ntc) %>%
    apply(., 2, function(value){
      (value / sum(value)) * 100
    }) %>%
    as.data.frame()
  
  value_ntc_top5 <- 
    data_ntc_asv_relab %>%
    arrange(desc(ntc)) %>%
    head(5) %>%
    rownames()
  
  value_category_1_asv_name <-
    value_ntc_top5
  
  value_category_2_asv_name <-
    rownames(data_ntc_asv_relab)[(data_ntc_asv_relab$ntc != 0) & !(rownames(data_ntc_asv_relab) %in% value_ntc_top5)]
  
  value_merged_ntc_positive_asv_name <- 
    c(value_category_1_asv_name, value_category_2_asv_name)
  
  data_tmp_category_check <- 
    data_merged_asv_count %>%
    mutate(true_signal_check = ifelse(D0 != 0, T, F)) %>%
    mutate(category_1_check = ifelse(asv_id %in% value_category_1_asv_name, T, F)) %>%
    mutate(category_2_check = ifelse(asv_id %in% value_category_2_asv_name, T, F)) %>%
    mutate(category_3_check = ifelse(!(asv_id %in% value_merged_ntc_positive_asv_name), T, F)) %>%
    mutate(D1_check = (D1_SCRuB != 0) & (D1_Green_Cleaner == 0)) %>%
    mutate(D2_check = (D2_SCRuB != 0) & (D2_Green_Cleaner == 0)) %>%
    mutate(D3_check = (D3_SCRuB != 0) & (D3_Green_Cleaner == 0)) %>%
    mutate(D4_check = (D4_SCRuB != 0) & (D4_Green_Cleaner == 0)) %>%
    mutate(D5_check = (D5_SCRuB != 0) & (D5_Green_Cleaner == 0)) %>%
    mutate(D6_check = (D6_SCRuB != 0) & (D6_Green_Cleaner == 0)) %>%
    mutate(D7_check = (D7_SCRuB != 0) & (D7_Green_Cleaner == 0))
  
  data_inter_contam_category_check <-
    data_tmp_category_check %>%
    select(asv_id, true_signal_check, category_1_check, category_2_check, category_3_check, D1_check, D2_check, D3_check, D4_check, D5_check, D6_check, D7_check) %>%
    filter(true_signal_check == F) %>%
    mutate(category = apply(., 1, func.category.check))
  
  data_category_count <- 
    data_inter_contam_category_check %>%
    filter(D1_check == T) %>%
    select(category) %>%
    mutate(type = "D1") %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D2_check == T) %>%
            select(category) %>%
            mutate(type = "D2")) %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D3_check == T) %>%
            select(category) %>%
            mutate(type = "D3")) %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D4_check == T) %>%
            select(category) %>%
            mutate(type = "D4")) %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D5_check == T) %>%
            select(category) %>%
            mutate(type = "D5")) %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D6_check == T) %>%
            select(category) %>%
            mutate(type = "D6")) %>%
    rbind(., data_inter_contam_category_check %>%
            filter(D7_check == T) %>%
            select(category) %>%
            mutate(type = "D7"))
  
  for (value_d_stage in unique(data_category_count$type)) {
    
    # value_d_stage <- data_category_count$type[1]
    
    data_target_table <- 
      data_category_count %>%
      filter(type == value_d_stage)
    
    value_total <- nrow(data_target_table)
    
    value_category_1 <- sum(data_target_table$category == "category_1")
    value_category_2 <- sum(data_target_table$category == "category_2")
    value_category_3 <- sum(data_target_table$category == "category_3")
    
    data_tmp_category_relab <- 
      data.frame(dataset = rep(value_tmp_dataset, 3),
                 d_stage = rep(value_d_stage, 3),
                 category = glue("category_{seq(3)}"),
                 value = c(value_category_1 / value_total, value_category_2 / value_total, value_category_3 / value_total))
    
    if (value_d_stage == unique(data_category_count$type)[1]) {
      
      data_merged_category_relab <- data_tmp_category_relab
      
    } else {
      
      data_merged_category_relab <- rbind(data_merged_category_relab, data_tmp_category_relab)
      
    }
    
  }
  
  data_category_count$dataset <- value_tmp_dataset
  
  if (value_tmp_dataset == input_dataset_list[1]) {
    
    data_dataset_merged_category_relab <- data_merged_category_relab
    data_merged_category_count <- data_category_count
    
  } else {
    
    data_dataset_merged_category_relab <- 
      rbind(data_dataset_merged_category_relab, data_merged_category_relab)
    data_merged_category_count <- 
      rbind(data_merged_category_count, data_category_count)
    
  }
  
}

for (value_d_stage in unique(data_merged_category_count$type)) {
  
  # value_d_stage <- data_category_count$type[1]
  
  data_target_table <- 
    data_merged_category_count %>%
    filter(type == value_d_stage)
  
  value_total <- nrow(data_target_table)
  
  value_category_1 <- sum(data_target_table$category == "category_1")
  value_category_2 <- sum(data_target_table$category == "category_2")
  value_category_3 <- sum(data_target_table$category == "category_3")
  
  data_tmp_category_relab <- 
    data.frame(d_stage = rep(value_d_stage, 3),
               category = glue("category_{seq(3)}"),
               value = c(value_category_1 / value_total, value_category_2 / value_total, value_category_3 / value_total))
  
  if (value_d_stage == unique(data_category_count$type)[1]) {
    
    data_merged_all_dataset_category_relab <- data_tmp_category_relab
    
  } else {
    
    data_merged_all_dataset_category_relab <- rbind(data_merged_all_dataset_category_relab, data_tmp_category_relab)
    
  }
  
}

data_merged_all_dataset_category_relab <-
  data_merged_all_dataset_category_relab %>%
  rename("Dilution series" = "d_stage") %>%
  rename("Category" = "category") %>%
  rename("Proportion" = "value") %>%
  mutate(Category = gsub("_", " ", Category)) %>%
  mutate(`Dilution series` = factor(`Dilution series`, levels = glue("D{seq(7)}")))

fig_category <- 
  ggplot(data_merged_all_dataset_category_relab, aes(fill = Category, x = `Dilution series`, y = Proportion)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()

png(glue("{path_output}/additional_file5_category_stacked_barplot_figure.png"), width = 500, height = 300)

print(fig_category)

dev.off()

write.table(data_merged_all_dataset_category_relab, glue("{path_output}/additional_file5_category_stacked_barplot_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)