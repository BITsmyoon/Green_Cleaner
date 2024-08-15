# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))

# distance
suppressPackageStartupMessages(library(proxy))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <- 
  glue("{path_data}/../output/batch_decontaminated")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

input_meta_data <- 
  fread(glue("{path_data}/input_sample_meta_data.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_asv_count <- 
  fread(glue("{path_data}/input_raw_asv_count.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_blacklist <- 
  fread(glue("{path_data}/Green_Cleaner_in_house_blacklist.txt"), 
        sep = "\t", quote = F, header = F, stringsAsFactors = F, data.table = F) %>%
  pull("V1")
input_bacdive <- 
  fread(glue("{path_data}/BacDive_not_human_source.txt"), 
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_taxa_blacklist <- 
  fread(glue("{path_data}/Green_Cleaner_taxa_blacklist.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_taxa_table <- 
  fread(glue("{path_data}/input_taxa_table.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

make.df <- function(x, y, z) {
  
  sapply(strsplit(x = x, split = y), '[[', z)
  
}


func.dist <- 
  function(value_tmp) {
    sqrt((as.numeric(value_tmp[5]) - as.numeric(value_tmp[2]))^2 + (as.numeric(value_tmp[6]) - as.numeric(value_tmp[3]))^2)
  }

func.make.bacdive.key.value <- function(value) {
  
  value_taxa_name <- value[2]
  value_type <- value[3]
  
  if (value_type == "genus") {
    
    value_key <- make.df(value_taxa_name, ";g_", 2)
    
  } else if (value_type == "species") {
    
    value_key <- 
      make.df(value_taxa_name, ";g_", 2) %>%
      make.df(., ";s_", 1)
    
  }
  
  return(value_key)
  
}

func.make.source.output <- 
  function(value) {
    
    value_type <- value[3]
    value_key <- value[4]
    
    value_output <- ifelse(value_key %in% input_bacdive$modi_genus, "not_human_source", "human_source")
    
    return(value_output)
    
  }

func.exclude.run <- function(value) {
  
  data_tmp_exclude <- 
    data_exclude_asv %>%
    filter(case == value)
  
  for (value_num in seq(nrow(data_tmp_exclude))) {
    
    value_tmp_sample <- data_tmp_exclude$sample[value_num]
    value_tmp_asv <- data_tmp_exclude$asv[value_num]
    
    data_final_asv[value_tmp_asv, value_tmp_sample] <<- 0
    
  }
  
}

# 2. pre processing ------------------------------------------------------------

for (input_decontamination_batch in glue("batch_{seq(6)}")) {
 
  data_selected_batch_meta_data <-
    input_meta_data %>%
    filter(decontamination_batch == input_decontamination_batch)
  
  value_green_cleaner_input_sample_id <- 
    data_selected_batch_meta_data %>%
    filter(dilution_series != "d0") %>%
    pull("sample_id")
  
  data_selected_asv_count <- 
    input_asv_count %>%
    select(asv_id, all_of(value_green_cleaner_input_sample_id))
  
  logic_zero_asv <- 
    data_selected_asv_count %>%
    column_to_rownames("asv_id") %>%
    rowSums() != 0
  
  data_filtered_asv_count <- 
    data_selected_asv_count[logic_zero_asv, ] %>%
    remove_rownames() %>%
    column_to_rownames("asv_id")
  
  value_ntc_name <- 
    data_selected_batch_meta_data %>%
    filter(ntc_check == "ntc") %>%
    pull("sample_id")
  
  data_filtered_asv_relab <- 
    data_filtered_asv_count %>%
    apply(., 2, function(value){
      (value / sum(value)) * 100
    }) %>%
    as.data.frame()
  
  value_ntc_negative_asv_name <- 
    rownames(data_filtered_asv_count)[data_filtered_asv_count[, value_ntc_name] == 0]
  value_ntc_positive_asv_name <- 
    rownames(data_filtered_asv_count)[data_filtered_asv_count[, value_ntc_name] != 0]
  
  data_ntc_asv <-
    data_filtered_asv_count %>%
    select(all_of(value_ntc_name))
  
  value_ntc_top5_asv <- 
    data_ntc_asv %>%
    arrange(desc(.[, 1])) %>%
    head(., n = 5) %>%
    rownames()
  
  data_ntc_top5_sample_relab <- 
    data_filtered_asv_relab[value_ntc_top5_asv, ]
  data_ntc_top5_sample_asv <- 
    data_filtered_asv_count[value_ntc_top5_asv, ]
  
  data_final_asv <- data_filtered_asv_count
  
  data_include_asv <- data.frame(sample = NULL,
                                 asv = NULL,
                                 case = NULL)
  data_exclude_asv <- data.frame(sample = NULL,
                                 asv = NULL,
                                 case = NULL)
  
  value_phylum_taxa_blacklist <- 
    input_taxa_blacklist %>%
    filter(level == "p") %>%
    pull("name")
  value_class_taxa_blacklist <- 
    input_taxa_blacklist %>%
    filter(level == "c") %>%
    pull("name")
  
  # 3. filtering -----------------------------------------------------------------
  # ... 3.1. case 2 --------------------------------------------------------------
  
  data_ntc_top5_colsums <- 
    colSums(data_ntc_top5_sample_relab) %>%
    t() %>%
    as.data.frame()
  
  value_case_2_sample <- 
    data_ntc_top5_colsums[, data_ntc_top5_colsums == 0] %>%
    colnames()
  
  if (length(value_case_2_sample) != 0) {
    
    for (value_sample_name in value_case_2_sample) {
      
      value_tmp_include_asv <- 
        data_filtered_asv_count %>%
        select(all_of(`value_sample_name`)) %>%
        filter(.[, 1] != 0) %>%
        rownames()
      
      data_include_asv <- 
        rbind(data_include_asv, data.frame(sample = value_sample_name,
                                           asv = value_tmp_include_asv,
                                           case = "case_2"))
      
    }
    
  }
  
  # ... 3.2. case 3 -----------------------------------------------------------
  
  value_case_3_sample <- 
    data_ntc_top5_colsums[, (data_ntc_top5_colsums < 5) & (data_ntc_top5_colsums != 0)] %>%
    colnames()
  
  if (length(value_case_3_sample) != 0) {
    
    for (value_sample_name in value_case_3_sample) {
      
      value_tmp_exclude_asv <- 
        c(value_ntc_top5_asv, data_filtered_asv_relab %>%
            select(all_of(`value_sample_name`)) %>%
            filter(.[, 1] != 0 & .[, 1] < .5) %>%
            rownames()) %>%
        unique()
      
      data_exclude_asv <-
        rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                           asv = value_tmp_exclude_asv,
                                           case = "case_3"))
      
    }
    
  }
  
  # ... 3.3. case 4-1 -----------------------------------------------------------
  
  value_ntc_top5_sum_upper_5_sample_name <- 
    data_ntc_top5_colsums[, data_ntc_top5_colsums >= 5] %>%
    colnames()
  
  if (!is.null(value_ntc_top5_sum_upper_5_sample_name)) {
    
    data_case_4_1_input_asv <- 
      data_ntc_top5_sample_asv[, value_ntc_top5_sum_upper_5_sample_name]
    
    data_4_1_input_relab <- 
      data_case_4_1_input_asv %>%
      apply(., 2, function(value){
        (value / sum(value)) * 100
      }) %>%
      as.data.frame()
    
    data_matrix <- 
      as.matrix(data_4_1_input_relab) %>%
      t()
    
    pca_result <- 
      prcomp(data_matrix, scale. = F)
    
    data_pca_euclidean <- 
      proxy::simil(pca_result$x[, 1:2], method = "euclidean") %>%
      as.matrix() %>%
      as.data.frame() %>%
      select(all_of(value_ntc_name))
    
    value_distance_cutoff_under_sample_name <- 
      rownames(data_pca_euclidean)[!is.na(data_pca_euclidean[, 1]) & data_pca_euclidean[, 1] < .019]
    
    colnames(data_pca_euclidean) <- "pca_euclidean"
    
    if (length(value_distance_cutoff_under_sample_name) == 0) {
      
      for (value_sample_name in colnames(data_4_1_input_relab)) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_ntc_top5_asv,
                                             case = "case_4_1"))
        
      }
      
    }
    
    if (length(value_distance_cutoff_under_sample_name) != 0) {
      
      data_arrow_point <- 
        data.frame(asv = colnames(data_matrix),
                   arrow_PC1 = apply(data_matrix, 2, function(value){
                     cor(value, pca_result$x[,1]) * 0.8 * sqrt(nrow(data_matrix) - 1)
                   }),
                   arrow_PC2 = apply(data_matrix, 2, function(value){
                     cor(value, pca_result$x[,2]) * 0.8 * sqrt(nrow(data_matrix) - 1)
                   }))
      
      data_pca_point <- 
        pca_result$x %>%
        as.data.frame() %>%
        select(PC1, PC2) %>%
        rename("data_PC1" = "PC1") %>%
        rename("data_PC2" = "PC2") %>%
        rownames_to_column("sample")
      
      data_sample_arrow_point_dist <- 
        merge(data_pca_point, data_arrow_point, by = NULL, all = TRUE) %>%
        mutate(dist = apply(., 1, func.dist))
      
      data_selected_sample_arrow_point_dist <- 
        data_sample_arrow_point_dist %>%
        filter(sample %in% value_distance_cutoff_under_sample_name)
      
      data_dist_min_asv <- data_selected_sample_arrow_point_dist %>%
        group_by(sample) %>%
        slice(which.min(dist))
      
      data_dist_min_asv_all_info <- 
        merge(data_pca_euclidean %>%
                rownames_to_column("sample") %>%
                filter(sample %in% value_distance_cutoff_under_sample_name), 
              data_dist_min_asv %>%
                select(sample, asv, dist), by = "sample")
      
      for (value_sample_name in data_dist_min_asv_all_info$sample) {
        
        value_tmp_include_asv <- 
          data_dist_min_asv_all_info %>%
          filter(sample == value_sample_name) %>%
          pull("asv")
        
        data_include_asv <-
          rbind(data_include_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_include_asv,
                                             case = "case_4_1"))
        
      }
      
      data_tmp_exclude <-
        data.frame(sample = rep(colnames(data_case_4_1_input_asv), each = 5),
                   asv = rep(rownames(data_case_4_1_input_asv), length(colnames(data_case_4_1_input_asv))),
                   case = "case_4_1")
      
      data_exclude_asv <-
        rbind(data_exclude_asv, anti_join(data_tmp_exclude, data_include_asv, by = c("sample", "asv", "case")))
      
    }
    
    # ... 3.4. case 4-2 -----------------------------------------------------------
    
    value_case_4_2_raw_asv_name <- 
      value_ntc_positive_asv_name[!(value_ntc_positive_asv_name %in% value_ntc_top5_asv)]
    
    data_case_4_2_raw_relab <- 
      data_filtered_asv_relab[value_case_4_2_raw_asv_name, value_ntc_top5_sum_upper_5_sample_name]
    
    value_case_4_2_asv_name <- 
      names(rowSums(data_case_4_2_raw_relab != 0))[rowSums(data_case_4_2_raw_relab != 0) >= 4]
    
    data_case_4_2_relab <- 
      data_case_4_2_raw_relab[value_case_4_2_asv_name, ]
    
    data_case_4_2_norm_relab <- data_case_4_2_relab
    
    if (length(data_case_4_2_norm_relab) >= 4) {
      
      for (value_tmp_sample in colnames(data_case_4_2_norm_relab)) {
        
        data_case_4_2_norm_relab[, value_tmp_sample] <- 
          data_case_4_2_norm_relab[, value_tmp_sample] / (data_ntc_top5_colsums[, value_tmp_sample])
        
      }
      
      data_z_score <- 
        apply(data_case_4_2_norm_relab
              , 1, function(value){
                value_not_zero <- value[value != 0]
                
                value_median_value <- median(value_not_zero)
                value_mad <- median(abs(value_not_zero - value_median_value))
                value_mean_ad <- mean(abs(value_not_zero - mean(value_median_value)))
                
                if (value_mad != 0) {
                  value_modi_z_score_denominator <- value_mad * 1.4826
                } else if (value_mad == 0) {
                  value_modi_z_score_denominator <- value_mean_ad * 1.2533
                }
                round((value - median(value_not_zero)) / value_modi_z_score_denominator, 4)
              }) %>%
        t() %>%
        as.data.frame()
      
      for (value_sample_name in colnames(data_z_score)) {
        
        data_tmp_z_score <-
          data_z_score %>%
          select(all_of(`value_sample_name`))
        
        data_tmp_relab <- 
          data_case_4_2_relab %>%
          select(all_of(`value_sample_name`))
        
        value_tmp_exclude_asv <-
          rownames(data_tmp_relab)[data_tmp_relab != 0]
        
        value_tmp_include_asv <- 
          data_tmp_z_score %>%
          filter(.[, 1] >= 8) %>%
          rownames()
        
        if (length(value_tmp_include_asv) != 0) {
          
          data_include_asv <-
            rbind(data_include_asv, data.frame(sample = value_sample_name,
                                               asv = value_tmp_include_asv,
                                               case = "case_4_2"))
          
          value_tmp_exclude_asv <-
            value_tmp_exclude_asv[!(value_tmp_exclude_asv %in% value_tmp_include_asv)]
          
          if (length(value_tmp_exclude_asv) != 0) {
            
            data_exclude_asv <-
              rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                                 asv = value_tmp_exclude_asv,
                                                 case = "case_4_2"))
            
          }
          
        } else if (length(value_tmp_exclude_asv) != 0) {
          
          data_exclude_asv <-
            rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                               asv = value_tmp_exclude_asv,
                                               case = "case_4_2"))
          
        }
        
      }
      
    }
    
    # ... 3.5. case 4-3 -----------------------------------------------------------
    
    value_case_4_3_raw_asv_name <- 
      value_ntc_positive_asv_name[!(value_ntc_positive_asv_name %in% value_ntc_top5_asv)]
    
    data_case_4_3_raw_relab <- 
      data_filtered_asv_relab[value_case_4_3_raw_asv_name, value_ntc_top5_sum_upper_5_sample_name]
    
    value_case_4_3_asv_name <- 
      names(rowSums(data_case_4_3_raw_relab != 0))[rowSums(data_case_4_3_raw_relab != 0) < 4]
    
    data_case_4_3_relab <- 
      data_case_4_3_raw_relab[value_case_4_3_asv_name, ]
    
    data_in_genus_taxa_table <-
      input_taxa_table[!grepl(pattern = ";g_;", x = input_taxa_table$Taxon), ] %>%
      mutate(Taxon = gsub(pattern = ";s_$", replacement = "", Taxon)) %>%
      select(asv_id, Taxon) %>%
      mutate(type = ifelse(grepl(pattern = ";s_", Taxon), "species", "genus"))
    
    data_in_genus_taxa_table$key_value <-
      apply(data_in_genus_taxa_table, 1, func.make.bacdive.key.value)
    data_in_genus_taxa_table$source <-
      apply(data_in_genus_taxa_table, 1, func.make.source.output)
    
    value_not_human_asv_list <-
      data_in_genus_taxa_table %>%
      filter(source == "not_human_source") %>%
      pull("asv_id")
    
    data_taxa_blacklist <- 
      data_filtered_asv_count %>%
      rownames_to_column("asv") %>%
      select(asv) %>%
      merge(., input_taxa_table %>%
              select(asv_id, Taxon) %>%
              rename("asv" = "asv_id") %>%
              rename("taxa" = "Taxon"), by = "asv") %>%
      mutate(phylum_name = make.df(taxa, ";p_", 2)) %>%
      mutate(phylum_name = make.df(phylum_name, ";c_", 1)) %>%
      mutate(class_name = make.df(taxa, ";c_", 2)) %>%
      mutate(class_name = make.df(class_name, ";o_", 1)) %>%
      mutate(order_name = make.df(taxa, ";o_", 2)) %>%
      mutate(order_name = make.df(order_name, ";f_", 1)) %>%
      mutate(taxa_blacklist = (phylum_name %in% value_phylum_taxa_blacklist) | (class_name %in% value_class_taxa_blacklist) | (phylum_name == "Proteobacteria" & class_name == "") | (class_name == "Gammaproteobacteria" & order_name == ""))
    
    value_taxa_blacklist_asv_name <-
      data_taxa_blacklist %>%
      filter(taxa_blacklist) %>%
      pull("asv")
    
    for (value_sample_name in colnames(data_case_4_3_relab)) {
      
      data_tmp_asv <- 
        data_case_4_3_relab %>%
        select(all_of(value_sample_name)) %>%
        rownames_to_column("asv")
      
      value_tmp_bacdive_input_asv_name <- 
        data_tmp_asv[(data_tmp_asv[, 2] != 0), ] %>%
        pull("asv")
      
      logic_asv_filter <- 
        (value_tmp_bacdive_input_asv_name %in% value_not_human_asv_list) | (value_tmp_bacdive_input_asv_name %in% value_taxa_blacklist_asv_name)
      
      value_tmp_exclude_asv_name <- 
        value_tmp_bacdive_input_asv_name[logic_asv_filter]
      
      if (length(value_tmp_exclude_asv_name) != 0) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_exclude_asv_name,
                                             case = "case_4_3"))
        
      }
      
    }
    
    # ... 3.6. case 4-4 -----------------------------------------------------------
    
    data_case_4_4_raw_relab <- 
      data_filtered_asv_relab[value_ntc_negative_asv_name, value_ntc_top5_sum_upper_5_sample_name]
    
    for (value_sample_name in colnames(data_case_4_4_raw_relab)) {
      
      data_tmp_asv <- 
        data_case_4_4_raw_relab %>%
        select(all_of(value_sample_name)) %>%
        rownames_to_column("asv")
      
      value_tmp_bacdive_input_asv_name <- 
        data_tmp_asv[(data_tmp_asv[, 2] >= 5), ] %>%
        pull("asv")
      
      logic_asv_filter <- 
        (value_tmp_bacdive_input_asv_name %in% value_not_human_asv_list) | (value_tmp_bacdive_input_asv_name %in% value_taxa_blacklist_asv_name)
      
      value_tmp_exclude_asv_name <- 
        value_tmp_bacdive_input_asv_name[logic_asv_filter]
      
      if (length(value_tmp_exclude_asv_name) != 0) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_exclude_asv_name,
                                             case = "case_4_4"))
        
      }
      
    }
    
    # ... 3.7. case 4-5 -----------------------------------------------------------
    
    data_case_4_5_raw_relab <- 
      data_filtered_asv_relab[value_ntc_negative_asv_name, value_ntc_top5_sum_upper_5_sample_name]
    
    for (value_sample_name in colnames(data_case_4_5_raw_relab)) {
      
      data_tmp_asv <-
        data_case_4_5_raw_relab %>%
        select(all_of(value_sample_name)) %>%
        rownames_to_column("asv")
      
      value_tmp_bacdive_input_asv_name <-
        data_tmp_asv[(data_tmp_asv[, 2] != 0) & (data_tmp_asv[, 2] < 5), ] %>%
        pull("asv")
      
      logic_asv_filter <-
        (value_tmp_bacdive_input_asv_name %in% value_not_human_asv_list) | (value_tmp_bacdive_input_asv_name %in% value_taxa_blacklist_asv_name) | (value_tmp_bacdive_input_asv_name %in% input_blacklist)
      
      value_tmp_exclude_asv_name <-
        value_tmp_bacdive_input_asv_name[logic_asv_filter]
      
      if (length(value_tmp_exclude_asv_name) != 0) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_exclude_asv_name,
                                             case = "case_4_5"))
        
      }
      
    }
    
    # ... 3.8 case 4-6 -----------------------------------------------------------
    
    data_case_4_6_raw_relab <- 
      data_filtered_asv_relab[value_ntc_negative_asv_name, value_ntc_top5_sum_upper_5_sample_name]
    
    data_exclude_case_4_5_asv <- 
      data_exclude_asv %>%
      filter(case == "case_4_5")
    
    for (value_sample_name in colnames(data_case_4_6_raw_relab)) {
      
      value_tmp_exclude_case_4_5_asv_name <- 
        data_exclude_case_4_5_asv %>%
        filter(sample == value_sample_name) %>%
        pull("asv")
      
      data_tmp_asv <-
        data_case_4_6_raw_relab %>%
        select(all_of(value_sample_name)) %>%
        rownames_to_column("asv")
      
      value_tmp_exclude_asv_name <-
        data_tmp_asv[(data_tmp_asv[, 2] != 0) & (data_tmp_asv[, 2] < .1), ] %>%
        pull("asv")
      
      value_tmp_exclude_asv_name <-
        value_tmp_exclude_asv_name[!(value_tmp_exclude_asv_name %in% value_tmp_exclude_case_4_5_asv_name)]
      
      if (length(value_tmp_exclude_asv_name) != 0) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_exclude_asv_name,
                                             case = "case_4_6"))
        
      }
      
    }
    
    # ... 3.9 case 4-7 -----------------------------------------------------------
    
    for (value_sample_name in value_ntc_top5_sum_upper_5_sample_name) {
      
      data_tmp_asv <-
        data_filtered_asv_count %>%
        select(all_of(value_sample_name)) %>%
        rownames_to_column("asv")
      
      value_tmp_exclude_asv_name <-
        data_tmp_asv[(data_tmp_asv[, 2] != 0) & (data_tmp_asv[, 2] < 10), ] %>%
        pull("asv")
      
      if (length(value_tmp_exclude_asv_name) != 0) {
        
        data_exclude_asv <-
          rbind(data_exclude_asv, data.frame(sample = value_sample_name,
                                             asv = value_tmp_exclude_asv_name,
                                             case = "case_4_7"))
        
      }
      
    }
    
  } else {
    
    print("no ntc top5 sum >= 5 sample")
    
  }
  
  # ... 3.10 case taxa filtering -----------------------------------------------------------
  
  data_taxa_filtering_asv <-
    merge(data_filtered_asv_count %>%
            rownames_to_column("asv"),
          input_taxa_table %>%
            select(asv_id, Taxon) %>%
            rename("asv" = "asv_id"),
          by = "asv")
  
  value_taxa_filtering_asv <-
    data_taxa_filtering_asv$asv[grepl(pattern = ";o_;f_;g_;s_$", x = data_taxa_filtering_asv$Taxon)]
  
  value_all_sample_name <- colnames(data_filtered_asv_count)[!(colnames(data_filtered_asv_count) %in% value_ntc_name)]
  
  for (value_case in names(table(data_exclude_asv$case))) {
    
    func.exclude.run(value_case)
    
  }
  
  data_final_asv[rownames(data_final_asv) %in% value_taxa_filtering_asv, ] <- 0
  
  # 4. save ----------------------------------------------------------------------
  
  data_tmp <- data_final_asv[, !(colnames(data_final_asv) %in% c(value_ntc_name))]
  
  if (is.null(dim(data_tmp))) {
    
    data_filter_out_asv <- 
      data.frame(otu = rownames(data_final_asv),
                 tmp = data_tmp)
    
    colnames(data_filter_out_asv) <- c("#OTU ID", 
                                       colnames(data_final_asv)[!(colnames(data_final_asv) %in% c(value_ntc_name))])
    
  } else {
    
    data_filter_out_asv <-
      data_final_asv[, !(colnames(data_final_asv) %in% c(value_ntc_name))] %>%
      rownames_to_column("#OTU ID")
    
  }
  
  logic_filtered_asv <- 
    data_filter_out_asv %>%
    column_to_rownames("#OTU ID") %>%
    rowSums() != 0
  
  data_decontam_asv_count <- 
    data_filter_out_asv[logic_filtered_asv, ] %>%
    rename("asv_id" = "#OTU ID")
  
  write.table(data_decontam_asv_count, glue("{path_output}/{input_decontamination_batch}_Green_Cleaner_decontaminated_asv_count.txt"),
              quote = F, sep = "\t", row.names = F, col.names = T)
   
}