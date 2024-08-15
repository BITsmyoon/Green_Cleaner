# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rlang))

# Plotting
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Polychrome))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_output <-
  glue("{path_data}/../output/paper_output/stacked_barplot")

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_species_count <-
  read.csv(glue("{path_data}/../output/qiime/qiime_barplot/level-7.csv"), 
           header = T, sep = ",", check.names = F, stringsAsFactors = F) %>%
  select(-dataset)
input_taxa_info <-
  fread(glue("{path_data}/input_rev_stacked_barplot_color_list.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_dataset_list <- glue("dataset_{seq(10)}")

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

make.df <- function(x, y, z) {
  
  sapply(strsplit(x = x, split = y), '[[', z)
  
}

# 2. pre processing ------------------------------------------------------------

value_color_taxa <- 
  input_taxa_info %>%
  filter(d0_in == T) %>%
  pull("taxa") %>%
  unique() 

value_tmp2_color <- palette36.colors(n = length(value_color_taxa) + 2)
value_color <- value_tmp2_color[-c(1, 2)]

data_taxa_color_list <- 
  data.frame(taxa = value_color_taxa,
             color = as.character(value_color))

input_dataset_list <- glue("dataset_{seq(10)}")

for (value_tmp_dataset in input_dataset_list) {
  
  # value_tmp_dataset <- input_dataset_list[1]
  
  value_tmp_ntc_name <- 
    input_meta %>%
    filter(ntc_check == "ntc") %>%
    filter(grepl(pattern = value_tmp_dataset, x = dataset)) %>%
    pull("sample_id")
  
  value_tmp_target_sample <- 
    input_meta %>%
    filter(dataset == value_tmp_dataset) %>%
    filter(dilution_series != "-") %>%
    pull("sample_id")
  
  if (length(value_tmp_ntc_name) >= 2) {
    
    value_tmp_ntc_name <- 
      value_tmp_ntc_name[1]
    
  }
  
  value_ntc_target_taxa <- 
    input_species_count %>%
    filter(index == value_tmp_ntc_name) %>%
    column_to_rownames("index") %>%
    t() %>%
    as.data.frame() %>%
    rename("ntc" = 1) %>%
    arrange(desc(ntc)) %>%
    head(10) %>%
    rownames() %>%
    gsub(pattern = ";s_$", replacement = "", .) %>%
    gsub(pattern = ";g_$", replacement = "", .) %>%
    gsub(pattern = ";f_$", replacement = "", .) %>%
    gsub(pattern = ";o_$", replacement = "", .) %>%
    gsub(pattern = ";c_$", replacement = "", .) %>%
    gsub(pattern = ";p_$", replacement = "", .)
  
  value_tmp_len <- 
    strsplit(x = value_ntc_target_taxa, split = ";") %>%
    sapply(., length)
  
  for (value_for_num in seq(length(value_tmp_len))) {
    
    # value_for_num <- 2
    
    value_target_len <- value_tmp_len[value_for_num]
    value_target_taxa <- value_ntc_target_taxa[value_for_num]
    
    if (value_target_len != 7) {
      
      value_tmp_taxa_name <- make.df(value_target_taxa, ";", value_target_len)
      value_tmp_taxa_name <- paste(toupper(substr(value_tmp_taxa_name, 1, 1)), substr(value_tmp_taxa_name, 2, nchar(value_tmp_taxa_name)), sep = "")
      
    } else {
      
      value_tmp_taxa_name <-
        glue('S_{make.df(value_target_taxa, ";", 6) %>%
      make.df(., "g_", 2)}_{make.df(value_target_taxa, ";", 7) %>%
      make.df(., "s_", 2)}')
      
    }
    
    if (value_for_num == 1) {
      
      value_processed_ntc_target_taxa_name <- value_tmp_taxa_name
      
    } else {
      
      value_processed_ntc_target_taxa_name <- c(value_processed_ntc_target_taxa_name, value_tmp_taxa_name)
      
    }
    
  }
  
  data_tmp <- 
    input_species_count %>%
    filter(index == value_tmp_ntc_name) %>%
    column_to_rownames("index") %>%
    t() %>%
    as.data.frame() %>%
    rename("ntc" = 1) %>%
    arrange(desc(ntc)) %>%
    apply(., 2, function(value){
      (value / sum(value)) * 100
    }) %>%
    as.data.frame() %>%
    head(10) %>%
    mutate(check = ifelse(ntc >= 5, T, F))
  
  data_tmp_taxa_count <- 
    input_species_count %>%
    filter(index %in% value_tmp_target_sample) %>%
    column_to_rownames("index")
  
  data_tmp_taxa_relab <-
    data_tmp_taxa_count[, colSums(data_tmp_taxa_count) != 0] %>%
    apply(., 1, function(value){
      (value / sum(value)) * 100
    }) %>%
    as.data.frame()
  
  colnames(data_tmp_taxa_relab) <-
    glue("D{seq(from = 0, to = 7)}")
  
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";s_$", replacement = "", rownames(data_tmp_taxa_relab))
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";g_$", replacement = "", rownames(data_tmp_taxa_relab))
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";f_$", replacement = "", rownames(data_tmp_taxa_relab))
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";o_$", replacement = "", rownames(data_tmp_taxa_relab))
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";c_$", replacement = "", rownames(data_tmp_taxa_relab))
  rownames(data_tmp_taxa_relab) <- gsub(pattern = ";p_$", replacement = "", rownames(data_tmp_taxa_relab))
  
  value_tmp_len <- 
    strsplit(x = rownames(data_tmp_taxa_relab), split = ";") %>%
    sapply(., length)
  
  for (value_for_num in seq(nrow(data_tmp_taxa_relab))) {
    
    # value_for_num <- 2
    # value_for_num <- 4
    
    value_target_len <- value_tmp_len[value_for_num]
    value_target_taxa <- rownames(data_tmp_taxa_relab)[value_for_num]
    
    if (value_target_len != 7) {
      
      value_tmp_taxa_name <- make.df(value_target_taxa, ";", value_target_len)
      value_tmp_taxa_name <- paste(toupper(substr(value_tmp_taxa_name, 1, 1)), substr(value_tmp_taxa_name, 2, nchar(value_tmp_taxa_name)), sep = "")
      
    } else {
      
      value_tmp_taxa_name <-
        glue('S_{make.df(value_target_taxa, ";", 6) %>%
      make.df(., "g_", 2)}_{make.df(value_target_taxa, ";", 7) %>%
      make.df(., "s_", 2)}')
      
    }
    
    rownames(data_tmp_taxa_relab)[value_for_num] <- value_tmp_taxa_name
    
  }
  
  value_empty_taxa <- 
    value_processed_ntc_target_taxa_name[!(value_processed_ntc_target_taxa_name %in% rownames(data_tmp_taxa_relab))]
  
  if (length(value_empty_taxa) != 0) {
    
    data_tmp_taxa_relab[value_empty_taxa, ] <- 0
    
  }
  
  data_arranged_taxa_relab <- 
    data_tmp_taxa_relab %>%
    arrange(desc(D0))
  
  data_taxa_info <- 
    data.frame(taxa = rownames(data_arranged_taxa_relab),
               d0_in = rownames(data_arranged_taxa_relab) %in% value_processed_ntc_target_taxa_name,
               d0_5_up = rownames(data_arranged_taxa_relab) %in% value_processed_ntc_target_taxa_name[data_tmp$check]) %>%
    arrange(desc(d0_in))
  
  value_viridis_palette <- rainbow(sum(data_taxa_info$d0_in)) %>%
    rev()
  value_gray_palette <- c("#6C6C6C", "#8D8D8D", "#A7A7A7", "#BDBDBD")
  
  value_nrow <- nrow(data_taxa_info)
  value_color_count <- sum(data_taxa_info$d0_in)
  value_rep_gray <- (value_nrow - value_color_count) %/% 4
  value_other <- value_nrow - (value_color_count + (value_rep_gray * 4))
  
  data_tmp2 <- 
    merge(data_taxa_info %>%
            filter(d0_in) %>%
            select(taxa),
          data_taxa_color_list,
          by = "taxa") %>%
    column_to_rownames("taxa")
  
  value_color_list <- 
    data_tmp2[data_taxa_info %>%
                filter(d0_in) %>%
                pull("taxa"), ]
  
  if (value_other != 0) {
    
    data_taxa_all_info <- 
      data_taxa_info %>%
      mutate(color = c(value_color_list,
                       rep(value_gray_palette, value_rep_gray),
                       value_gray_palette[seq(value_other)]))
    
  } else {
    
    data_taxa_all_info <- 
      data_taxa_info %>%
      mutate(color = c(value_color_list,
                       rep(value_gray_palette, value_rep_gray)))
    
  }
  
  data_arranged_taxa_relab <-
    data_arranged_taxa_relab[data_taxa_all_info$taxa, ]
  
  for (tmp_count in seq(ncol(data_arranged_taxa_relab))) {
    
    # tmp_count <- 1
    
    if (tmp_count == 1) {
      
      value_total_relab <- data_arranged_taxa_relab[, tmp_count]
      
    } else {
      
      value_total_relab <- c(value_total_relab, data_arranged_taxa_relab[, tmp_count])
      
    }
    
  }
  
  data_plot_input <- 
    data.frame(taxa = rep(rownames(data_arranged_taxa_relab), ncol(data_arranged_taxa_relab)),
               d_series = rep(colnames(data_arranged_taxa_relab), each = nrow(data_arranged_taxa_relab)),
               value = value_total_relab) %>%
    mutate(batch = value_tmp_dataset) %>%
    mutate(taxa = factor(taxa, levels = rownames(data_arranged_taxa_relab)))
  
  data_legend <-
    data_taxa_all_info %>%
    filter(d0_in) %>%
    select(taxa, color)
  
  value_selected_legend <-
    data_legend %>%
    select(color) %>%
    pull() %>%
    set_names(data_legend$taxa)
  
  value_legend_vector <- setNames(data_taxa_all_info$color, data_taxa_all_info$taxa)
  
  data_figure <- 
    ggplot(data_plot_input, aes(fill = taxa, y = value, x = d_series)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual("Bacteria", 
                      values = value_legend_vector, 
                      labels = data_legend$taxa,
                      breaks = data_legend$taxa) +
    theme_bw() + 
    labs(y = "Relative abundance", x = "Dilution series and Type") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
  png(glue("{path_output}/{value_tmp_dataset}_rev_barplot.png"), width = 1200, height = 600)
  
  print(data_figure)
  
  dev.off()
  
  if (value_tmp_dataset == input_dataset_list[1]) {
    
    data_merged_plot_input <- data_plot_input
    data_merged_taxa_info <- data_taxa_info
    
  } else {
    
    data_merged_plot_input <- rbind(data_merged_plot_input, data_plot_input)
    data_merged_taxa_info <- rbind(data_merged_taxa_info, data_taxa_info)
    
  }
  
}

# write.table(data_merged_taxa_info,
#             "/storm/User/ysm/Project/06.Microbiome/21.Green_Cleaner_paper/output/paper_output/processed_barplot/bar_plot/240429_merged_rev_taxa_info.txt",
#             quote = F, sep = "\t", row.names = F, col.names = T)