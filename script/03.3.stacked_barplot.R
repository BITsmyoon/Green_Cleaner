# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

# Plotting
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Polychrome))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)
input_taxa_color_list <-
  fread(glue("{path_data}/input_stacked_barplot_color_list.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

path_qiime <-
  glue("{path_data}/../output/qiime")
path_output <-
  glue("{path_data}/../output/paper_output/stacked_barplot")

input_dataset_list <-
  glue("dataset_{seq(10)}")
input_all_level_7_count <-
  read.csv(glue("{path_qiime}/qiime_barplot/level-7.csv"), sep = ",", header = T, stringsAsFactors = F, check.names = F) %>%
  select(-dataset)

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

make.df <- function(x, y, z) {
  
  sapply(strsplit(x = x, split = y), '[[', z)
  
}

# 2. pre processing and plot ---------------------------------------------------

value_color <- input_taxa_color_list$color

value_num <- 0
for (value_dataset in input_dataset_list) {
  
  data_dataset_meta <-
    input_meta %>%
    filter(dataset == value_dataset) %>%
    filter(dilution_series != "-")
  
  data_level_7_count <-
    input_all_level_7_count %>%
    filter(index %in% data_dataset_meta$sample_id) %>%
    column_to_rownames("index")
  
  data_raw_relab <- 
    apply(data_level_7_count, 1, function(value){
    (value / sum(value)) * 100
  }) %>%
    as.data.frame()

  colnames(data_raw_relab) <- 
    glue("D{0:7}")
  
  rownames(data_raw_relab) <- gsub(pattern = ";s_$", replacement = "", rownames(data_raw_relab))
  rownames(data_raw_relab) <- gsub(pattern = ";g_$", replacement = "", rownames(data_raw_relab))
  rownames(data_raw_relab) <- gsub(pattern = ";f_$", replacement = "", rownames(data_raw_relab))
  rownames(data_raw_relab) <- gsub(pattern = ";o_$", replacement = "", rownames(data_raw_relab))
  rownames(data_raw_relab) <- gsub(pattern = ";c_$", replacement = "", rownames(data_raw_relab))
  rownames(data_raw_relab) <- gsub(pattern = ";p_$", replacement = "", rownames(data_raw_relab))
    
  value_taxa_len <- 
    strsplit(rownames(data_raw_relab), ";") %>%
    sapply(., length)
  
  for (value_row_num in seq(nrow(data_raw_relab))) {
    
    value_target_taxa <- 
      rownames(data_raw_relab)[value_row_num]
    value_target_len <-
      value_taxa_len[value_row_num]
    
    if (value_target_len == 1) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "k_", replacement = "K_", .)
      
    } else if (value_target_len == 2) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "p_", replacement = "P_", .)
      
    } else if (value_target_len == 3) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "c_", replacement = "C_", .)
      
    } else if (value_target_len == 4) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "o_", replacement = "O_", .)
      
    } else if (value_target_len == 5) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "f_", replacement = "F_", .)
      
    } else if (value_target_len == 6) {
      
      rownames(data_raw_relab)[value_row_num] <-
        make.df(value_target_taxa, ";", value_target_len) %>%
        gsub(pattern = "g_", replacement = "G_", .)
      
    } else if (value_target_len == 7) {
      
      rownames(data_raw_relab)[value_row_num] <-
        glue('{make.df(value_target_taxa, ";", value_target_len - 1)}_{make.df(value_target_taxa, ";", value_target_len) %>%
        make.df(., "s_", 2)}') %>%
        gsub(pattern = "g_", replacement = "S_", .)
      
    }
    
  }
  
  data_arranged_taxa_relab <- 
    data_raw_relab %>%
    arrange(desc(D0))
  
  data_plot_input <- 
    data_arranged_taxa_relab %>%
    rownames_to_column("taxa") %>%
    pivot_longer(cols = starts_with("D"), 
                 names_to = "d_series", 
                 values_to = "value") %>%
    arrange(d_series) %>%
    mutate(dataset = value_dataset) %>%
    mutate(taxa = factor(taxa, levels = rownames(data_arranged_taxa_relab))) %>%
    mutate(d_series = factor(d_series, levels = colnames(data_arranged_taxa_relab)))
  
  data_taxa_info <- 
    data.frame(taxa = rownames(data_arranged_taxa_relab),
               d0_in = ifelse(data_arranged_taxa_relab$D0 != 0, T, F),
               d0_5_up = ifelse(data_arranged_taxa_relab$D0 >= 5, T, F)) %>%
    mutate(the_others_check = ifelse((d0_in == T) & (d0_5_up == F), T, F))

  data_taxa_info %>%
    filter(d0_5_up) %>%
    mutate(index = seq(nrow(.)))
  
  value_viridis_palette <- 
    merge(data_taxa_info %>%
            filter(d0_5_up) %>%
            mutate(index = seq(nrow(.))),
          input_taxa_color_list,
          by = "taxa") %>%
    arrange(index) %>%
    pull("color")
  value_gray_palette <- c("#6C6C6C", "#8D8D8D", "#A7A7A7", "#BDBDBD")
  value_the_others_num <- sum(data_taxa_info$the_others_check)
  
  value_nrow <- nrow(data_taxa_info)
  value_color_count <- sum(data_taxa_info$d0_5_up) + value_the_others_num
  value_rep_gray <- (value_nrow - value_color_count) %/% 4
  value_other <- value_nrow - (value_color_count + (value_rep_gray * 4))
  
  if (value_other != 0) {
    
    data_taxa_all_info <- 
      data_taxa_info %>%
      mutate(color = c(value_viridis_palette,
                       rep("#3B00FB", value_the_others_num),
                       rep(value_gray_palette, value_rep_gray),
                       value_gray_palette[seq(value_other)]))
    
  } else {
    
    data_taxa_all_info <- 
      data_taxa_info %>%
      mutate(color = c(value_viridis_palette,
                       rep("#3B00FB", value_the_others_num),
                       rep(value_gray_palette, value_rep_gray)))
    
  }
  
  value_true_signal_num <- 
    data_taxa_all_info %>%
    filter(d0_5_up) %>%
    nrow()
  
  value_the_other_num <- 
    sum(data_taxa_all_info$the_others_check)
  
  value_legend_vector <- 
    setNames(data_taxa_all_info$color, data_taxa_all_info$taxa)
  
  value_true_signal_taxa_name <- 
    names(value_legend_vector)[seq(value_true_signal_num)]
  value_contaminant_taxa_name <- 
    names(value_legend_vector)[seq(from = value_true_signal_num + value_the_other_num + 1, to = value_true_signal_num + value_the_other_num + 4)]
  
  if (value_true_signal_num != 0) {
    
    value_the_others_taxa_name <- names(value_legend_vector)[value_true_signal_num + 1]
    value_full_in_taxa_name <- c(value_true_signal_taxa_name, value_the_others_taxa_name, value_contaminant_taxa_name)
    
  } else {
    
    value_full_in_taxa_name <- c(value_true_signal_taxa_name, value_contaminant_taxa_name)
    
  }
  
  value_all_target_taxa <- 
    data_taxa_all_info %>%
    filter(d0_5_up == T) %>%
    pull("taxa")
  
  for (value_target_taxa in value_all_target_taxa) {
    
    value_target_taxa <- value_all_target_taxa[1]
    
    value_target_color <- 
      input_taxa_color_list %>%
      filter(taxa == value_target_taxa) %>%
      pull("color")
    
    data_taxa_all_info[data_taxa_all_info$taxa == value_target_taxa, ]$color <- value_target_color
    
  }
  
  data_figure <-
    ggplot(data_plot_input, aes(fill = taxa, y = value, x = d_series)) + 
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual("Bacteria",
                      values = value_legend_vector, 
                      labels = c(value_true_signal_taxa_name, "The others", rep("Contaminant", 4)),
                      breaks = value_full_in_taxa_name) +
    theme_bw() + 
    labs(y = "Relative abundance", x = "Dilution series and Type") +
    scale_x_discrete(guide = guide_axis(n.dodge = 2))
  
  png(glue("{path_output}/{value_dataset}_stacked_barplot.png"), width = 1200, height = 600)
  
  print(data_figure)
  
  dev.off()

  write.table(data_plot_input, glue("{path_output}/{value_dataset}_stacked_barplot_table.txt"),
              quote = F, sep = "\t", row.names = F, col.names = T)  
  
}