# Basic
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

# Plotting
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggpubr))

rm(list = ls())

# 1. variables and function ----------------------------------------------------

path_data <- commandArgs(T)[1]

path_alpha <-
  glue("{path_data}/../output/qiime/alpha_diversity")
path_output <-
  glue("{path_data}/../output/paper_output/alpha_diversity")

input_meta <-
  fread(glue("{path_data}/input_sample_meta_data.txt"),
        sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F)

if (!dir.exists(path_output)) {
  
  dir.create(path = path_output, recursive = T)
  
}

make.df <- function(x, y, z) {
  
  sapply(strsplit(x = x, split = y), '[[', z)
  
}

input_dir_list <-
  c("raw_chao1_ci", "scrub_chao1_ci", glue("dataset_{seq(10)}_gc_chao1_ci"))

# 2. pre processing ------------------------------------------------------------

for (value_tmp_dir in input_dir_list) {

  # value_tmp_dir <- input_dir_list[1]
  
  data_tmp_chao1 <- 
    fread(glue("{path_alpha}/{value_tmp_dir}/alpha-diversity.tsv"), 
          sep = "\t", quote = F, header = T, stringsAsFactors = F, data.table = F) %>%
    rename("sample_id" = "V1") %>%
    rename("value" = 2) %>%
    mutate(value = make.df(value, ",", 1)) %>%
    mutate(value = make.df(value, "\\(", 2)) %>%
    mutate(value = as.numeric(value))
  
  if (grepl(pattern = "raw", x = value_tmp_dir)) {
    
    data_tmp_chao1 <-
      data_tmp_chao1 %>%
      mutate(type = "Undiluted")
    
  } else if (grepl(pattern = "scrub", x = value_tmp_dir)) {
    
    data_tmp_chao1 <-
      data_tmp_chao1 %>%
      mutate(type = "SCRuB")
    
  } else if (grepl(pattern = "gc", x = value_tmp_dir)) {
    
    data_tmp_chao1 <-
      data_tmp_chao1 %>%
      mutate(type = "Green Cleaner")
    
  }
  
  if (value_tmp_dir == input_dir_list[1]) {
    
    data_merged_chao1 <- data_tmp_chao1
    
  } else {
    
    data_merged_chao1 <- rbind(data_merged_chao1, data_tmp_chao1)
    
  }
  
}

data_merged_table <-
  merge(data_merged_chao1,
        input_meta %>%
          select(sample_id, dilution_series),
        by = "sample_id") %>%
  filter(dilution_series != "-") %>%
  filter(dilution_series != "d0") %>%
  mutate(dilution_series = factor(toupper(dilution_series), levels = glue("D{seq(7)}"))) %>%
  mutate(Type = factor(type, levels = c("Undiluted", "SCRuB", "Green Cleaner"))) %>%
  rename("Dilution Series" = "dilution_series") %>%
  mutate(value = log10(value))

# 3. box plot ------------------------------------------------------------

my_comparisons <- list(c("SCRuB", "Green Cleaner"))

fig_chao <- ggboxplot(data_merged_table, 
          x = "Type", y = "value",
          color = "Type", palette = "npg",
          add = "jitter",
          facet.by = "Dilution Series", short.panel.labs = T) +
  stat_compare_means(comparisons = my_comparisons, 
                     label.y = 2.35, method = "wilcox", paired = TRUE, label = "p.signif") + 
  theme_bw() + 
  ylim(0, 2.6) +
  scale_color_manual(values = c("Undiluted" = "red3", "SCRuB" = "blue3", "Green Cleaner" = "green3")) + 
  ylab("Chao1 value")

png(glue("{path_output}/fig4_alpha_diversity_figure.png"), width = 600, height = 500)

print(fig_chao)

dev.off()

write.table(data_merged_table, glue("{path_output}/fig4_alpha_diversity_table.txt"), 
            quote = F, sep = "\t", row.names = F, col.names = T)
