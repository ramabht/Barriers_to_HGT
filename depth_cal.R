#! /home/rama/R-4.0.0/bin/Rscript
library(data.table)
library(fs)
library(dplyr)

setwd("/home/rama/lmamr_uol_NGS/DeepSeq/RawReads_copy_02/Gene_depths/depth_files_for_Rscript")

file_path = fs::dir_ls("/home/rama/lmamr_uol_NGS/DeepSeq/RawReads_copy_02/Gene_depths/depth_files_for_Rscript")
file_path

file_contents <- list()
for(i in seq_along(file_path)){
  file_contents[[i]] <- read.table(file = file_path[[i]], header = T)
}

file_contents <- setNames(file_contents, file_path)

names(file_contents) <- basename(file_path)

names(file_contents) <- tools::file_path_sans_ext(names(file_contents))

print(names(file_contents))

for (i in 1:length(file_contents)){
  name = names(file_contents[i])
  medianfile = file_contents[[i]] %>%
    group_by(gene) %>%
    summarise(median = median(depth, na.rm = TRUE))
  filename <- paste("genedepths_",name,".csv")
  write.table(medianfile, filename, quote = F, sep = ",", row.names = F)
    }


