#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3)
{
  print("Please provide arguments: input_file, output_file, sim_count")
  stop()
}
library(data.table)
library(tidyr)
library(stats)
library(dplyr)

#Load in the files

  if_name <- args[1]
  annot_set <- fread(if_name)
  c <- colnames(annot_set)
    annot_set <- data.frame(scale(annot_set))
    colnames(annot_set) <- c
  
  iter <- as.integer(args[3])
  pcs_list_ <- vector("list", iter)
  print(iter)
    for (i in 1:iter)
  {
    if(i %% (iter*0.1) == 0) { print(i) }
     shuffled <- as.data.frame(lapply(annot_set, sample))
shuffled_pcs <- prcomp(shuffled, scale. = F, center = F, rank. = 100)
shuffled_cor <- data.frame(cor(data.frame(shuffled_pcs$x), shuffled))
    pcs_list_[[i]] <- data.frame(shuffled_cor)
  }
  save(pcs_list_, file =args[2])


