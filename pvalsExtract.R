#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3)
{
  print("Please provide arguments: pcs_rdata_file, true_annotations_data, output_file")
  stop()
}
library(data.table)
library(tidyr)
library(readr)
library(stats)
library(stringr)
getPC <- function(df, pc_num)
{
  return(df[pc_num,])
}

loadSims <- function(sims_list)
{
   simdat <- str_split_fixed(sims_list, ",", n = Inf) 
    base <- load(simdat[1])
    for (i in 2:length(simdat))
    {
        add <- load(simdat[i])
        base <- rbind(base, add)
    }
    return(base)
}


print("Loading in simulation data...")
load(args[1]) #pcs_list_
print("Reading in annotation data.")
annots_data <- fread(args[2])
print("Calculating annotation PCs now...")
pca_ <- prcomp(annots_data, scale. = T, center = T, rank. = 100)
print("Calculating PC correlation with annotations...")
true_pcs <- data.frame(cor(data.frame(pca_$x), annots_data))
print("Correlation calculation complete")
print(paste("Correlation matrix has dimensions",dim(true_pcs)[1], dim(true_pcs)[2]))

npcs = 100 #Default value we have been using.
#get the first one out to get the dimensions
all_of_one_pc <- as.matrix(t(mapply(getPC, pcs_list_,  pc_num = 1)))
nfeats <- dim(all_of_one_pc)[2]
nsim <- dim(all_of_one_pc)[1]
pvals <- matrix(1,npcs,nfeats)
print("Beginning extraction of simulation results...")
for(i in 1:npcs)
{
  all_of_one_pc <- as.matrix(t(mapply(getPC, pcs_list_,  pc_num = i)))
  true_pc <- true_pcs[i,] #same number of annotations
  sims_squared <- data.frame(apply(all_of_one_pc, 2, function(x) unlist(x)^2))
  true_squared <- true_pc*true_pc
  true_squared <- rbind(true_squared, true_squared[rep(1,nsim-1),])
  significance_sq <- colSums(sims_squared >= true_squared)/nsim #how often is a permuted one greater than the true one?
  pvals[i,] = significance_sq
  if(i %% (npcs*0.1) == 0) { print(i)}
}
PC <- paste0("PC", 1:npcs)
pvals_ <- data.frame(pvals)
names(pvals_) <- names(true_pc)
pvals_ <- cbind(PC, pvals_)
write_tsv(pvals_, args[3])


