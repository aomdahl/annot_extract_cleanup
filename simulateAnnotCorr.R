#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 5)
{
  print("Please provide arguments: input_file, output_file, sim_count, num_svs number_cores")
  stop()
}
library(data.table)
library(tidyr)
library(stats)
library(dplyr)
library(RSpectra)
library(parallel)
library(MASS)


shuffleIteration <- function(iter, annot_set, sv_count)
{
    if(iter %% 100 == 0) { print(iter) }
    shuffled <- as.data.frame(lapply(annot_set, sample))
    sing_val_d <- svds(as.matrix(shuffled), sv_count)
    shuffled_cor <- data.frame(cor(data.frame(sing_val_d$u), shuffled))
    data.frame(shuffled_cor)
}
#Load in the files
print(args)
print("Arguments as above...")
  if_name <- args[1]
  annot_set <- fread(if_name)
  c <- colnames(annot_set)
    annot_set <- data.frame(scale(annot_set))
    colnames(annot_set) <- c

  iter <- as.integer(args[3])
  pcs_list_ <- vector("list", iter)
  print(iter)
  reps <- 1:iter
#results <- lapply(reps, function (x)  shuffleIteration(x, annot_set, 150))
#Paralleleized version
numCores <- detectCores()
if(numCores > args[5])
{
    numCores = args[5]
}

pcs_list_ <- mclapply(reps, function (x)  shuffleIteration(x, annot_set, 150), mc.cores = numCores)

save(pcs_list_, file =args[2])

