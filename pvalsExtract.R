#!/usr/bin/env Rscript

library(data.table)
library(tidyr)
library(readr)
library(stats)
library(stringr)
library(Xmisc)
library(RSpectra)
library(dplyr)
parser <- ArgumentParser$new()
parser$add_description("Script for calculating p-values based on simulated annotation data. Also calculates PCs for annotation data.")
parser$add_argument("--input_annots", type = 'character', help = "File containing annotations in a tsv.")
parser$add_argument("--input_ids", type = 'character', help = "File containing the SNP IDs for the input annots file.")
parser$add_argument("--output", type = 'character', help = "Path and name for output file.")
parser$add_argument("--sim_data", type = "character", help = "Path to simluation results from other script, .RData file.")
parser$add_argument("--pca_only", type = "logical", help = "specify to just return PC data from annotation inputs.", default =F)
parser$add_argument("--num_svs", type = "numeric", help = "Specify the number of SVs to include in each run. Default is 100", default = 100)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()

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


print("Reading in annotation data.")
annots_data <- fread(args$input_annots) %>% scale()
print("Calculating annotation PCs now...")
#Do this with Rspectra.

sing_val_d <- svds(as.matrix(annots_data), args$num_svs)
pcs <- data.frame(sing_val_d$u)
col_names <- paste0("PC", 1:args$num_svs)
  names(pcs) <- col_names
#pcs_v <- data.frame(sing_val_d$v)

if (args$pca_only)
{
     ids <- fread(args$input_ids,header = F) %>% rename("SNP" = V1) %>% separate(SNP, into = c("Chr", "Pos", "Ref", "Alt"), sep = ":", remove = F) %>% select(-Ref, -Alt)
    #Assuming the order is the same
    pcs_ <- data.frame(cbind(ids, round(pcs, digits = 6)))
    write_tsv(pcs_, paste0(args$output, "_pcs.tsv"))
    pve <- (sing_val_d$d)^2/(sum((sing_val_d$d)^2))
    png(paste0(args$output, "_annots_pve.png"))
    plot(pve, col = "deepskyblue", ylab = "Percent Variance Explained", xlab = "PC", pch = 19)
    dev.off()
    print("PCs written to output and PVE plot saved. Program will terminate.")
    quit()
}

print("Loading in simulation data...")
load(args$sim_data)

print("Calculating PC correlation with annotations...")
true_pcs <- data.frame(cor(data.frame(pca_$x), annots_data))
print("Correlation calculation complete")
print(paste("Correlation matrix has dimensions",dim(true_pcs)[1], dim(true_pcs)[2]))

npcs = args$num_svs  #Default value we have been using.
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
write_tsv(pvals_, args$output)


