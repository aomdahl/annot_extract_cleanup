#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
fpath <- args[1]
files <- list.files(path= fpath, pattern="*.1.effect_sizes.tsv", full.names=TRUE, recursive=FALSE)
bonf_thresh <- args[2]
library(readr)
library(magrittr)
library(dplyr)

col_types <- cols(
  annot = col_character(),
  R2 = col_double(),
  beta = col_double(),
  t_stat = col_double(),
  pval = col_double(),
  std_error = col_double(),
  abs_beta = col_double(),
  abs_t = col_double(),
  bonf_p = col_double()
)



p <- read_tsv(files[1])
comp <- p[0,]
num_pcs = dim(p)[1]
num_samps = length(files)
no_sigs <-  c()
  
for (f in files)
{
  f_name <- basename(f)
  n <- sub('.1.effect_sizes.tsv', '', f_name) 
  t <- read_tsv(f, col_types = col_types) %>% mutate(Source = n) %>% mutate(bonf_p_true = pval * (num_pcs*num_samps))
  fin <- filter(t, bonf_p_true < bonf_thresh)
    if(dim(fin)[1] == 0)
  {
    no_sigs <- c(no_sigs, n)
  }
  comp <- rbind(comp,t)
}
print(comp)
print("We have compiled all the files now...")
print(typeof(bonf_thresh))
library(ggplot2)
corrected_threshold <- as.numeric(bonf_thresh) / (num_pcs * num_samps)
print("Bonf-adjusted threshold:")
print(corrected_threshold)

#Look at the traits
print("Trait-specific analysis")
bonf_ <- comp %>% filter(bonf_p_true < bonf_thresh)  %>% count(Source) %>% arrange(-n) %>% mutate(has = "has")
missingno <- comp %>% filter(!(Source %in% bonf_$Source)) %>% count(Source) %>% mutate(n = 0, has = "hasnot") #Select all those that had NO below our threshold, set them to 0
if (any(bonf_$Source %in% missingno$Source))
{
    print(bonf_[which(bonf_$Source %in% missingno$Source),])
}
bonf_ <- rbind(bonf_, missingno)
write_tsv(bonf_, "counts_per_trait.bonf_p.tsv") #phenotype in one column, number of PCs in the other
ggplot(data=bonf_, aes(n)) + geom_histogram(fill="skyblue") + xlab("number of PCs per trait") + ggtitle("Histogram of PCs per trait (bonf-p)")
ggsave("bonf_p_pcs_histogram.png")


pval_ <- comp %>% filter(pval < bonf_thresh)  %>% count(Source) %>% arrange(-n)
missingno <- comp %>% filter(!(Source %in% pval_$Source)) %>% count(Source) %>% mutate(n = 0)
pval_ <- rbind(pval_, missingno)
write_tsv(pval_, "counts_per_trait.pval.tsv")
ggplot(data=pval_, aes(n)) + geom_histogram(fill="skyblue") + xlab("number of PCs per trait") + ggtitle("Histogram of PCs per trait (pval)")
ggsave("pval_pcs_histogram.png")


ggplot(data = bonf_, aes(x = reorder(Source, -n), y = n)) + geom_bar(stat = "identity", fill = "skyblue") + theme(axis.text.x = element_text(angle=90)) +
  xlab("Trait") + ylab("Count") + ggtitle(paste("Significant PC associations per trait (p < ",round(corrected_threshold, digits = 8), ")"))
ggsave("bonf_p_pcs_per_trait_all_traits.png")

ggplot(data = pval_, aes(x = reorder(Source, -n), y = n)) + geom_bar(stat = "identity", fill = "skyblue") + theme(axis.text.x = element_text(angle=90)) +
  xlab("Trait") + ylab("Count") + ggtitle("Significant PC associations per trait, pval")
ggsave("pval_pcs_per_trait_all_traits.png")


#Focus on the PCs
print("PC-specific analysis")
bonf_pc <- comp %>% filter(bonf_p_true < bonf_thresh)  %>% count(annot) %>% arrange(n)
missingno <- comp %>% filter(!(annot %in% bonf_pc$annot)) %>% count(annot) %>% mutate(n = 0)
bonf_pc <- rbind(bonf_pc, missingno)
write_tsv(bonf_pc, "PCs_across_traits.bonf_pval.tsv") #PCs in one column, number of traits in next
ggplot(data=bonf_pc, aes(n)) + geom_histogram( fill = "skyblue") + xlab("number of Traits per PC") + ggtitle("Histogram of Traits per PC (bonf-pval)")
ggsave("bonf-p_traits_histogram.png")


pval_pc <- comp %>% filter(pval < bonf_thresh)  %>% count(annot) %>% arrange(n)
missingno <- comp %>% filter(!(annot %in% pval_pc$annot)) %>% count(annot) %>% mutate(n = 0)
pval_pc <- rbind(pval_pc, missingno)
write_tsv(pval_pc, "PCs_across_traits.pval.tsv")

ggplot(data=pval_pc, aes(n)) + geom_histogram( fill = "skyblue") + xlab("number of Traits per PC") + ggtitle("Histogram of Traits per PC (pval)")
ggsave("p_traits_histogram.png")

ggplot(data = bonf_pc, aes(x = reorder(annot, -n), y = n)) + geom_bar(stat = "identity", fill = "skyblue") + theme(axis.text.x = element_text(angle=90)) +
  xlab("PCs") + ylab("Count") + ggtitle(paste0("Frequency of p < ", round(corrected_threshold, digits = 8),  "associated annotation PCs with GWAS SumStats across ", length(files), " traits"))
ggsave("bonf_p_gwas_vs_pcs_all_traits.png")
ggplot(data = pval_pc, aes(x = reorder(annot, -n), y = n)) + geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle=90))+
  xlab("PCs") + ylab("Count") + ggtitle(paste0("Frequency of p <= 0.05 associated annotation PCs with GWAS Sum Stats across", length(files), " traits"))
ggsave("pval_gwas_vs_pcs_all_traits.png")


  
