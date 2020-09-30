#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
fpath <- args[1]
files <- list.files(path= fpath, pattern="*.1.effect_sizes.tsv", full.names=TRUE, recursive=FALSE)

library(readr)
library(magrittr)
library(dplyr)

p <- read_tsv(files[1])
comp <- p[0,]
num_pcs = dim(p)[1]
num_samps = length(files)
no_sigs <-  c()
  
for (f in files)
{
  f_name <- basename(f)
  n <- sub('.1.effect_sizes.tsv', '', f_name) 
  t <- read_tsv(f) %>% mutate(Source = n) %>% mutate(bonf_p_true = pval * (num_pcs*num_samps))
  fin <- filter(t, bonf_p_true < 0.05)
    if(dim(fin)[1] == 0)
  {
    no_sigs <- c(no_sigs, n)
  }
  comp <- rbind(comp,t)
}

library(ggplot2)
print("Trait-specific analysis")
bonf_ <- comp %>% filter(bonf_p_true < 0.05)  %>% count(Source) %>% arrange(-n)
write_tsv(bonf_, "counts_per_trait.bonf_p.tsv") #phenotype in one column, number of PCs in the other
ggplot(data=bonf_, aes(n)) + geom_histogram() + xlab("number of PCs per trait") + ggtitle("Histogram of PCs per trait (bonf-p)")
ggsave("bonf_p_pcs_histogram.png")


pval_ <- comp %>% filter(pval < 0.05)  %>% count(Source) %>% arrange(-n)
write_tsv(pval_, "counts_per_trait.pval.tsv")
ggplot(data=pval_, aes(n)) + geom_histogram() + xlab("number of PCs per trait") + ggtitle("Histogram of PCs per trait (pval)")
ggsave("pval_pcs_histogram.png")


ggplot(data = bonf_, aes(x = reorder(Source, -n), y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90)) +
  xlab("Trait") + ylab("Count") + ggtitle("Significant PC associations per trait, bonfp")
ggsave("bonf_p_pcs_per_trait_all_traits.png")

ggplot(data = pval_, aes(x = reorder(Source, -n), y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90)) +
  xlab("Trait") + ylab("Count") + ggtitle("Significant PC associations per trait, pval")
ggsave("pval_pcs_per_trait_all_traits.png")


#Focus on the PCs
print("PC-specific analysis")
bonf_pc <- comp %>% filter(bonf_p_true < 0.05)  %>% count(annot) %>% arrange(n)
write_tsv(bonf_pc, "PCs_across_traits.bonf_pval.tsv") #PCs in one column, number of traits in next
ggplot(data=bonf_pc, aes(n)) + geom_histogram() + xlab("number of Traits per PC") + ggtitle("Histogram of Traits per PC (bonf-pval)")
ggsave("bonf-p_traits_histogram.png")


pval_pc <- comp %>% filter(pval < 0.05)  %>% count(annot) %>% arrange(n)
write_tsv(pval_pc, "PCs_across_traits.pval.tsv")
ggplot(data=pval_pc, aes(n)) + geom_histogram() + xlab("number of Traits per PC") + ggtitle("Histogram of Traits per PC (pval)")
ggsave("p_traits_histogram.png")

ggplot(data = bonf_pc, aes(x = reorder(annot, -n), y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle=90)) +
  xlab("PCs") + ylab("Count") + ggtitle(paste0("Frequency of bonf_p <= 0.05 associated annotation PCs with GWAS SumStats across ", length(files), " traits"))
ggsave("bonf_p_gwas_vs_pcs_all_traits.png")
ggplot(data = pval_pc, aes(x = reorder(annot, -n), y = n)) + geom_bar(stat = "identity")+ theme(axis.text.x = element_text(angle=90))+
  xlab("PCs") + ylab("Count") + ggtitle(paste0("Frequency of p <= 0.05 associated annotation PCs with GWAS Sum Stats across", length(files), " traits"))
ggsave("pval_gwas_vs_pcs_all_traits.png")


  
