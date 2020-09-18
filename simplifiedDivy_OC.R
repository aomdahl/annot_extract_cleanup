library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)




annot_names <- scan("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/header.txt", what = character())
max_annots <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/all_annotations.250-kb-win_0.5_r2.tsv")
print("Now dividing by MAF and coding/non-coding. For this simple analysis, we only look at variants with MAF > 0.05")
names(max_annots) <- annot_names
tog <- merge(max_annots, comb, by="F_ID")

#Clean up the annotation data a little bit more.
binarized_maf_common <- (numeric_annots %>% select(MAFbin_frequent_1,MAFbin_frequent_2,MAFbin_frequent_3,MAFbin_frequent_4,MAFbin_frequent_5,MAFbin_frequent_6,MAFbin_frequent_7,MAFbin_frequent_8,MAFbin_frequent_9,MAFbin_frequent_10)) %>% rowSums()

common_freq <- numeric_annots[as.logical(binarized_maf_common),] %>% select(-MAFbin_lowfreq_1,-MAFbin_lowfreq_2,-MAFbin_lowfreq_3,-MAFbin_lowfreq_4,-MAFbin_lowfreq_5,-MAFbin_lowfreq_6,-MAFbin_lowfreq_7,-MAFbin_lowfreq_8,-MAFbin_lowfreq_9,-MAFbin_lowfreq_10,-MAFbin_frequent_1,-MAFbin_frequent_2,-MAFbin_frequent_3,-MAFbin_frequent_4,-MAFbin_frequent_5,-MAFbin_frequent_6,-MAFbin_frequent_7,-MAFbin_frequent_8,-MAFbin_frequent_9,-MAFbin_frequent_10) 
#Remove ambifuous variants
common_freq <- common_freq %>% mutate(ambig_check = paste0(Ref,Alt)) %>% filter(!(ambig_check %in% ambig)) %>% select(-ambig_check)

n_coding <- common_freq %>% filter(Coding_UCSC_common == 1)
common_noncoding <- common_freq %>% filter(Coding_UCSC_common == 0)
#filter out weird ones from CADD
annot_names <- names(common_freq)
indx <- grepl('.na', annot_names)
common_freq <- common_freq[,!indx]
#Write out to file.
write_tsv(common_coding, "common_coding_oc_join.sept2020.simple.tsv")
write_tsv(common_noncoding, "common_noncoding_oc_join.sept2020.simple.tsv")
