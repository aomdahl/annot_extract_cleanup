#This script joins the annotation data sources, and extracts only common variants.
#As of Sept 30, this does not retain maf info.
library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)

annot_names <- scan("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/header.txt", what = character())
max_annots <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/all_annotations.250-kb-win_0.5_r2.tsv")
comb <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/openCravat/openCravat_250kb-0.5r2_clean.tsv")
names(max_annots) <- annot_names

tog <- merge(max_annots, comb, by="F_ID")

print("Total combined LDSC, CADD and OpenCravat data...")
print(dim(tog)[1])
print(dim(tog)[2])




#Merge in the MAF data:

print("Accounting for LDSC MAF bins....")

#Clean up the annotation data a little bit more.
binarized_maf_common <- (tog %>% select(MAFbin_frequent_1,MAFbin_frequent_2,MAFbin_frequent_3,MAFbin_frequent_4,MAFbin_frequent_5,MAFbin_frequent_6,MAFbin_frequent_7,MAFbin_frequent_8,MAFbin_frequent_9,MAFbin_frequent_10))



#Put them into those windows

common_freq <- tog[as.logical(rowSums(binarized_maf_common)),]



#Get maf info
maf_tab <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/annot_extract_cleanup/ldsc_maf_bins.2018.txt") %>% mutate("avg" = (lower+upper)/2)

locale <- apply(binarized_maf_common, 1, function(x) (which(x == 1)))
locale <- unlist(locale)
maf_tab_lookup <- maf_tab[6:nrow(maf_tab),]
maf_bin_data <- maf_tab_lookup[locale, 4]
common_freq$`ldsc_maf` <- maf_bin_data

common_freq <- common_freq %>% select(-MAFbin_lowfreq_1,-MAFbin_lowfreq_2,-MAFbin_lowfreq_3,-MAFbin_lowfreq_4,-MAFbin_lowfreq_5,-MAFbin_lowfreq_6,-MAFbin_lowfreq_7,-MAFbin_lowfreq_8,-MAFbin_lowfreq_9,-MAFbin_lowfreq_10,-MAFbin_frequent_1,-MAFbin_frequent_2,-MAFbin_frequent_3,-MAFbin_frequent_4,-MAFbin_frequent_5,-MAFbin_frequent_6,-MAFbin_frequent_7,-MAFbin_frequent_8,-MAFbin_frequent_9,-MAFbin_frequent_10)

print(common_freq)



#filter out weird ones from CADD
annot_names <- names(common_freq)
indx <- grepl('.na', annot_names)
#common_freq <- common_freq[,!indx]
common_freq <- common_freq[,!indx, with = F]



print("Limited to common variants only:")
print(dim(common_freq)[1])
print(dim(common_freq)[2])

#Remove ambiguous variant
print("Removing ambiguous variants....")
ambig = c("AT","TA", "GC", "CG")
common_freq <- common_freq %>% mutate(ambig_check = paste0(Ref, Alt)) %>% filter(!(ambig_check %in% ambig) | ldsc_maf >= 0.3) %>% select(-ambig_check)

print("AFter removing ambiguous variants:")
print(dim(common_freq)[1])
print(dim(common_freq)[2])


common_coding <- common_freq %>% filter(Coding_UCSC_common == 1) %>% select(-ldsc_maf)
common_noncoding <- common_freq %>% filter(Coding_UCSC_common == 0) %>% select(-ldsc_maf)
#Write out to file.
write_tsv(common_coding, "common_coding_no-maf.sept2020.simple.tsv")

write_tsv(common_noncoding, "common_noncoding_no-maf.sept2020.simple.tsv")

