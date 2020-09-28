#This script joins the annotation data sources, and extracts only common variants.
#It also appends MAF data from a few different sources
#As of Sept 25, it keeps in MAF bins.
library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)

annot_names <- scan("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/header.txt", what = character())
max_annots <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/all_annotations.250-kb-win_0.5_r2.tsv")
comb <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/openCravat/openCravat_250kb-0.5r2_clean.tsv")
names(max_annots) <- annot_names


thousandG_maf <- fread("/work-zfs/abattle4/ashton/prs_dev/1000genomes_refLD/plink.frq") %>% rename("S_ID" = "SNP") %>% mutate("FF_ID" = paste0(S_ID, ":", A2, ":", A1), "FB_ID" = paste0(S_ID, ":", A1, ":", A2))
print("Now dividing by MAF and coding/non-coding. For this simple analysis, we only look at variants with MAF > 0.05")

names(max_annots) <- annot_names
tog <- merge(max_annots, comb, by="F_ID")

print("Total combined LDSC, CADD and OpenCravat data...")
print(dim(tog)[1])
print(dim(tog)[2])




#Merge in the MAF data:
print("Adding in MAF data from UKBB and 1000G")
FF_ID <- tog %>% filter(F_ID %in% thousandG_maf$FF_ID) #try different possible IDs, 
FB_ID <- tog %>% filter(F_ID %in% thousandG_maf$FB_ID)
tog_maf <- rbind(FB_ID, FF_ID) 
keep_tg1 <-thousandG_maf[which(thousandG_maf$FF_ID %in% tog_maf$F_ID),] %>% select(FF_ID, MAF) %>% rename("F_ID" = FF_ID)
keep_tg2 <- thousandG_maf[which(thousandG_maf$FB_ID %in% tog_maf$F_ID),] %>% select(FB_ID, MAF) %>% rename("F_ID" = FB_ID)
thousandG_merge <- rbind(keep_tg1, keep_tg2)
tog_maf <- merge(tog_maf, thousandG_merge) #adding it in actually
rm(thousandG_merge, FF_ID,FB_ID,keep_tg1,keep_tg2)

ukbb_maf <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/ukbb_bmi-subset_maf.txt")
ukbb_maf <- ukbb_maf %>% select(-low_confidence_variant, -minor_allele) %>% rename("F_ID" = variant)
tog_maf <- merge(tog_maf, ukbb_maf, by = "F_ID")
rm(ukbb_maf)

print("Merged with UKBB and  1000G data:")
print(dim(tog_maf)[1])
print(dim(tog_maf)[2])
print("Note that loss in variants indicates that the Minor/major alleles in 1000G didn't quite line up with the Ref/Alt given in our UKBB reference")

print("Accounting for LDSC MAF bins....")

#Clean up the annotation data a little bit more.
binarized_maf_common <- (tog_maf %>% select(MAFbin_frequent_1,MAFbin_frequent_2,MAFbin_frequent_3,MAFbin_frequent_4,MAFbin_frequent_5,MAFbin_frequent_6,MAFbin_frequent_7,MAFbin_frequent_8,MAFbin_frequent_9,MAFbin_frequent_10)) 



#Put them into those windows

common_freq <- tog_maf[as.logical(rowSums(binarized_maf_common)),] 

#I don't think we want to remove these anymore... useful for assesing things, which variants are most likely to be interested, etc.
#%>% select(-MAFbin_lowfreq_1,-MAFbin_lowfreq_2,-MAFbin_lowfreq_3,-MAFbin_lowfreq_4,-MAFbin_lowfreq_5,-MAFbin_lowfreq_6,-MAFbin_lowfreq_7,-MAFbin_lowfreq_8,-MAFbin_lowfreq_9,-MAFbin_lowfreq_10,-MAFbin_frequent_1,-MAFbin_frequent_2,-MAFbin_frequent_3,-MAFbin_frequent_4,-MAFbin_frequent_5,-MAFbin_frequent_6,-MAFbin_frequent_7,-MAFbin_frequent_8,-MAFbin_frequent_9,-MAFbin_frequent_10) 

maf_tab <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/annot_extract_cleanup/ldsc_maf_bins.2018.txt") %>% mutate("avg" = (lower+upper)/2)

locale <- apply(binarized_maf_common, 1, function(x) (which(x == 1)))
locale <- unlist(locale)
maf_tab_lookup <- maf_tab[6:nrow(maf_tab),]

#maf_bin_data <- apply(locale, 1, function(x) maf_tab_lookup[x,4])
maf_bin_data <- maf_tab_lookup[locale, 4]
common_freq$`ldsc_maf` <- maf_bin_data


#filter out weird ones from CADD
annot_names <- names(common_freq)
indx <- grepl('.na', annot_names)
common_freq <- common_freq[,!indx]
#common_freq <- common_freq[,!indx, with = F]

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


common_coding <- common_freq %>% filter(Coding_UCSC_common == 1)
common_noncoding <- common_freq %>% filter(Coding_UCSC_common == 0)
#Write out to file.
write_tsv(common_coding, "common_coding_oc_join.sept2020.simple.tsv")

write_tsv(common_noncoding, "common_noncoding_oc_join.sept2020.simple.tsv")
