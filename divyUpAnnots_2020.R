library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)



#Remove columns that have 0 variance
varClear <- function(tab)
{
  zv <- names(which(apply(tab, 2, var) == 0))
  tab %>% select(-zv)
}


numeric_annots <- fread("cleaned_up_numeric_annotations.tsv")
print("Now dividing by MAF and coding/non-coding")
#Clean up the annotation data a little bit more.
binarized_maf_lowfreq <- (numeric_annots %>% select(MAFbin_lowfreq_1,MAFbin_lowfreq_2,MAFbin_lowfreq_3,MAFbin_lowfreq_4,MAFbin_lowfreq_5,MAFbin_lowfreq_6,MAFbin_lowfreq_7,MAFbin_lowfreq_8,MAFbin_lowfreq_9,MAFbin_lowfreq_10)) %>% rowSums()
binarized_maf_common <- (numeric_annots %>% select(MAFbin_frequent_1,MAFbin_frequent_2,MAFbin_frequent_3,MAFbin_frequent_4,MAFbin_frequent_5,MAFbin_frequent_6,MAFbin_frequent_7,MAFbin_frequent_8,MAFbin_frequent_9,MAFbin_frequent_10)) %>% rowSums()

low_freq <- numeric_annots[as.logical(binarized_maf_lowfreq),] %>% select(-MAFbin_lowfreq_1,-MAFbin_lowfreq_2,-MAFbin_lowfreq_3,-MAFbin_lowfreq_4,-MAFbin_lowfreq_5,-MAFbin_lowfreq_6,-MAFbin_lowfreq_7,-MAFbin_lowfreq_8,-MAFbin_lowfreq_9,-MAFbin_lowfreq_10,-MAFbin_frequent_1,-MAFbin_frequent_2,-MAFbin_frequent_3,-MAFbin_frequent_4,-MAFbin_frequent_5,-MAFbin_frequent_6,-MAFbin_frequent_7,-MAFbin_frequent_8,-MAFbin_frequent_9,-MAFbin_frequent_10) 
common_freq <- numeric_annots[as.logical(binarized_maf_common),] %>% select(-MAFbin_lowfreq_1,-MAFbin_lowfreq_2,-MAFbin_lowfreq_3,-MAFbin_lowfreq_4,-MAFbin_lowfreq_5,-MAFbin_lowfreq_6,-MAFbin_lowfreq_7,-MAFbin_lowfreq_8,-MAFbin_lowfreq_9,-MAFbin_lowfreq_10,-MAFbin_frequent_1,-MAFbin_frequent_2,-MAFbin_frequent_3,-MAFbin_frequent_4,-MAFbin_frequent_5,-MAFbin_frequent_6,-MAFbin_frequent_7,-MAFbin_frequent_8,-MAFbin_frequent_9,-MAFbin_frequent_10) 

#Coding
coding_tracks_c <- data.frame(Common = common_freq$Coding_UCSC_common > 0, LowFreq = common_freq$Coding_UCSC_lowfreq > 0)
coding_indices_c <- coding_tracks_c$Common + coding_tracks_c$LowFreq
coding_tracks_l <- data.frame(Common = low_freq$Coding_UCSC_common > 0, LowFreq = low_freq$Coding_UCSC_lowfreq > 0)
coding_indices_l <- coding_tracks_l$Common + coding_tracks_l$LowFreq
coding_low <- low_freq[as.logical(coding_indices_l),] %>% select(-Coding_UCSC_lowfreq, -Coding_UCSC_common)
coding_common <- common_freq[as.logical(coding_indices_c),] %>% select(-Coding_UCSC_lowfreq, -Coding_UCSC_common)
noncoding_low <- low_freq[as.logical(1-unlist(coding_indices_l)),] %>% select(-Coding_UCSC_lowfreq, -Coding_UCSC_common)
noncoding_common <- common_freq[as.logical(1-unlist(coding_indices_c)),] %>% select(-Coding_UCSC_lowfreq, -Coding_UCSC_common)

print("Check: are we covering all the annotations?")
tot <- dim(coding_low)[1]/dim(numeric_annots)[1] + dim(coding_common)[1]/dim(numeric_annots)[1] + dim(noncoding_low)[1]/dim(numeric_annots)[1] + dim(noncoding_common)[1]/dim(numeric_annots)[1]
print(tot)

dim(coding_low)
dim(coding_common)
dim(noncoding_low)
dim(noncoding_common)
if(tot != 1)
{
  print("We may have an error! Do not proceed.")
}


#Remove 0 variance columns
output_dir <- "./"
coding_low <- varClear(coding_low)
coding_common <- varClear(coding_common)
noncoding_low <- varClear(noncoding_low)
noncoding_common <- varClear(noncoding_common)
write_tsv(coding_low, paste0(output_dir, "/coding_lowMAF.annots"))
write_tsv(coding_common, paste0(output_dir, "/coding_commonMAF.annots"))
write_tsv(noncoding_low, paste0(output_dir, "/noncoding_lowMAF.annots"))
write_tsv(noncoding_common, paste0(output_dir, "/noncoding_commonMAF.annots"))
