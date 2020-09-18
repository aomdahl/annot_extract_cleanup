library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)


output_dir =  ""
annot_names <- scan("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/header.txt", what = character())
max_annots <- fread("/work-zfs/abattle4/ashton/prs_dev/scratch/annotation_curation/redoing_annots_may/full_annotations_list/all_annotations.250-kb-win_0.5_r2.tsv")
colnames(max_annots) <- annot_names
indx <- grepl('.na', annot_names) #Special ones imputed with uniqe cases, to be done separately
max_annots <- max_annots[,!indx, with = F]
print("Annotations read in")

#Helper functions
clearMissing <- function(dat_in)
{
  remove_cols <- apply(dat_in, 2, function(col)sum(is.na(col))/length(col)) <= 0.3
  dat_in <- dat_in[,remove_cols, with = F]
  #Remove columns with no variance
  #Or here, keep them if the largest and smallest values are different
  remove_cols <- apply(dat_in, 2, function(col) max(col, na.rm = T) != min(col, na.rm = T))
  dat_in <- dat_in[,remove_cols, with = F]
  dat_in
}

#Return it sorted by 
selectNumeric <- function(dat_in)
{
  #note that there may be some others in there that look numeric but we don't want, like 
  dat_in %>% arrange(F_ID) %>% select(-Chr, -Pos, -Ref, -Alt, -SNP) %>% select_if(is.numeric)
}
selectLabels <- function(dat_in)
{
  dat_in$Chr <- as.character(dat_in$Chr)
  dat_in$Pos <- as.character( dat_in$Pos)
  dat_in %>% arrange(F_ID) %>% select_if(is.character)
}
#Remove columns that have 0 variance
varClear <- function(tab)
{
  zv <- names(which(apply(tab, 2, var) == 0))
  tab %>% select(-zv)
}

max_annots <- clearMissing(max_annots) %>% drop_na() #the
numeric_annots <- selectNumeric(max_annots)
numeric_annots <- varClear(numeric_annots)
label_annots <- selectLabels(max_annots)
write_tsv(numeric_annots, "./cleaned_up_numeric_annotations.tsv") 
write_tsv(label_annots, "./cleaned-up_annot_labels.tsv")
print("Annotations cleaned up (removed missing,0 variance, etc.")

