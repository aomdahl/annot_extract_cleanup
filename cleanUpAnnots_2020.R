#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyr)
library(data.table)
library(readr)
library(dplyr)
library(stringr)
library(Xmisc)

parser <- ArgumentParser$new()
parser$add_description("Script that cleans up input annotation file by removing high missingness columns, removing columns with no variance, and spliting into numeric and non-numeric annotations")
parser$add_argument("--input_annots", type = 'character', help = "File containing annotations in a tsv.")
parser$add_argument("--output", type = "character", help = "Output file path and handle (i.e. /this/location/filename")
parser$add_argument("--impute", type = "logical", help = "If you would like to impute missing data", default =F)
parser$add_argument("--r2_thresh", type = "numeric", help = "Specify an R2 threshold for annotation similarity. If 2 features with R2 above this threshold exist, one will be dropped based on the order it appears in the list. The default is R2 = 0.95", default = 0.95)
parser$add_argument("--na_thresh", type = "numeric", help = "Specify a threshold of % NAs for dropping an annotation column. Default is 0.1", default = 0.1)
parser$add_argument('--help',type='logical',action='store_true',help='Print the help page')
parser$helpme()
args <- parser$get_args()

output_dir =  args$output
max_annots <- fread(args$input_annots)
impute <- args$impute
threshold <- args$na_thresh

print("Annotations read in of size:")
print(paste(dim(max_annots)[1], dim(max_annots)[2]))
#Helper functions
#Get rid of columns with >10% NAs, and any columns where all values are the same.
clearMissing <- function(dat_in, na_thresh)
{
  keep_cols <- apply(dat_in, 2, function(col) (sum(is.na(col))/length(col)) < na_thresh)
  t <- dat_in[,keep_cols, with = F]

#Or here, keep them if the largest and smallest values are different
  keep_cols <- apply(t, 2, function(col) max(col, na.rm = T) != min(col, na.rm = T))
  t[,keep_cols, with = F]
}

#Return it sorted by 
selectNumeric <- function(dat_in)
{
  #note that there may be some others in there that look numeric but we don't want, like 
  dat_in %>% arrange(F_ID) %>% select_if(is.numeric)
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
#alternative
    #Remove columns with no variance
}


#Remove annotations that are highly correlated (corr > 0.95)
highCorClear <- function(annots, drop_thresh, outdir)
{
    r2 <- cor(annots)^2
    name_list <- names(annots)
    dropped_names <- c()
    drop_sto <- c()
    for ( i in 1:length(name_list))
    {
      if (!is.na(name_list[i]))
      {
        rem <- which(r2[,i] > drop_thresh & r2[,i] < 1)
        if(length(rem) > 0)
        {
          name_list[rem] <- NA
        dropped_names <- c(dropped_names, names(rem))
        drop_sto <- c(drop_sto, paste("Dropped", names(rem), "kept", name_list[i]))
        }
          
      }
    }
    final_drops <- unique(dropped_names)
    #Write out the drops
    sink(paste0(outdir, "_dropped_features.txt"))
    print(drop_sto)
    sink()
    return(annots %>% select(-final_drops))

}
print("Removing features with high missingness and no variance")
f_annots <- clearMissing(max_annots, threshold) 
print("New size of data after removing columns with high missingness and with no variance")
print(paste(dim(f_annots)[1], dim(f_annots)[2]))
#How many NAs per column?
#missing_count <- apply(t, 2, function(x) sum(is.na(x)))/5455
if (impute)
{
        #Some imputetation scheme.
        print("Do nothing yet.")
}else
{
   f_annots <- drop_na(f_annots)
}

#Data for features dropped because of NAs
print("New size of data after dropping NAs.")
print(paste(dim(f_annots)[1], dim(f_annots)[2]))
dropped <- names(max_annots)[!(names(max_annots) %in% names(f_annots))]
prop_nas <- apply(max_annots[, ..dropped], 2, function(x) sum(is.na(x)))/nrow(max_annots)
dropped <- data.frame("dropped_annots" = dropped, "proportion_of_na" = prop_nas)

numeric_annots <- selectNumeric(f_annots) %>% select(-Pos, -Chr)
if("biogrid:id" %in% names(numeric_annots))
{
numeric_annots <- select(numeric_annots, -`biogrid:id`)
}
numeric_annots <- varClear(numeric_annots) %>% highCorClear(., args$r2_thresh, args$output)
label_annots <- selectLabels(f_annots)
write_tsv(numeric_annots, paste0(output_dir, "_numeric_annotations.tsv"))
write_tsv(label_annots, paste0(output_dir, "_annotation_labels.tsv"))

write_tsv(dropped, paste0(output_dir, "_removed_annots.tsv"))
print(paste0("Dropped annotations written out to ", output_dir, "_removed_annots.tsv"))
#print("Dropped variants written out to ", output_dir, "_removed_vars.tsv")
print("Annotations cleaned up and split into numeric/text ones. (removed missing, removed columns with 0 variance, removed highly correlated annotations)")
