import os
from subprocess import check_call
import sys
import argparse

def whichCol(file, annot):
    i = 1
    t = None
    with open(file, 'r') as istream:
        for line in istream:
            line = line.strip()
            if line == annot:
                print("Found it!",i)
                return i
            else:
                print(line)
                i+=1
    if i == 2:
        tab = line.split()
        t= tab.index(annot)    
    if not t:
        print("Unable to find match for", annot, ", please try again")
        sys.exit()
    return t
parser = argparse.ArgumentParser()
parser.add_argument("--search", help = "Specify the list of snps you'd like search for in a single file with chr:pos:ref:alt ")
parser.add_argument("--output", help = "Specify where to write the output")
parser.add_argument("--database", help = "Specify if you want LDSC annots (for CADD, use  /work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/CADD-scripts-master/src/scripts/extract_scored.py", default = "LDSC")
parser.add_argument("--annot", help = "Specify which annotation you want to pull out")
args = parser.parse_args()

if args.database == "LDSC": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/head.tmp"
elif args.database == "DANN": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/dann_header.txt" 
#elif args.database == "LD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ld_scores_only/ldscore_header.txt"
elif args.database == "LD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ldsc_ukbb_3.2020/baselineLF_v2.2.UKB/ldheader.txt"
elif args.database == "CADD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/imputed_header.tsv"
elif args.database == "genocanyon": header = "chr\tpos\tresult"
elif args.database == "OC" or args.database == "OpenCravat": header = "NOT YET IMPLEMENTED"
else:
    print("Specified database doesn't exist, program will terminate")
    sys.exit()
    #header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/head_lines.tsv"

match_dex = 1
search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/imputed_cadd_annots_ukbb.tsv"
if args.database == "LD":
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ldsc_ukbb_3.2020/baselineLF_v2.2.UKB/weights.UKB." + str(args.chr) + ".l2.ldscore"
if args.database == "DANN":
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/DANN_whole_genome_SNVs_hg19.tsv"
if args.database == "LDSC":
    match_dex = 2
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/baselineLF." + str(args.chr) + ".annot"
if args.database == "genocanyon":
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/genocanyon_annotations/GenoCanyon_Chr" + str(args.chr) + "_hg19/"


i = whichCol(header, args.annot)
if args.database == "CADD":
  command = """awk '(FNR == NR) {arr[$1];next} ($2":"$3":"$4":"$5 in arr) {print $2":"$3":"$4":"$5"\t" $ """ + str(i) + "}' " + args.search + " " + search_db + " > " + args.output
  check_call(command, shell = True)
elif args.database == "LDSC:":
  print("Warning- there are likely multiple columns for many annotations in LDSC (low vs high freq).")
  for filename in os.listdir("/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots"):
    if filename.endswith(".annot"):
      f_name = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/" + filename
      command = """awk '(FNR == NR) {arr[$1];next} ($3 in arr) {print $3"\t" $ """ + str(i) + "}' " + args.search + " " + f_name + " >> " + args.output
      check_call(command, shell = True)
else:
  print("We have not implemented this for other databases, sorry.")
