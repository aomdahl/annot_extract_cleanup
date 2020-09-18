import os
from subprocess import check_call
import sys
import argparse

def writeOut(retdict, diro, db):
    #with open(h, 'r') as istream:
    #    f = istream.readline().strip()
    #dat = f.split()

    nofind = list()
    with open(diro, 'a') as ostream:
        for entry in retdict:
            if retdict[entry] == None:
                nofind.append(entry)
            else:
                ostream.write(retdict[entry] + '\n')
    print("Unable to find ", str(len(nofind)), "queries. These will be written out to unfound_" + diro +".txt")
    fname = 'unfound_' + diro + '.txt'
    with open(fname, 'a') as ostream:
        #Write the header
        #ostream.write('\t'.join(dat) + '\n')
        for entry in nofind:
            ostream.write(entry + '\n')


parser = argparse.ArgumentParser()
parser.add_argument("--build_queries", help = "Build chromosome based queries using a list of queries to search.")
parser.add_argument("--chr", help = "Specify which chromosome to search for")
parser.add_argument("--out", help = "Specify where the results should go to", default = "query_hits.tsv")
parser.add_argument("--search", help = "Specify the list of snps you'd like to search for in a VCF file, where column 3 is the ID ")
parser.add_argument("--database", help = "Specify if you want LDSC annots (for CADD, use  /work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/CADD-scripts-master/src/scripts/extract_scored.py", default = "LDSC")
parser.add_argument("--printo", help = "Specify this if you wish to print directly to console instead of writing out to file", action = "store_true")
parser.add_argument("--separate_files", help = "Specify this argument if you'd like to keep each file separate by chromosome, default is to merge", action = "store_true")
args = parser.parse_args()

if args.database == "LDSC": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/head.tmp"
elif args.database == "DANN": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/dann_header.txt" 
#elif args.database == "LD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ld_scores_only/ldscore_header.txt"
elif args.database == "LD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ldsc_ukbb_3.2020/baselineLF_v2.2.UKB/ldheader.txt"
elif args.database == "CADD": header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/imputed_header.tsv"
elif args.database == "genocanyon": header = "chr\tpos\tresult"
else:
    print("Specified database doesn't exist, program will terminate")
    sys.exit()
    #header = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/head_lines.tsv"
if args.database == "CADD" and args.chr:
    print("CADD is not chr specific, chr argument will be ignored")
if args.build_queries:
    if args.database == "CADD":
        print("Current version of CADD is not broken down by chromosme. Please run the following command")
        print("python annot_extraction -py --out outfile --database CADD --search search_file") 
    command = '''FIN="''' + str(args.build_queries) + '''";for i in {1..22}; do grep -e "^${i}:" $FIN > chr${i}_query.tmp ; done'''
    check_call(command, shell = True)
    script_spot = __file__
    for i in range(1,23):
        app = ""
        print("python " + script_spot + " --chr " + str(i) + " --out chr" + str(i) + "db_hits.tsv " + app + " --database " + args.database + " --search " + "chr" + str(i) + "_query.tmp &")
    print("wait")
    print("cat " + header + " chr1db_hits.tsv > running")
    print("for i in {2..22}; do cat running chr${i}db_hits.tsv > tmp; mv tmp running;done")
    print("mv running chralldb.tsv")
    sys.exit()
#awk '(FNR == NR) {q[$1]; next} {split($3,d, ":")} (d[1]":"d[2] in q'

if args.search and not args.chr and args.database != "CADD" and args.database != "DANN" and args.database != "LD":
    print("Please specify a chromosome to search for. Program will terminate")
    sys.exit()
#Into python

match_dex = 1
#/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/cadd_extraction_approved/
search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/imputed_cadd_annots_ukbb.tsv"
#search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/cadd_extraction_approved/chr" + str(args.chr) + "_cadd_annots.tsv"

if args.database == "LD":
    #search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ld_scores_only/all_ldscores.tsv"
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/ldsc_ukbb_3.2020/baselineLF_v2.2.UKB/weights.UKB." + str(args.chr) + ".l2.ldscore"
if args.database == "DANN":
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/DANN_whole_genome_SNVs_hg19.tsv"
if args.database == "LDSC":
    match_dex = 2
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/ldsc_annotations/baselineLF_v2.2.UKB/annots/baselineLF." + str(args.chr) + ".annot"
if args.database == "genocanyon":
    search_db = "/work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/genocanyon_annotations/GenoCanyon_Chr" + str(args.chr) + "_hg19/"
#Read in the queries for it
q_in = dict()
with open(args.search) as istream:
    for q in istream:
        if q[0] == "#":
            continue
        q_in[q.strip()] = None
            #q_in[q.strip().split()[2]] = None
if args.database == "genocanyon":
    search_list = [search_db + x for x in os.listdir(search_db)]
else: 
    search_list = [search_db]
for db in search_list:
    with open(db) as istream:
        line_counter = 0
        dups = set()
        for line in istream:
            line_counter += 1
            line = line.strip()
            dat = line.split()
            if args.database == "CADD": #first arguments in file are y chr pos ref alt
                #if dat[0] != args.chr and line_counter > 1:
                    #print("Invalid dataline, chr", args.chr, dat[0:3])
                    #print("Line num", line_counter)
                    #continue #just continue on
                try:
                    search_id = dat[1] + ":" + dat[2]
                except IndexError:
                    print("Unexpected error in id, where")
                    print(dat[0:10])
                    continue
            if args.database == "LD":
                #chr = dat[1].split(":")
                try:
                    search_id = (dat[0] + ':' + dat[2])
                except IndexError:
                    continue
            if args.database == "DANN":
                chr = dat
                try:
                    search_id = (chr[0] + ':' + chr[1])
                except IndexError:
                    continue
            if args.database == "LDSC":
                if dat[0] == "CHR":
                    continue
                chr = dat[match_dex].split(":")
                try:
                    search_id = dat[match_dex]
                    #search_id = (chr[0] + ':' + chr[1])
                except IndexError:
                    continue
            if args.database == "genocanyon":
                if dat[0] == "chr":
                    continue
                try:
                    search_id = (dat[0] +":"+ dat[1])[3:]
                except IndexError:
                    continue
            if search_id in q_in:
                if q_in[search_id] == None and search_id not in dups: #i.e. we haven't encountered it yet.
                    q_in[search_id] = line
                    if args.printo: 
                        print(line)
                else: #something ha been put there already
                    dups.add(search_id)
                    print("Duplicate entries found for ", search_id, ". We will skip this.")
                    #del q_in[search_id]
                    #no longer remove the entry
                    continue #keep only the original entry.
print("Number of duplicate entries (skipped):", len(dups))
                
if not args.printo: 
    writeOut(q_in, args.out, args.database)
