#!/bin/env/python
import argparse

def writeQueryFile(l, c):
    with open(c+"_queries.q", 'w') as ostream:
        for q in l:
            ostream.write(q + '\n')
    return c+"_queries.q"

parser = argparse.ArgumentParser()
parser.add_argument("fin", help = "Specify the list of genomic locations to extract...")
parser.add_argument("out", help = "Specify where the results should go to", default = "cadd_annot_data.tsv")
args = parser.parse_args()
chr_dict = dict()
with open(args.fin,'r',) as istream:
    for line in istream:
        snp = line.strip()
        chr = snp.split(":")[0]
        nucl = snp.split(":")[1]
        if chr not in chr_dict:
            chr_dict[chr] = list()
        chr_dict[chr].append(nucl)

print("set -e")
for chr in chr_dict:
    fname = writeQueryFile(chr_dict[chr],chr)
    print("grep -a -F -f " + fname + " /work-zfs/abattle4/lab_data/genomic_annotation_data/hg19/cadd_annotations/chr" + str(chr) + "_cadd_annots.tsv >> "  + args.out)
print("rm *_queries.q")
