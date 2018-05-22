#!/usr/bin/env python
import re, sys, os

try:
    files = sys.argv
    txtfile = files[1]
except:
    print "usage: python dbNSFP35_geneAnnotator_macaron.py macaron_output.txt"
    sys.exit(1)

dbNSFPfile = open("/gpfs/work/SIGU_Tarta18_0/db/dbNSFP3.5/dbNSFP3.5_gene_comma", "r")

diz = {}
for line in dbNSFPfile.readlines():
    if "Gene_name" in line:
        next
    else:
        spline = re.split("\t", line.rstrip())
	gene = "_"+spline[0]+"_"
        diz[gene]="\t"+spline[19]+"\t"+spline[20]+"\t"+spline[22]+"\t"+spline[39]+"\t"+spline[40]+"\t"+spline[42]+"\t"+spline[43]+"\t"+spline[55]

dbNSFPfile.close()

txt = open(txtfile,"r")
out = open(txtfile[:-4]+"_dbNSFP3_gene.tsv", "w")

for txtline in txt.readlines():
    if txtline.startswith("CHROM"):
        out.write(txtline.rstrip()+"\tFunction_description(Uniprot)\tDisease_description\tMIM_disease=\tRVIS_ExAC\tRVIS_percentile_ExAC\tExAC_pLI\tExAC_pRec\tGDI\n")
    else:
        spline2 = re.split("\t", txtline.rstrip())
        gene2="_"+spline2[5]+"_"
        if gene2 in diz:
            out.write(txtline.rstrip()+diz[gene2]+"\n")
        else:
            out.write(txtline)
txt.close()
out.close()

