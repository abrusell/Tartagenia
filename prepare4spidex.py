#!/usr/bin/env python
import re, sys, os, gzip

try:
    files = sys.argv
    vcfile = files[1]
    id = files[2]
    outdir = files[3]
except:
    print "usage: python prepare4spidex.py file.vcf ID /path/to/outdir"
    sys.exit(1)

if vcfile.endswith(".gz"):
    vcf = gzip.open(vcfile,"r")
    filename = vcfile[:-7]
    out = open(vcfile[:-7]+"_spidex_job", "w")
else:
    vcf = open(vcfile,"r")
    filename = vcfile[:-4]
    out = open(vcfile[:-4]+"_spidex_job", "w")

out.write("#!/bin/bash\n#PBS -l select=1:ncpus=1:mem=32GB\n#PBS -N spidex_prep\n#PBS -W group_list=SIGU_Tarta17\n#PBS -A SIGU_Tarta17\n#PBS -q parallel\n#PBS -l walltime=24:00:00\n\numask 0002\nspidex=/pico/work/IscrC_FoRWArDS_1/db/spidex_public_noncommercial_v1.0/spidex_public_noncommercial_v1_0.tab.gz\nmodule load tabix\n\n")

for vcfline in vcf.readlines():
	if vcfline.startswith("##"):
		next
	elif vcfline.startswith("#CHROM"):
		splithead = re.split("\t", vcfline.rstrip())
		for i in range(0, len(splithead)):
			if splithead[i] == id:
                		genofield = int(i)
	else:
		spline = re.split("\t", vcfline)
		if (spline[6] == "PASS") or (spline[6] == "SnpCluster"):
			if (str(spline[genofield][0:3]) != "./.") and (str(spline[genofield][0:3]) != "0/0"):
				if ((len(spline[3]) == 1) and (len(spline[4]) == 1)):
					out.write("tabix $spidex "+str(spline[0])+":"+str(spline[1])+"-"+str(spline[1])+" | grep \"\t"+str(spline[3])+"\t"+str(spline[4])+"\t\" >> "+filename+"_spidex_output.txt\n")
																				# prima: >> "+str(outdir)+"/"+id+"_spidex_output.txt\n")
# l'output si chiamera' come il vcf + "_spidex_output.txt" per distinguere tra quello normale e quello di muTect2
vcf.close()
out.close()
