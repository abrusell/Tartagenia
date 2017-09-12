#!/usr/bin/env python
import re, sys, os, gzip

try:
	files = sys.argv
	vcfile = files[1]
	txt = files[2]
except:
	print "usage: python gnomAD_wgsa_annotator.py file.vcf wgsa_javaclass_output.txt"
	sys.exit(1)

txt_file = open(txt, "r")

diz = {}

for line in txt_file.readlines():
	spline = re.split("\t", line.rstrip())
	if spline[0] == "#CHROM":
		for i in range(1,len(spline)):
			if spline[i].startswith("gnomAD"):
				startCol = int(i)
				break
	else:
		ID = spline[0]+"_"+spline[1]+"_"+spline[3]+"_"+spline[4]
		if spline[startCol] == ".":
			next
		else:
			diz[ID]=spline[startCol:len(spline)]
#			print ID, diz[ID]

txt_file.close()


if vcfile.endswith(".gz"):
    vcf = gzip.open(vcfile,"r")
    out = open(vcfile[:-7]+"_gnomAD.vcf", "w")
else:
    vcf = open(vcfile,"r")
    out = open(vcfile[:-4]+"_gnomAD.vcf", "w")


for vcfline in vcf.readlines():
    if vcfline.startswith("#"):
        out.write(vcfline)
    else:
        splitvcfline = re.split("\t", vcfline.rstrip())
        var = splitvcfline[0]+"_"+splitvcfline[1]+"_"+splitvcfline[3]+"_"+splitvcfline[4]
        if var in diz:
        	out.write(splitvcfline[0]+"\t"+splitvcfline[1]+"\t"+splitvcfline[2]+"\t"+splitvcfline[3]+"\t"+splitvcfline[4]+"\t"+splitvcfline[5]+"\t"+splitvcfline[6]+"\t"+splitvcfline[7]+";dbNSFP_gnomAD_AF="+diz[var][2]+";dbNSFP_gnomAD_AFR_AF="+diz[var][5]+";dbNSFP_gnomAD_AMR_AF="+diz[var][8]+";dbNSFP_gnomAD_ASJ_AF="+diz[var][11]+";dbNSFP_gnomAD_EAS_AF="+diz[var][14]+";dbNSFP_gnomAD_FIN_AF="+diz[var][17]+";dbNSFP_gnomAD_NFE_AF="+diz[var][20]+";dbNSFP_gnomAD_SAS_AF="+diz[var][23]+";dbNSFP_gnomAD_OTH_AF="+diz[var][26])
		for i in range(8,len(splitvcfline)):
			out.write("\t"+splitvcfline[i])
		out.write("\n")
        else:
		out.write(vcfline)
vcf.close()
out.close()

