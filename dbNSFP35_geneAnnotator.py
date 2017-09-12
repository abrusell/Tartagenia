#!/usr/bin/env python
import re, sys, os, gzip

try:
    files = sys.argv
    vcfile = files[1]
except:
    print "usage: python dbNSFP35_geneAnnotator.py file.vcf"
    sys.exit(1)

dbNSFPfile = open("/pico/work/IscrC_FoRWArDS_1/db/dbNSFP3.5/dbNSFP3.5_gene_comma", "r")

diz = {}
for line in dbNSFPfile.readlines():
    if "Gene_name" in line:
        next
    else:
        spline = re.split("\t", line.rstrip())
	gene = "_"+spline[0]+"_"
        diz[gene]=";Gene_old_names="+spline[3]+";Gene_other_names="+spline[4]+";Entrez_gene_id="+spline[7]+";CCDS_id="+spline[8]+";Refseq_id="+spline[9]+";ucsc_id="+spline[10]+";Gene_full_name="+spline[12]+";Pathway(ConsensusPathDB)="+spline[16]+";Function_description(Uniprot)="+spline[19]+";Disease_description(Uniprot)="+spline[20]+";MIM_disease="+spline[22]+";GO_biological_process="+spline[24]+";Tissue_specificity(uniprot)="+spline[27]+";RVIS_ExAC="+spline[39]+";RVIS_percentile_ExAC="+spline[40]+";ExAC_pLI="+spline[42]+";ExAC_pRec="+spline[43]+";GDI="+spline[55]+";Gene_damage_prediction="+spline[57]+";SORVA_LOF_MAF0.001_HetOrHom="+spline[70]+";SORVA_LOF_MAF0.001_HomOrCompHet="+spline[71]+";MGI_mouse_phenotype="+spline[78]+";ZFIN_zebrafish_phenotype_tag="+spline[82]

dbNSFPfile.close()

if vcfile.endswith(".gz"):
    stdout_handle= os.popen("zgrep ^## "+vcfile, "r")
    header1 = stdout_handle.read()
    stdout_handle= os.popen("zgrep ^# "+vcfile+" | grep -v ^##", "r")
    header2 = stdout_handle.read()
    vcf = gzip.open(vcfile,"r")
    out = open(vcfile[:-7]+".dbNSFP3_gene.vcf", "w")
else:
    stdout_handle= os.popen("grep ^## "+vcfile, "r")
    header1 = stdout_handle.read()
    stdout_handle= os.popen("grep ^# "+vcfile+" | grep -v ^##", "r")
    header2 = stdout_handle.read()
    vcf = open(vcfile,"r")
    out = open(vcfile[:-4]+".dbNSFP3_gene.vcf", "w")

geneKey = ""
found = ""

out.write(header1+"##dbNSFP_gene_file="+dbNSFPfile.name+"\n##dbNSFP_gene annotation fields:\n##INFO=<ID=Gene_old_names,Number=.,Type=String,Description=\"Old gene symbol (from HGNC)\">\n##INFO=<ID=Gene_other_names,Number=.,Type=String,Description=\"Other gene names (from HGNC)\">\n##INFO=<ID=Entrez_gene_id,Number=.,Type=String,Description=\"Entrez gene id (from HGNC)\">\n##INFO=<ID=CCDS_id,Number=.,Type=String,Description=\"CCDS id (from HGNC)\">\n##INFO=<ID=Refseq_id,Number=.,Type=String,Description=\"Refseq gene id (from HGNC)\">\n##INFO=<ID=ucsc_id,Number=.,Type=String,Description=\"UCSC gene id (from HGNC)\">\n##INFO=<ID=Gene_full_name,Number=.,Type=String,Description=\"Gene full name (from HGNC)\">\n##INFO=<ID=Pathway(ConsensusPathDB),Number=1,Type=String,Description=\"Pathway(s) the gene belongs to (from ConsensusPathDB)\">\n##INFO=<ID=Function_description(Uniprot),Number=.,Type=String,Description=\"Function description of the gene (from Uniprot)\">\n##INFO=<ID=MIM_disease,Number=.,Type=String,Description=\"MIM_disease: MIM disease name(s) with MIM id(s) between square brackets (from Uniprot)\">\n##INFO=<ID=Disease_description,Number=.,Type=String,Description=\"Disease(s) the gene caused or associated with (from Uniprot)\">\n##INFO=<ID=GO_biological_process,Number=.,Type=String,Description=\"GO terms for biological process\">\n##INFO=<ID=Tissue_specificity(Uniprot),Number=.,Type=String,Description=\"Tissue specificity description from Uniprot\">\n##INFO=<ID=RVIS_ExAC,Number=1,Type=Float,Description=\"ExAC-based RVIS; setting 'common' MAF filter at 0.05% in at least one of the six individuals ethnic strata from ExAC. Cited from RVIS document\">\n##INFO=<ID=RVIS_percentile_ExAC,Number=1,Type=Float,Description=\"Genome-Wide percentile for the new ExAC-based RVIS; setting 'common' MAF filter at 0.05% in at least one of the six individual ethnic strata from ExAC. Cited from RVIS document\">\n##INFO=<ID=ExAC_pLi,Number=1,Type=Float,Description=\"the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants). Based on ExAC r0.3 data\">\n##INFO=<ID=ExAC_pRec,Number=1,Type=Float,Description=\"the probability of being intolerant of homozygous, but not heterozygous lof variants. Based on ExAC r0.3 data\">\n##INFO=<ID=GDI,Number=1,Type=Float,Description=\"Gene damage index score, 'a genome-wide, gene-level metric of the mutational damage that has accumulated in the general population' from doi: 10.1073/pnas.1518646112. The higher the score the less likely the gene is to be responsible for monogenic diseases\">\n##INFO=<ID=Gene_damage_prediction,Number=1,Type=String,Description=\"(all disease-causing genes): gene damage prediction (low/medium/high) by GDI\">\n##INFO=<ID=SORVA_LOF_MAF0.001_HetOrHom,Number=1,Type=Float,Description=\"the fraction of individuals in the 1000 Genomes Project data (N=2504) who are either Compound Heterozygote or Homozygote of LOF SNVs whose MAF<0.005. doi: 10.1101/103218\">\n##INFO=<ID=SORVA_LOF_MAF0.001_HomOrCompHetSORVA_LOF_MAF0.001_HomOrCompHet,Number=1,Type=Float,Description=\"the fraction of individuals in the 1000 Genomes Project data (N=2504) who are either Compound Heterozygote or Homozygote of LOF or missense SNVs whose MAF<0.001. doi: 10.1101/103218\">\n##INFO=<ID=MGI_mouse_phenotype,Number=.,Type=String,Description=\"Phenotype description for the homolog mouse gene from MGI\">\n##INFO=<ID=ZFIN_zebrafish_phenotype_tag,Number=.,Type=String,Description=\"Phenotype tag for the homolog zebrafish gene from ZFIN\">\n"+header2)
for vcfline in vcf.readlines():
    if vcfline.startswith("#"):
        next
    else:
        spline = re.split("\t", vcfline.rstrip())
        gene2=re.search('geneName=(.+?);', vcfline.rstrip())
        if gene2:
            found = gene2.group(1)
	else:
	    found = ""
        geneKey = "_"+found+"_"
        if geneKey in diz:
            out.write(spline[0]+"\t"+spline[1]+"\t"+spline[2]+"\t"+spline[3]+"\t"+spline[4]+"\t"+spline[5]+"\t"+spline[6]+"\t"+spline[7]+diz[geneKey])
	    for i in range(8,len(spline)):
	        out.write("\t"+spline[i])
	    out.write("\n")
        else:
            out.write(vcfline)
vcf.close()
out.close()

