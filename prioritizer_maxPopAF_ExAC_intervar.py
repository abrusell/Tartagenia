#!/usr/bin/env python
import re, sys, os, gzip

try:
    files = sys.argv
    phenolist = files[1]
    vcfile = files[2]
except:
    print "usage: python prioritizer.py phenotype_list.txt file.vcf"
    sys.exit(1)


if vcfile.endswith(".gz"):
    stdout_handle= os.popen("zgrep ^# "+vcfile+" | grep -v ^##", "r")
    header = stdout_handle.read()
    vcf = gzip.open(vcfile,"r")
    out = open(vcfile[:-7]+".scored.tsv", "w")
else:
    stdout_handle= os.popen("grep ^# "+vcfile+" | grep -v ^##", "r")
    header = stdout_handle.read()
    vcf = open(vcfile,"r")
    out = open(vcfile[:-4]+".scored.tsv", "w")

with open(phenolist) as f:
	pheno_terms = f.read().splitlines()

out.write("##Phenolyzer terms used: "+", ".join(pheno_terms)+"\n")
out.write(header.rstrip()+"\tVar_nt_pos\tVar_aa_pos\tEffect\tgeneName\tDisease_description\tGene_full_name\tCADD_phred\tdbNSFP_MetaSVM_score\tRECURRENCE_IN_OUR_DB\tPhenolyzer_score\tFinal_score\tInterVar\tImpact\tFeature_id\tANN\tdbSNP.GMAF\tGO_biological_process\tdbNSFP_maxPop_AF_ExAC_AF\tRVIS_ExAC_0.05_percentile\tGDI\tWarning\n")

found_exac			= ""
found_db			= ""
found_impact			= ""
found_cadd			= ""
found_svm			= ""
found_rvis			= ""
found_gdi			= ""
found_phlyz			= ""
found_var_nt_pos		= ""
found_var_aa_pos		= ""
found_maineffect		= ""
found_genename			= ""
found_feature_id		= ""
found_ann			= ""
found_dbsnp			= ""
found_disease_desc		= ""
found_gene_full_name		= ""
found_go_biol_proc		= ""
found_intervar			= ""

for vcfline in vcf.readlines():
	if vcfline.startswith("#"):
		next
	else:
		exac_score		= 0
		db_score		= 0
		impact_score		= 0
		Var_effect_score	= 0
		Gene_predict_score	= 0
		Phenolyzer_score	= 0
		Final_score		= 0
		warn			= 0

		spline = re.split("\t", vcfline.rstrip())
		var_nt_pos=re.search('Var_nt_pos=(.+?);', vcfline.rstrip())
		if var_nt_pos:
			found_var_nt_pos = var_nt_pos.group(1)
			if found_var_nt_pos.startswith(";"):
				found_var_nt_pos = "."
		else:
			found_var_nt_pos = "."

		var_aa_pos=re.search('Var_aa_pos=(.+?);', vcfline.rstrip())
		if var_aa_pos:
			found_var_aa_pos = var_aa_pos.group(1)
			if found_var_aa_pos.startswith(";"):
				found_var_aa_pos = "."
		else:
			found_var_aa_pos = "."

		maineffect=re.search('MainEffect=(.+?);', vcfline.rstrip())
		if maineffect:
			found_maineffect = maineffect.group(1)
			if found_maineffect.startswith(";"):
				found_maineffect = "."
		else:
			found_maineffect = "."

		genename=re.search('geneName=(.+?);', vcfline.rstrip())
		if genename:
			found_genename = genename.group(1)
			if found_genename.startswith(";"):
				found_genename = "."
		else:
			found_genename = "."

		feature_id=re.search('Feature_id=(.+?);', vcfline.rstrip())
		if feature_id:
			found_feature_id = feature_id.group(1)
			if found_feature_id.startswith(";"):
				found_feature_id = "."
		else:
			found_feature_id = "."

		ann=re.search('ANN=(.+?);', vcfline.rstrip())
		if ann:
			found_ann = ann.group(1)
			if found_ann.startswith(";"):
				found_ann = "."
		else:
			found_ann = "."

		dbsnp=re.search('dbSNP1..\.GMAF=(.+?);', vcfline.rstrip())
		if dbsnp:
			found_dbsnp = dbsnp.group(1)
			if found_dbsnp.startswith(";"):
				found_dbsnp = "."
		else:
			found_dbsnp = "."

		disease_desc=re.search('Disease_description=(.+?);', vcfline.rstrip())
		if disease_desc:
			found_disease_desc = disease_desc.group(1)
			if found_disease_desc.startswith(";"):
				found_disease_desc = "."
		else:
			found_disease_desc = "."

		gene_full_name=re.search('Gene_full_name=(.+?);', vcfline.rstrip())
		if gene_full_name:
			found_gene_full_name = gene_full_name.group(1)
			if found_gene_full_name.startswith(";"):
				found_gene_full_name = "."
		else:
			found_gene_full_name = "."

		go_biol_proc=re.search('GO_biological_process=(.+?);', vcfline.rstrip())
		if go_biol_proc:
			found_go_biol_proc = go_biol_proc.group(1)
			if found_go_biol_proc.startswith(";"):
				found_go_biol_proc = "."
		else:
			found_go_biol_proc = "."

                afr_exac=re.search('dbNSFP_ExAC_AFR_AF=(.+?);', vcfline.rstrip())
                if afr_exac:
                        found_afr_exac = afr_exac.group(1)
                        if "," in found_afr_exac:
                                split_found_afr_exac = re.split(",", found_afr_exac)
                                found_afr_exac = split_found_afr_exac[0]
                                warn = 1
                else:
                        found_afr_exac = "."

		if found_afr_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC = []
			all_ExAC.append(found_afr_exac)

                amr_exac=re.search('dbNSFP_ExAC_AMR_AF=(.+?);', vcfline.rstrip())
                if amr_exac:
                        found_amr_exac = amr_exac.group(1)
                        if "," in found_amr_exac:
                                split_found_amr_exac = re.split(",", found_amr_exac)
                                found_amr_exac = split_found_amr_exac[0]
                                warn = 1
                else:
                        found_amr_exac = "."
		if found_amr_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC.append(found_amr_exac)

                eas_exac=re.search('dbNSFP_ExAC_EAS_AF=(.+?);', vcfline.rstrip())
                if eas_exac:
                        found_eas_exac = eas_exac.group(1)
                        if "," in found_eas_exac:
                                split_found_eas_exac = re.split(",", found_eas_exac)
                                found_eas_exac = split_found_eas_exac[0]
                                warn = 1
                else:
                        found_eas_exac = "."
		if found_eas_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC.append(found_eas_exac)

                fin_exac=re.search('dbNSFP_ExAC_FIN_AF=(.+?);', vcfline.rstrip())
                if fin_exac:
                        found_fin_exac = fin_exac.group(1)
                        if "," in found_fin_exac:
                                split_found_fin_exac = re.split(",", found_fin_exac)
                                found_fin_exac = split_found_fin_exac[0]
                                warn = 1
                else:
                        found_fin_exac = "."
		if found_fin_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC.append(found_fin_exac)

                nfe_exac=re.search('dbNSFP_ExAC_NFE_AF=(.+?);', vcfline.rstrip())
                if nfe_exac:
                        found_nfe_exac = nfe_exac.group(1)
                        if "," in found_nfe_exac:
                                split_found_nfe_exac = re.split(",", found_nfe_exac)
                                found_nfe_exac = split_found_nfe_exac[0]
                                warn = 1
                else:
                        found_nfe_exac = "."
		if found_nfe_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC.append(found_nfe_exac)

                sas_exac=re.search('dbNSFP_ExAC_SAS_AF=(.+?);', vcfline.rstrip())
                if sas_exac:
                        found_sas_exac = sas_exac.group(1)
                        if "," in found_sas_exac:
                                split_found_sas_exac = re.split(",", found_sas_exac)
                                found_sas_exac = split_found_sas_exac[0]
                                warn = 1
                else:
                        found_sas_exac = "."
		if found_sas_exac == ".":
			all_ExAC = "."
		else:
			all_ExAC.append(found_sas_exac)

		# convert each value in data_list to an integer using a list comprehension and then return the maximum value
		if not "." in all_ExAC:
			max_ExAC_AF = max([float(i) for i in all_ExAC])
		#print found_afr_exac,found_amr_exac,found_eas_exac,found_fin_exac,found_nfe_exac,found_sas_exac,max(all_ExAC)

		db=re.search('RECURRENCE_IN_OUR_DB=(.+?);', vcfline.rstrip())
		if db:
			found_db = db.group(1)
			subfound_db = re.split("/", found_db)
		else:
			found_db = "."

                impact=re.search('Impact=(.+?);', vcfline.rstrip())
                if impact:
                        found_impact = impact.group(1)
			if found_impact.startswith(";"):
				found_impact = "."
                else:
                        found_impact = "."

                cadd=re.search('CADD_phred=(.+?);', vcfline.rstrip())
                if cadd:
                        found_cadd = cadd.group(1)
			if "," in found_cadd:
				split_found_cadd = re.split(",", found_cadd)
				found_cadd = split_found_cadd[0]
				warn = 1
                else:
                        found_cadd = "."

                svm=re.search('dbNSFP_MetaSVM_score=(.+?);', vcfline.rstrip())
                if svm:
                        found_svm = svm.group(1)
			if "," in found_svm:
				split_found_svm = re.split(",", found_svm)
				found_svm = split_found_svm[0]
				warn = 1
                else:
                        found_svm = "."

                rvis=re.search('RVIS_ExAC_0\.05_percentile=(.+?);', vcfline.rstrip())
                if rvis:
                        found_rvis = rvis.group(1)
			if found_rvis.startswith(";"):
				found_rvis = "."
                else:
                        found_rvis = "."

                gdi=re.search('GDI=(.+?);', vcfline.rstrip())
                if gdi:
                        found_gdi = gdi.group(1)
			if found_gdi.startswith(";"):
				found_gdi = "."
                else:
                        found_gdi = "."

		phlyz=re.search('Phenolyzer_score=([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)', vcfline.rstrip())
		if phlyz:
			found_phlyz = phlyz.group(1)
			if str(found_phlyz).startswith(";"):
				found_phlyz = "."
		else:
			found_phlyz = "."

		intervar=re.search('InterVar=(.+?);', vcfline.rstrip())
		if intervar:
			found_intervar = intervar.group(1)
			if str(found_intervar).startswith(";"):
				found_intervar = "."
		else:
			found_intervar = "."

		if (str(all_ExAC) == ".") or (max_ExAC_AF < 0.005):
                        exac_score = 200
                elif (max_ExAC_AF >= 0.005) and (max_ExAC_AF < 0.01):
                        exac_score = 80
                elif (max_ExAC_AF >= 0.01) and (max_ExAC_AF < 0.03):
                        exac_score = 20
                else:
                        exac_score = 0

		if ( float(subfound_db[0]) == 0 ) or ( float(subfound_db[0]) == 1 ):
		        db_score = 200
		elif float(subfound_db[0])/float(subfound_db[1]) < 0.01:
		        db_score = 80
		elif float(subfound_db[0])/float(subfound_db[1]) < 0.03:
		        db_score = 20
		else:
	        	db_score = 0

		if str(found_impact) == "HIGH":
			impact_score = 200
		elif str(found_impact) == "MODERATE":
			impact_score = 80
		elif str(found_impact) == "LOW":
			impact_score = 20
		else:
			impact_score = 0

		if found_cadd == "." :
			if found_svm == ".":
				Var_effect_score = 40
			elif float(found_svm) > 0:
				Var_effect_score = 100
			else:
				Var_effect_score = 10
		elif float(found_cadd) > 15: 
			if found_svm == ".":
				Var_effect_score = 100
			elif float(found_svm) > 0:
				Var_effect_score = 100
			else:
				Var_effect_score = 40
		elif float(found_cadd) <= 15:
			if found_svm == ".":
				Var_effect_score = 10
			elif float(found_svm) > 0:
				Var_effect_score = 40
			else:
				Var_effect_score = 10

		if found_rvis == ".":
			if found_gdi == ".":
				Gene_predict_score = 40
			elif float(found_gdi) <= 13.84:
				Gene_predict_score = 100
			else:
				Gene_predict_score = 10
		elif float(found_rvis) < 25:
			if found_gdi == ".":
				Gene_predict_score = 100
			elif float(found_gdi) <= 13.84:
				Gene_predict_score = 100
			else:
				Gene_predict_score = 40
		elif float(found_rvis) >= 25:
			if found_gdi == ".":
				Gene_predict_score = 10
			elif float(found_gdi) <= 13.84:
				Gene_predict_score = 40
			else:
				Gene_predict_score = 10

		if found_phlyz == ".":
			Phenolyzer_score = 100
		elif float(found_phlyz) < 0.25:
			Phenolyzer_score = 100
		elif float(found_phlyz) < 0.5:
			Phenolyzer_score = 200
		elif float(found_phlyz) <= 1:
			Phenolyzer_score = 300
		else:
			Phenolyzer_score = 100

		Final_score = 0.5*float(exac_score)+0.5*float(db_score)+0.25*float(impact_score)+0.25*float(Var_effect_score)+0.25*float(Gene_predict_score)+0.25*float(Phenolyzer_score)
		if (float(subfound_db[0])/float(subfound_db[1]) >= 0.25):
			Final_score = Final_score - 50
		if "," in spline[4]:
			Final_score = Final_score - 20

#		print spline[0], spline[1], "exac: ",found_exac,"\tdb: ",str(float(subfound_db[0])/float(subfound_db[1])),"\timpact: ",found_impact,"\tcadd: ",found_cadd,"\tsvm: ",found_svm,"\trvis: ",found_rvis,"\tgdi: ",found_gdi,"\tphenolyzer score: ",found_phlyz
#		print "exac_score:\t",0.5*(float(exac_score))
#		print "db_score:\t",0.5*(float(db_score))
#		print "impact_score:\t",0.25*(float(impact_score))
#		print "Var_effect_score:\t",0.25*(float(Var_effect_score))
#		print "Gene_predict_score:\t",0.25*(float(Gene_predict_score))
#		print "Phenolyzer_score:\t",0.25*(float(Phenolyzer_score))
#		print spline[0]+"\t"+spline[1]+"\tFinal_score:\t"+str(Final_score)
		if not "." in all_ExAC:
			out.write(vcfline.rstrip()+"\t"+found_var_nt_pos+"\t"+found_var_aa_pos+"\t"+found_maineffect+"\t"+found_genename+"\t"+found_disease_desc+"\t"+found_gene_full_name+"\t"+str(found_cadd)+"\t"+str(found_svm)+"\t="+found_db+"\t"+str(found_phlyz)+"\t"+str(Final_score)+"\t"+found_intervar+"\t"+found_impact+"\t"+found_feature_id+"\t"+found_ann+"\t"+found_dbsnp+"\t"+found_go_biol_proc+"\t"+str(max_ExAC_AF)+"\t"+str(found_rvis)+"\t"+str(found_gdi))
			if warn == 1:
				out.write("\tMultiple Exac/CADD/SVM values\n")
			else:
				out.write("\t.\n")
		else:
			out.write(vcfline.rstrip()+"\t"+found_var_nt_pos+"\t"+found_var_aa_pos+"\t"+found_maineffect+"\t"+found_genename+"\t"+found_disease_desc+"\t"+found_gene_full_name+"\t"+str(found_cadd)+"\t"+str(found_svm)+"\t="+found_db+"\t"+str(found_phlyz)+"\t"+str(Final_score)+"\t"+found_intervar+"\t"+found_impact+"\t"+found_feature_id+"\t"+found_ann+"\t"+found_dbsnp+"\t"+found_go_biol_proc+"\t"+str(all_ExAC)+"\t"+str(found_rvis)+"\t"+str(found_gdi))
			if warn == 1:
                                out.write("\tMultiple Exac/CADD/SVM values\n")
                        else:
                                out.write("\t.\n")

vcf.close()
out.close()

