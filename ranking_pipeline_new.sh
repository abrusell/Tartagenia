#!/bin/bash -e

if [ $# -ne 3 ]
then
	echo -e "Usage: $0 ID_phenotype_list.txt ID outdir\n\tN.B. In the outdir the file ID_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd_ACMG_DDG2P_Mendel_intervar_db.vcf.gz must be present!"
	exit 1
fi

## definizione delle variabili di input
phenolist=$1
ID=$2
outdi=$3
outdir=$outdi"/"

cd $outdir

if [ -s $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd_spidex_ACMG_DDG2P_Mendel_intervar_db.vcf.gz" ] 
then
        spid=_spidex
elif [ ! -s $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd_ACMG_DDG2P_Mendel_intervar_db.vcf.gz" ] 
then
        echo -e "Check your input files. In the outdir, ID_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd(_spidex)_ACMG_DDG2P_Mendel_intervar_db.vcf.gz must be present."
        exit 1
else
        spid=
fi


# filtro per le CDS e i SS
python /pico/work/IscrC_FoRWArDS_1/NGS_tools/vcfSifter.py $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db.vcf.gz" 2> $outdir$ID"_ranking_log"

# estraggo la lista dei geni nel vcf
python /pico/work/IscrC_FoRWArDS_1/NGS_tools/vcf2genelist.py $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.vcf" 2>> $outdir$ID"_ranking_log"

# ordino e elimino duplicati
sort -u $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.genelist.txt" > $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered_genelist_sorted.txt"
rm $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.genelist.txt"

# qui uso la lista con i max 4 termini fenotipici che avro' precedentemente creato
perl /pico/home/userexternal/aciolfi0/phenolyzer/disease_annotation.pl -f -p --gene $ID"_raw_snps-indels_HapCall_genotype_filtered_genelist_sorted.txt" -ph $phenolist -logistic -out "Pheno/"$ID -addon DB_DISGENET_GENE_DISEASE_SCORE,DB_GAD_GENE_DISEASE_SCORE -addon_weight 0.25 2>> $outdir$ID"_ranking_log"

# aggiungo il phenolyzer score al vcf
python /pico/work/IscrC_FoRWArDS_1/NGS_tools/phenolyzer_score_annotator.py $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.vcf" $outdir"Pheno/"$ID".final_gene_list" 2>> $outdir$ID"_ranking_log"

# applico il nostro schema di scoring alle varianti
python /pico/work/IscrC_FoRWArDS_1/NGS_tools/prioritizer_maxPopAF_gnomAD_MCAP.py $phenolist $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.phenolyzer.vcf" 2>> $outdir$ID"_ranking_log"

# creo dei file rinominati in maniera piu' sintetica da consegnare ai biologi che hanno windows...
mkdir $outdir"RENAMED"
cp $outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.phenolyzer.scored.tsv" $outdir"RENAMED/"
# rinomino l'output per il singleton
rename "_raw_snps-indels_HapCall_genotype_filtered.g.SnpEff_UCSCAnnot_dbNSFP3.dbNSFP3_gene_cadd"$spid"_ACMG_DDG2P_Mendel_intervar_db_filtered.phenolyzer.scored" "_raw_snps-indels_HapCall_genotype_filtered.phenolyzer.scored" $outdir"RENAMED/"*"tsv"
