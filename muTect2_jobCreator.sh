#!/bin/bash

echo "Sample ID: "
read ID

if [ ${#ID} -gt 11 ]
then
	shortID=$(echo $ID | cut -c 1-11)
fi

#echo  -e ">Select reference: [Type hg19/hg38] (default hg19): \c"
#read genomeref
#if [ -z "$genomeref" ]
#then
#        genomeref=hg19
#        echo -e $genomeref
#fi

echo -e "number of Cores (default 12): \c"
read CPU
if [ -z "$CPU" ]
then
	CPU=12
	echo -e $CPU
fi

echo -e "Amount of Memory (default 120Gb): \c"
read RAM
if [ -z "$RAM" ]
then
	RAM=120
	echo -e $RAM
fi

echo -e "Name of the queue (default "parallel"): \c"
read QUEUE
if [ -z "$QUEUE" ]
then
	QUEUE=parallel
	echo $QUEUE
fi

echo -e "#SBATCH -N 1 -n "$CPU" --mem "$RAM"GB" >> $ID"_muTect2_Job"
echo -e "#SBATCH -J "$shortID			>> $ID"_muTect2_Job"
echo -e "#SBATCH -A "$ACCOUNT			>> $ID"_muTect2_Job"
echo -e "#SBATCH --time 100:00:00"		>> $ID"_muTect2_Job"
echo -e "#SBATCH --error %j.%x.err"		>> $ID"_muTect2_Job"
echo -e "#SBATCH --output %j.%x.err"		>> $ID"_muTect2_Job"
echo -e "#SBATCH --partition bdw_usr_prod\n"	>> $ID"_muTect2_Job"

echo -e "umask 0002\n" >> $ID"_muTect2_Job" 

echo "Control Sample Fastq Directory: "
read ctrl_fastq_di
ctrl_fastq_dir=$ctrl_fastq_di"/"
if [ ! -d $ctrl_fastq_dir ]
then
        echo $ctrl_fastq_dir "does not exist!"
        exit 1
fi

echo -e "ctrl_fastq_dir="$ctrl_fastq_dir >> $ID"_muTect2_Job"


echo "How many fastq couples for the control sample? "
read ctrl_fq_couple

echo -e "ctrl_fq_couple="$ctrl_fq_couple >> $ID"_muTect2_Job"

declare -a ctrl_names=( $(seq 1 $ctrl_fq_couple) )
for i in ${ctrl_names[@]}; do
        echo "First-pair fastqfile (e.g. _1.fq.gz) of couple "$i
        read ctrl_fastqname1
	echo -e "ctrl_fastq_name_original_1["$i"]="$ctrl_fastqname1 >> $ID"_muTect2_Job"
        echo "Second-pair fastqfile (e.g. _2.fq.gz) of couple "$i
        read ctrl_fastqname2
	echo -e "ctrl_fastq_name_original_2["$i"]="$ctrl_fastqname2 >> $ID"_muTect2_Job"
	echo "Lane ID for fastq couple "$i
	read ctrl_laneid
	echo -e "ctrl_lane_id["$i"]="$ctrl_laneid >> $ID"_muTect2_Job"
	ctrl_fastqlane1[$i]=$ctrl_fastqname1
        ctrl_fastqlane2[$i]=$ctrl_fastqname2
        if [ ! -f ${ctrl_fastq_dir}${ctrl_fastqlane1[$i]} ]
        then
                echo ${ctrl_fastq_dir}${ctrl_fastqlane1[$i]} " does not exist!"
                exit 1
        fi
        if [ ! -f ${ctrl_fastq_dir}${ctrl_fastqlane2[$i]} ]
        then
                echo ${ctrl_fastq_dir}${ctrl_fastqlane2[$i]} " does not exist!"
                exit 1
        fi
done

echo "Tumor Sample Fastq Directory: "
read tum_fastq_di
tum_fastq_dir=$tum_fastq_di"/"
if [ ! -d $tum_fastq_dir ]
then
        echo $tum_fastq_dir "does not exist!"
        exit 1
fi

echo -e "tum_fastq_dir="$tum_fastq_dir >> $ID"_muTect2_Job"

echo "How many fastq couples for the tumor sample? "
read tum_fq_couple

echo -e "tum_fq_couple="$tum_fq_couple >> $ID"_muTect2_Job"

declare -a tum_names=( $(seq 1 $tum_fq_couple) )
for i in ${tum_names[@]}; do
        echo "First-pair fastqfile (e.g. _1.fq.gz) of couple "$i
        read tum_fastqname1
        echo -e "tum_fastq_name_original_1["$i"]="$tum_fastqname1 >> $ID"_muTect2_Job"
        echo "Second-pair fastqfile (e.g. _2.fq.gz) of couple "$i
        read tum_fastqname2
        echo -e "tum_fastq_name_original_2["$i"]="$tum_fastqname2 >> $ID"_muTect2_Job"
        echo "Lane ID for fastq couple "$i
        read tum_laneid
        echo -e "tum_lane_id["$i"]="$tum_laneid >> $ID"_muTect2_Job"
        tum_fastqlane1[$i]=$tum_fastqname1
        tum_fastqlane2[$i]=$tum_fastqname2
        if [ ! -f ${tum_fastq_dir}${tum_fastqlane1[$i]} ]
        then
                echo ${tum_fastq_dir}${tum_fastqlane1[$i]} " does not exist!"
                exit 1
        fi
        if [ ! -f ${tum_fastq_dir}${tum_fastqlane2[$i]} ]
        then
                echo ${tum_fastq_dir}${tum_fastqlane2[$i]} " does not exist!"
                exit 1
        fi
done

echo ">Target regions path: "
read targetRegChipDir
while [ ! -d $targetRegChipDir ]; do
	echo $targetRegChipDir " does not exist!"
	echo ">Target regions path: "
	read targetRegChipDir
done

echo "Output directory: "
read outdi
outdir=$outdi"/"
if [ ! -d $outdir ]
then
        echo $outdir "does not exist!"
        exit 1
fi

echo "Temporary directory: "
read tmpdi
tmpdir=$tmpdi"/"
if [ ! -d $tmpdir ]
then
        echo $tmpdir "does not exist!"
        exit 1
fi

echo "Do you want FastQC data pre-processing? [Type y/n]"
read fqc

echo "Do you have Illumina (Q PHRED+64) or Sanger fastq (Q PHRED+33)? [Type i/s]"
read qscore

echo ">Do you want splice sites analysis? [Type y/n]"
read spid
echo -e "spidex="$spid >> $ID"_muTect2_Job"	


echo "Do you want to calculate global coverage? [Type y/n]"
read coverBedTools

echo "Do you want to calculate HsMetrics? [Type y/n]"
read metrics

if [ $fqc == "y" ]; then
	if [ $fq_couple -gt 1 ]; then
		echo -e "\n  concatenating fastq files..."
		CWD=$(pwd)
		cd  ${ctrl_fastq_dir}
		cat $( echo ${ctrl_fastqlane1[@]} ) > $ID"_ctrl_all_fq_reads_1.fq.gz"
		cat $( echo ${ctrl_fastqlane2[@]} ) > $ID"_ctrl_all_fq_reads_2.fq.gz"
		cd  ${tum_fastq_dir}
		cat $( echo ${tum_fastqlane1[@]} ) > $ID"_"$ID"_tum_all_fq_reads_1.fq.gz"
		cat $( echo ${tum_fastqlane2[@]} ) > $ID"_"$ID"_tum_all_fq_reads_2.fq.gz"
		cd ${CWD}
	else
		CWD=$(pwd)
		cd  ${ctrl_fastq_dir}
		ln -s ${ctrl_fastq_dir}${ctrl_fastqname1} $outdir$ID"_all_fq_reads_1.fq.gz"
		ln -s ${ctrl_fastq_dir}${ctrl_fastqname2} $outdir$ID"_all_fq_reads_2.fq.gz"
		cd  ${tum_fastq_dir}
		ln -s ${tum_fastq_dir}${tum_fastqname1} $outdir$ID"_tum_all_fq_reads_1.fq.gz"
		ln -s ${tum_fastq_dir}${tum_fastqname2} $outdir$ID"_tum_all_fq_reads_2.fq.gz"
		cd ${CWD}
	fi
fi

echo -e "outdir="$outdir			>> $ID"_muTect2_Job"
echo -e "ID="$ID				>> $ID"_muTect2_Job"
echo -e "tmpdir="$tmpdir			>> $ID"_muTect2_Job"
echo -e "fqc="$fqc				>> $ID"_muTect2_Job"
echo -e "qscore="$qscore			>> $ID"_muTect2_Job"
echo -e "coverBedTools="$coverBedTools		>> $ID"_muTect2_Job"
echo -e "metrics="$metrics			>> $ID"_muTect2_Job"
echo -e "targetRegChipDir="$targetRegChipDir	>> $ID"_muTect2_Job"

cat "/pico/work/IscrC_FoRWArDS_1/NGS_tools/Pipeline/Tartagenia/v3.1/muTect2_pipeline_part.txt" >> $ID"_muTect2_Job"
mv -f $ID"_muTect2_Job" $outdir

##Creo il Job per l'annotazione del file muTect2
echo -e "#!/bin/bash"				 > $ID"_muTect2_annot_Job"
echo -e "#SBATCH -N 1 -n 6 --mem 120GB"		>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH -J "$shortID"_ann"		>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH -A "$ACCOUNT			>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH --time 12:00:00"		>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH --error %j.%x.err"		>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH --output %j.%x.err"		>> $ID"_muTect2_annot_Job"
echo -e "#SBATCH --partition bdw_usr_prod"	>> $ID"_muTect2_annot_Job"
echo -e "\numask 0002\n"			>> $ID"_muTect2_annot_Job"
echo "## definizione delle variabili di input" 	>> $ID"_muTect2_annot_Job"
echo -e "ID="$ID				>> $ID"_muTect2_annot_Job"
echo -e "outdir="$outdir			>> $ID"_muTect2_annot_Job"
echo -e "vcfinput="$outdir$ID"_muTect2_raw_snps-indels.vcf.gz" >> $ID"_muTect2_annot_Job"

cat "/pico/work/IscrC_FoRWArDS_1/NGS_tools/Pipeline/Tartagenia/v3.1/muTect2_custom_annotation_part.txt" >> $ID"_muTect2_annot_Job"
mv -f $ID"_muTect2_annot_Job" $outdir

##Creo il Job per l'annotazione del file ctrl
echo -e "#!/bin/bash"				 > $ID"_annot_Job"
echo -e "#SBATCH -N 1 -n 6 --mem 120GB"		>> $ID"_annot_Job"
echo -e "#SBATCH -J "$shortID"_ann"		>> $ID"_annot_Job"
echo -e "#SBATCH -A "$ACCOUNT			>> $ID"_annot_Job"
echo -e "#SBATCH --time 12:00:00"		>> $ID"_annot_Job"
echo -e "#SBATCH --error %j.%x.err"		>> $ID"_annot_Job"
echo -e "#SBATCH --output %j.%x.err"		>> $ID"_annot_Job"
echo -e "#SBATCH --partition bdw_usr_prod"	>> $ID"_annot_Job"
echo -e "\numask 0002\n"			>> $ID"_annot_Job"
echo "## definizione delle variabili di input"	>> $ID"_annot_Job"
echo -e "ID="$ID				>> $ID"_annot_Job"
echo -e "outdir="$outdir			>> $ID"_annot_Job"
echo -e "vcfinput="$outdir$ID"_raw_snps-indels_HapCall_genotype_filtered.g.vcf.gz" >> $ID"_annot_Job"

cat "/pico/work/IscrC_FoRWArDS_1/NGS_tools/Pipeline/Tartagenia/v3.1/muTect2_ctrl_tum_custom_annotation_part.txt" >> $ID"_annot_Job"
mv -f $ID"_annot_Job" $outdir

##Creo il Job per l'annotazione del file tum
echo -e "#!/bin/bash" > $ID"_tum_annot_Job"
echo -e "#SBATCH -N 1 -n 6 --mem 120GB"		>> $ID"_tum_annot_Job"
echo -e "#SBATCH -J "$shortID"_t_ann"		>> $ID"_tum_annot_Job"
echo -e "#SBATCH -A "$ACCOUNT			>> $ID"_tum_annot_Job"
echo -e "#SBATCH --time 12:00:00"		>> $ID"_tum_annot_Job"
echo -e "#SBATCH --error %j.%x.err"		>> $ID"_tum_annot_Job"
echo -e "#SBATCH --output %j.%x.err"		>> $ID"_tum_annot_Job"
echo -e "#SBATCH --partition bdw_usr_prod"	>> $ID"_tum_annot_Job"
echo -e "\numask 0002\n" >> $ID"_tum_annot_Job"
echo -e "## definizione delle variabili di input" >> $ID"_tum_annot_Job"
echo -e "ID="$ID >> $ID"_tum_annot_Job"
echo -e "outdir="$outdir >> $ID"_tum_annot_Job"
echo -e "vcfinput="$outdir$ID"_tum_raw_snps-indels_HapCall_genotype_filtered.g.vcf.gz" >> $ID"_tum_annot_Job"

cat "/pico/work/IscrC_FoRWArDS_1/NGS_tools/Pipeline/Tartagenia/v3.1/muTect2_ctrl_tum_custom_annotation_part.txt" >> $ID"_tum_annot_Job"
mv -f $ID"_tum_annot_Job" $outdir

echo -e "\n\tFor the somatic variants annotation: sbatch\t"$outdir$ID"_muTect2_annot_Job\n\n\tFor the control variants annotation: sbatch\t"$outdir$ID"_annot_Job\n\n\tFor the tumor variants annotation: sbatch\t"$outdir$ID"_tum_annot_Job\n" > $outdir$ID"_muTect2_custom_annotation.INFO" 
echo -e "\n\tsbatch this file: "$outdir$ID"_muTect2_Job\n\n\tinstructions for the annotation of the final vcf can be found in "$outdir$ID"_muTect2_custom_annotation.INFO\n"
