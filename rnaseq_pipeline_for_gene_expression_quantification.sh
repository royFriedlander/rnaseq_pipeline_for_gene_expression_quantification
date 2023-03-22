#! /bin/bash

# rsem-script-5.0.sh

Readme()
{
cat << EOF


Introduction:

rnaseq_pipeline_for_gene_expression_quantification.sh is a bash script that automates the bioinformatic processes required for estimating gene and isoform expression levels, given RNA-Seq data obtained from organisms with a reference genome and annotation.
Its purpose is to simplify things as much as possible.


Software Required:

1. RSEM    (https://github.com/deweylab/RSEM)
2. STAR    (https://github.com/alexdobin/STAR)
3. fastp   (https://github.com/OpenGene/fastp)
4. MultiQC (https://github.com/ewels/MultiQC)

Hardware Required:

A Linux server.


Usage Instructions:
(if you feel lost in the Linux command line environment, please check out this friendly tutorial: https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)

First,

    prepare the following dirs:
    1. Genome data dir
       contains:
       - a reference genome FASTA file  (file extension: *.fa or *.fasta)
       - an annotation GTF file         (file extension: *.gtf)
    2. RNA-Seq data dir
       contains:
       - RNA-Seq FASTQ files  (file extension: *.fastq or *.fastq.gz
       please notice:
       for paired-end sequencing data, filenames should be marked with 'R1' or 'R2' in one of the following formats:
       1. sample_name_R#.fastq ('R?' is between '_' and '.')
       2. sample_name.R#.fastq ('R?' is between '.' and '.')
       as '#' can be '1' or '2'
       
    for example:

    human
    └── LSLNGS2015
        ├── GENOME_data
        │   ├── Homo_sapiens.GRCh38.86.gtf
        │   ├── Homo_sapiens.GRCh38.dna.primary_assembly.fa
        └── RNASEQ_data
	    ├── GM12878.rep1.R1.fastq.gz
	    ├── GM12878.rep1.R2.fastq.gz
	    ├── GM12878.rep2.R1.fastq.gz
	    ├── GM12878.rep2.R2.fastq.gz
	    ├── K562.rep1.R1.fastq.gz
	    ├── K562.rep1.R2.fastq.gz
	    ├── K562.rep2.R1.fastq.gz
	    ├── K562.rep2.R2.fastq.gz


    in case your samples are split into several fastq files per sample
    (sequencing across multiple lanes might generate a separate file for each lane),
    you can use the 'cat' command to concatenate them together, as specified below:

    cat [filename 1] [filename 2] ... [filename N] > [new filename]

    * this will create a new file containing the concatenated text from all files given before the '>' character.

    for example, let's look at the forward read (Read 1) of sample 'K562.rep1'.
    in case it's split into the following 3 fastq files:
    K562.rep1.R1_L001.fastq.gz
    K562.rep1.R1_L002.fastq.gz
    K562.rep1.R1_L003.fastq.gz
    we can use the 'cat' command as followed:

    cat K562.rep1.R1_L001.fastq.gz K562.rep1.R1_L002.fastq.gz K562.rep1.R1_L003.fastq.gz > K562.rep1.R1.fastq.gz

    * this will create the file: 'K562.rep1.R1.fastq.gz',
      containing the sequencing data from all 3 files by the order given in the command,
      and ready to be analyzed in the pipeline.

Second,

    run the execution command specified below.


Execution Command:

nohup bash [path to rsem-script-5.0.sh] [path to genome data dir] [path to RNA-Seq data dir] [no. of threads] &

* 'nohup'  means: the process won't be killed in case the connection to the server fails.
* 'bash'   means: run the script specified next.
* '&'      means: run the script in the background (allows you to keep working with the command line while the pipeline is running).

for example:

nohup \
bash /home/roy/rnaseq-pipeline/rsem-script-5.0.sh \
/home/roy/data/human/LSLNGS2015/GENOME_data \
/home/roy/data/human/LSLNGS2015/RNASEQ_data \
10 \
&

('\' allows newlines within a command)

please notice:
the script needs 3 arguments to run properly:
1. a path to the reference genome data dir 
2. a path to the RNA-Seq data dir
3. no. of threads for multithread processing (if you're not sure - enter '10')



Outputs and results:

for every execution of the pipeline's script, a new directory will be created in your RNA-Seq data dir, under 'out/'.
the directory name will be determined by the time of execution (format: yymmdd-HHMM).
the new directory will contain:
- log.txt:       log file describing the pipeline's processes and activity
- output-data:   directory with all data created by the pipeline
- results:       copy of the main results (including rsem.genes.results, rsem.isoforms.results, and star's ReadsPerGene.out)




EOF
}


# arguments first check:
if [ $# -ne 3 ]
then
	echo >> $log_file
	echo "$# arguments received, but the pipeline's script needs 3 to operate properly. printing readme file and exiting..." >> $log_file
	Readme
	exit 0
fi


# read arguments:
genome_data_path=$1
RNASeq_data_path=$2
n_threads=$3


# getting reference genome fasta file name:
reference_file_name=$(find $genome_data_path -maxdepth 1 -name *.fa -printf "%f\n")
length=${#reference_file_name}
if [ $length -eq 0 ]
then
	reference_file_name=$(find $genome_data_path -maxdepth 1 -name *.fasta -printf "%f\n")
	length=${#reference_file_name}
	if [ $length -eq 0 ]
	then
		echo "couldn't find fasta file in $genome_data_path. exiting..." >> $log_file
		exit 0
	fi
fi
reference_name=${reference_file_name:0:7}


# getting reference annotation gtf/gff file name:
annotation_file_name=$(find $genome_data_path -maxdepth 1 -name *.gtf -printf "%f\n")
annotation_type='gtf'
length=${#annotation_file_name}
if [ $length -eq 0 ]
then
	annotation_file_name=$(find $genome_data_path -maxdepth 1 -name *.gff -printf "%f\n")
	annotation_type='gff3'
	length=${#annotation_file_name}
	if [ $length -eq 0 ]
	then
		echo "couldn't find gtf/gff file in $genome_data_path. exiting..." >> $log_file
		exit 0
	fi
	echo >> $log_file
	echo "-------------------------------------------------------------------------------------------" >> $log_file
	echo "please note:" >> $log_file
	echo >> $log_file
	echo "couldn't find a standard GTF annotation file. the pipeline will try to convert from gff to gtf, but errors may occur..." >> $log_file
	echo "if RSEM fails because the BAM file declares more reference sequences than RSEM knows (specified in the rsem.log file)," >> $log_file
	echo "please try to use a standard GTF annotation file. check out this tool: https://github.com/gpertea/gffread" >> $log_file
	echo "-------------------------------------------------------------------------------------------" >> $log_file
	echo >> $log_file
fi


# getting RNA-Seq fastq file names:
fastq_file_names=( $(ls $RNASeq_data_path | grep .fastq) )
length=${#fastq_file_names[@]}
if [ $length -eq 0 ]
then
	echo "couldn't find fastq files in $RNASeq_data_path. exiting..." >> $log_file
	exit 0
fi


# n_threads input checks:
if ((n_threads < 1 || n_threads > 1000))
then
	echo "no. of threads received: $n_threads. allowed values: between '1' and '1000'. if you're not sure, please enter '10'. exiting..." >> $log_file
	exit 0
fi


# extracting paired/single end option, and sample names:
is_paired_end=0
sample_names=()
for (( i=0; i<$length; i++ ));
do
	name=${fastq_file_names[i]}
	if [[ "$name" == *".R2."* ]] || [[ "$name" == *"_R2."* ]]
	then
		is_paired_end=1
	else
		# remove suffix and save sample name:
		name=${name%".fastq"*}
		name=${name%"_R1"}
		name=${name%".R1"}
		sample_names+=($name)
	fi
done

n_samples=${#sample_names[@]}


# paired/single end conclusion check:
if [ $is_paired_end -eq 1 ]
then
	let double_samples=2*${#sample_names[@]}
	if [ ${#fastq_file_names[@]} -ne $double_samples ]
	then
		echo "no. of fastq filenames do not match no. of sample names for paired-end seq data (1:2 ratio expected). exiting..." >> $log_file
		exit 0
	fi
else
	if [ ${#fastq_file_names[@]} -ne ${#sample_names[@]} ]
	then
		echo "no. of fastq filenames do not match no. of sample names for single-end seq data (1:1 ratio expected). exiting..." >> $log_file
		exit 0
	fi
fi



# prepare output destinations

timestamp=$(date +%y%m%d-%H%M)
mkdir $RNASeq_data_path/out
mkdir $RNASeq_data_path/out/$timestamp
mkdir $RNASeq_data_path/out/$timestamp/output-data
mkdir $RNASeq_data_path/out/$timestamp/results
log_file=$RNASeq_data_path/out/$timestamp/log.txt


# configuration summary

echo rsem-script-5.0.sh >> $log_file
date +%d-%m-%y >> $log_file
date +%T >> $log_file
echo >> $log_file

cat << EOF >> $log_file

Configurations:

genome data path:
$genome_data_path

RNASeq data path:
$RNASeq_data_path

reference file name:
$reference_file_name

reference name:
$reference_name

reference annotation file name:
$annotation_file_name

annotation type:
$annotation_type

EOF

echo "fastq file names:" >> $log_file
printf '%s\n' ${fastq_file_names[@]} >> $log_file
echo >> $log_file

if [ $is_paired_end -eq 1 ]
then
	end="paired-end"
else
	end="single-read"
fi

echo "sample names ($n_samples, $end):" >> $log_file
printf '%s\n' ${sample_names[@]} >> $log_file
echo >> $log_file

echo "aligner:" >> $log_file
echo "STAR" >> $log_file
echo >> $log_file

echo "no. of threads: $n_threads" >> $log_file
echo >> $log_file



# Create mapping indices:

if [ ! -d "$genome_data_path/star" ] 
then
	# mapping indices for STAR:
	: '
	usage:
	STAR --runMode genomeGenerate \
	     --genomeDir path_to_genomedir \
	     --genomeFastaFiles reference_fasta_file(s)
	'
	echo "Creating mapping indices for STAR in $genome_data_path/star..." >> $log_file
	date +"%T" >> $log_file
	mkdir $genome_data_path/star
	STAR --runThreadN $n_threads \
	     --runMode genomeGenerate \
	     --genomeDir $genome_data_path/star \
	     --genomeFastaFiles $genome_data_path/$reference_file_name \
	     --sjdbGTFfile $genome_data_path/$annotation_file_name
	echo "Creating mapping indices for STAR is done." >> $log_file
	echo >> $log_file
else
	echo "Mapped indices for STAR were found in $genome_data_path/star" >> $log_file
        echo >> $log_file
fi


if [ ! -d "$genome_data_path/rsem" ] 
then
	# mapping indices for RSEM:
	: '
	usage:
	rsem-prepare-reference [options] reference_fasta_file(s) reference_name
	'
	echo "Creating mapping indices for RSEM in $genome_data_path/rsem..." >> $log_file
	date +"%T" >> $log_file
	mkdir $genome_data_path/rsem
	rsem-prepare-reference \
		--$annotation_type \
		$genome_data_path/$annotation_file_name \
		$genome_data_path/$reference_file_name \
		$genome_data_path/rsem/$reference_name
	echo "Creating mapping indices for RSEM is done." >> $log_file
	echo >> $log_file
else
	echo "Mapped indices for RSEM were found in $genome_data_path/rsem" >> $log_file
        echo >> $log_file
fi



# Trim with fastp:

echo "Trimming and producing quality reports with fastp..." >> $log_file
mkdir $RNASeq_data_path/out/$timestamp/output-data/trimmed
mkdir $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports

if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		r1_fastq_file=${fastq_file_names[2*i]}
		r2_fastq_file=${fastq_file_names[2*i+1]}
		echo >> $log_file
		echo "trimming sample $sample_name" >> $log_file
		date +"%T" >> $log_file
		echo "input files: $r1_fastq_file and $r2_fastq_file" >> $log_file
		fastp -i $RNASeq_data_path/$r1_fastq_file \
		      -I $RNASeq_data_path/$r2_fastq_file \
		      -o $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.R1.P.fq.gz \
		      -O $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.R2.P.fq.gz \
		      -w 16

		mkdir $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/$sample_name
		mv fastp.html fastp.json $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/$sample_name
		echo "done." >> $log_file
                date +"%T" >> $log_file

	done
else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		fastq_file=${fastq_file_names[i]}
		echo >> $log_file
		echo "trimming sample $sample_name" >> $log_file
		date +"%T" >> $log_file
		echo "input file: $fastq_file" >> $log_file
		fastp -i $RNASeq_data_path/$fastq_file \
		      -o $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.fq.gz \
		      -w 16

		mkdir $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/$sample_name
		mv fastp.html fastp.json $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/$sample_name
		echo "done." >> $log_file
                date +"%T" >> $log_file

	done
fi
echo >> $log_file
echo "Trimming with fastp is done." >> $log_file
echo >> $log_file
echo "trimmed files are here:   $RNASeq_data_path/out/$timestamp/output-data/trimmed" >> $log_file
echo "quality reports are here: $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports" >> $log_file
echo >> $log_file


# Alignment:

# Mapping with STAR:

echo "Alignment: mapping with STAR..." >> $log_file
if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo >> $log_file
		echo "mapping sample $sample_name" >> $log_file
		date +"%T" >> $log_file

		mkdir $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name
		STAR --genomeDir $genome_data_path/star \
		     --readFilesCommand zcat \
		     --readFilesIn $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.R1.P.fq.gz \
		                   $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.R2.P.fq.gz \
		     --outSAMtype BAM SortedByCoordinate \
		     --limitBAMsortRAM 16000000000 \
		     --outSAMunmapped Within \
		     --twopassMode Basic \
		     --outFilterMultimapNmax 1 \
		     --quantMode TranscriptomeSAM GeneCounts\
		     --runThreadN $n_threads \
		     --outFileNamePrefix "$RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/"
		mkdir $RNASeq_data_path/out/$timestamp/results/star_$sample_name
                cp $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/ReadsPerGene.out.tab \
		   $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/Log.out \
                   $RNASeq_data_path/out/$timestamp/results/star_$sample_name
		echo "done." >> $log_file
                date +"%T" >> $log_file

	done
else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo >> $log_file
		echo "mapping sample $sample_name" >> $log_file
		date +"%T" >> $log_file

		mkdir $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name
		STAR --genomeDir $genome_data_path/star \
		     --readFilesCommand zcat \
		     --readFilesIn $RNASeq_data_path/out/$timestamp/output-data/trimmed/$sample_name.fq.gz \
		     --outSAMtype BAM SortedByCoordinate \
		     --limitBAMsortRAM 16000000000 \
		     --outSAMunmapped Within \
		     --twopassMode Basic \
		     --outFilterMultimapNmax 1 \
		     --quantMode TranscriptomeSAM GeneCounts\
		     --runThreadN $n_threads \
		     --outFileNamePrefix "$RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/"
		mkdir $RNASeq_data_path/out/$timestamp/results/star_$sample_name
                cp $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/ReadsPerGene.out.tab \
		   $RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/Log.out \
                   $RNASeq_data_path/out/$timestamp/results/star_$sample_name
		echo "done." >> $log_file
                date +"%T" >> $log_file

	done
fi
echo >> $log_file
echo "Mapping with STAR is done." >> $log_file
echo >> $log_file


# Quantification with RSEM:

: '
usage:
rsem-calculate-expression [optioins] upstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
rsem-calculate-expression [options] --sam/--bam [--paired-end] input reference_name sample_name
'

echo "Quantifying gene expression with RSEM..." >> $log_file
echo >> $log_file
if [ $is_paired_end -eq 1 ]
then # for paired-end seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo "quantifying sample $sample_name" >> $log_file
		date +"%T" >> $log_file

		mkdir $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name
		rsem-calculate-expression \
			--bam \
			--no-bam-output \
			-p $n_threads \
			--paired-end \
			--forward-prob 0 \
			$RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/Aligned.toTranscriptome.out.bam \
			$genome_data_path/rsem/$reference_name \
			$RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem >& \
			$RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.log

		echo "done." >> $log_file
		date +"%T" >> $log_file
		echo "the resaults are here: $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name" >> $log_file
		echo >> $log_file
		mkdir $RNASeq_data_path/out/$timestamp/results/rsem_$sample_name
		cp $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.genes.results \
		   $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.isoforms.results \
		   $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.log \
		   $RNASeq_data_path/out/$timestamp/results/rsem_$sample_name

	done

else # for single-read seq
	for (( i=0; i<$n_samples; i++ ));
	do
		sample_name=${sample_names[i]}
		echo "quantifying sample $sample_name" >> $log_file
		date +"%T" >> $log_file

		mkdir $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name

		rsem-calculate-expression \
			--bam \
			--no-bam-output \
			-p $n_threads \
			--forward-prob 0 \
			$RNASeq_data_path/out/$timestamp/output-data/star_$sample_name/Aligned.toTranscriptome.out.bam \
			$genome_data_path/rsem/$reference_name \
			$RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem >& \
			$RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.log

		echo "done." >> $log_file
		date +"%T" >> $log_file
		echo "the resaults are here: $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name" >> $log_file
		echo >> $log_file
		mkdir $RNASeq_data_path/out/$timestamp/results/rsem_$sample_name
		cp $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.genes.results \
		   $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.isoforms.results \
		   $RNASeq_data_path/out/$timestamp/output-data/rsem_$sample_name/rsem.log \
		   $RNASeq_data_path/out/$timestamp/results/rsem_$sample_name
	done
fi
echo "Quantification with RSEM is done." >> $log_file
echo >> $log_file

echo "Summarizing qc reports with multiqc..." >> $log_file
date +"%T" >> $log_file

multiqc -o $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/ $RNASeq_data_path/out/$timestamp
mkdir $RNASeq_data_path/out/$timestamp/results/qc-reports
cp $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/multiqc_report.html \
   $RNASeq_data_path/out/$timestamp/results/qc-reports
cp -r $RNASeq_data_path/out/$timestamp/output-data/trimmed/quality-reports/multiqc_data \
      $RNASeq_data_path/out/$timestamp/results/qc-reports

echo "done." >> $log_file
date +"%T" >> $log_file

echo >> $log_file
echo >> $log_file
echo "you'll find the pipelin's log file and a copy of the results here: $RNASeq_data_path/out/$timestamp" >> $log_file
echo >> $log_file
echo >> $log_file

