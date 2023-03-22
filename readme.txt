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

nohup bash [path to rnaseq_pipeline_for_gene_expression_quantification.sh] [path to genome data dir] [path to RNA-Seq data dir] [no. of threads] &

* 'nohup'  means: the process won't be killed in case the connection to the server fails.
* 'bash'   means: run the script specified next.
* '&'      means: run the script in the background (allows you to keep working with the command line while the pipeline is running).

for example:

nohup \
bash /home/roy/rnaseq-pipelines/rnaseq_pipeline_for_gene_expression_quantification.sh \
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

