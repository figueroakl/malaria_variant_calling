#!/bin/bash

source /broad/software/scripts/useuse

use Java-1.8
use Picard-Tools
use Anaconda

tool="/seq/plasmodium/tools/malaria_variant_calling"
run_dir="/seq/plasmodium/test/test2"
temp_dir=${run_dir}/other_files
raw_fq=1
metafile=${run_dir}/bam_meta.tsv
app="/cil/shed/apps/external"
fastqc="${app}/FastQC/fastqc"


# Parse metafile
myline=$(sed -n "${SGE_TASK_ID}"p ${metafile})
read -ra INFO <<<"$myline"

sampleid=${INFO[1]} ## Check against metafile
bampath=${INFO[0]} ## check against metafile

path_to_fq="${temp_dir}/fastq/${sampleid}"

# Bam to fastq
java -Xmx12G -jar $PICARD SamToFastq INPUT=${bampath} \
FASTQ=${path_to_fq}.1.fq \
SECOND_END_FASTQ=${path_to_fq}.2.fq \
VALIDATION_STRINGENCY=LENIENT

# Run fastQC on obtained fq files (optional)
if [ ${raw_fq} -eq 1 ]; then
	${fastqc} ${path_to_fq}.*.fq
else
	echo "Skipping fastqc step as per User request"
fi
