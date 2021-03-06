#!/bin/bash

source /broad/software/scripts/useuse

use BWA
use Java-1.8
use Picard-Tools
use Anaconda
use Samtools
use BEDTools

#sampleid=$1
#temp_dir=$2
#hg19_ref=$3
#ref_path=$4
#runid=$5
#gatk=$6
#dict=$7
#tool=$8

cores=8
tool="/seq/plasmodium/tools/malaria_variant_calling"
app="/cil/shed/apps/external"
gatk="/humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar"
run_dir="/seq/plasmodium/test/test2"
runid="testrun1"
ref_path="/gsap/garage-protistvector/U19_Aim4/REF/PlasmoDB-28_Pfalciparum3D7_Genome.fasta"
dict="/gsap/garage-protistvector/U19_Aim4/REF/PlasmoDB-28_Pfalciparum3D7_Genome.dict"
temp_dir=${run_dir}/other_files
raw_fq=1
hg19_ref="/seq/plasmodium/data/refs/hg19_db/hg19_db"
metafile=${run_dir}/bam_meta.tsv

# Parse metafile
myline=$(sed -n "${SGE_TASK_ID}"p ${metafile})
read -ra INFO <<<"$myline"

sampleid=${INFO[1]} ## Check against metafile
bampath=${INFO[0]} ## check against metafile

### hg19 alignment - host removal
# map to BWA hg19 db index
echo "Mapping of sample FASTQs to BWA hg19 db started"
use BWA
bwa mem -t ${cores} ${hg19_ref} ${temp_dir}/fastq/${sampleid}.1.fq ${temp_dir}/fastq/${sampleid}.2.fq | samtools view -bS -> ${temp_dir}/${sampleid}.mapAndUnmapped.hg19.bam
echo "Mapping of sample FASTQs to hg19 db finished"
# filter out reads unmapped to hg19
echo "Filtering of unmapped reads started"
samtools view -b -f 12 -F 256 ${temp_dir}/${sampleid}.mapAndUnmapped.hg19.bam > ${temp_dir}/${sampleid}.unmapped.bam
echo "Filtering of unmapped reads finished"
# sort BAM file by read name (-n) to have reads next to each other [required by bedtools]
echo "Unmapped BAM sorting started"
samtools sort -n ${temp_dir}/${sampleid}.unmapped.bam -o ${temp_dir}/${sampleid}.unmapped.sorted.bam
echo "Unmapped BAM sorting finished"
# BAM to FASTQ files
use BEDTools
echo "BAM to FASTQ started"
bedtools bamtofastq -i ${temp_dir}/${sampleid}.unmapped.sorted.bam -fq ${temp_dir}/fastq/${sampleid}.hostRemoved.1.fq -fq2 ${temp_dir}/fastq/${sampleid}.hostRemoved.2.fq
echo "BAM to FASTQ finished"

### PFal Alignment
echo "PFal Alignment started"
bwa mem -t ${cores} \
-R "@RG\\tID:FLOWCELL_${sampleid}\\tSM:${sampleid}\\tPL:ILLUMINA\\tLB:LIB_${sampleid}" ${ref_path} \
${temp_dir}/fastq/${sampleid}.hostRemoved.1.fq ${temp_dir}/fastq/${sampleid}.hostRemoved.2.fq | samtools view -bS -> ${temp_dir}/${sampleid}.aligned.bam
echo "PFal Alignment finished"

### Sort BAM
java -Xmx8G -jar $PICARD SortSam I=${temp_dir}/${sampleid}.aligned.bam \
O=${temp_dir}/${sampleid}.sorted.bam \
SO=coordinate

### Mark Duplicates
java -Xmx8G -jar $PICARD MarkDuplicates I=${temp_dir}/${sampleid}.sorted.bam \
O=${temp_dir}/${sampleid}.marked_duplicates.bam \
M=${temp_dir}/${sampleid}.marked_duplicates.metrics

### Re-order BAM
java -Xmx8G -jar $PICARD ReorderSam I=${temp_dir}/${sampleid}.marked_duplicates.bam \
O=${temp_dir}/${sampleid}.reordered.bam R=${ref_path} SD=${dict}
samtools index ${temp_dir}/${sampleid}.reordered.bam

### Read Metrics / Picard Metrics / DepthOfCoverage
# Multiple Metrics
Picard Multiple Metrics
java -jar $PICARD CollectMultipleMetrics I=${temp_dir}/${sampleid}.reordered.bam \
O=${temp_dir}/${sampleid}_multiple_metrics \
R=${ref_path}

# Parse Multiple metrics and add some new ones
python ${tool}/QC_subProc.py --sample_path ${temp_dir}/${sampleid}.reordered.bam \
--sampleid ${sampleid} --reference ${ref_path} --path_to_picard $PICARD \
--path_to_gatk3 ${gatk} \
--path_to_dir ${temp_dir} --output_file ${run_dir}/${runid}_qc_summary.txt

### GATK Variant calling
# GATK RealignerTargetCreator
java -Xmx8G -jar ${gatk} -T RealignerTargetCreator -nct 1 -nt ${cores} \
-R ${ref_path} -I ${temp_dir}/${sampleid}.reordered.bam \
-o ${temp_dir}/${sampleid}.interval_list

# GATK IndelRealigner
java -Xmx4G -jar ${gatk} -T IndelRealigner -nct 1 -nt 1 \
-R ${ref_path} -I ${temp_dir}/${sampleid}.reordered.bam \
-targetIntervals ${temp_dir}/${sampleid}.interval_list \
-o ${temp_dir}/${sampleid}.indels_realigned.bam
samtools index ${temp_dir}/${sampleid}.indels_realigned.bam

# GATK BaseRecalibrator
java -Xmx4G -jar ${gatk} -T BaseRecalibrator -nct 8 -nt 1 \
-R ${ref_path} -I ${temp_dir}/${sampleid}.indels_realigned.bam \
-knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz \
-knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz \
-knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz \
-o ${temp_dir}/${sampleid}_recal_report.grp

# GATK Print Reads (BaseRecalibrator)
java -Xmx4G -jar ${gatk} -T PrintReads -nct 8 -nt 1 \
-R ${ref_path} -I ${temp_dir}/${sampleid}.indels_realigned.bam \
-BQSR ${temp_dir}/${sampleid}_recal_report.grp \
-o ${temp_dir}/${sampleid}.bqsr.bam
samtools index ${temp_dir}/${sampleid}.bqsr.bam

# GATK HaplotypeCaller
java -Xmx8G -jar ${gatk} -T HaplotypeCaller -nt 1 \
-R ${ref_path} --input_file ${temp_dir}/${sampleid}.bqsr.bam \
-ERC GVCF -ploidy 2 --interval_padding 100 -o ${run_dir}/gvcf/${sampleid}.g.vcf \
-variant_index_type LINEAR -variant_index_parameter 128000


### Cleanup
rm ${temp_dir}/fastq/${sampleid}.1.fq
rm ${temp_dir}/fastq/${sampleid}.2.fq
rm ${temp_dir}/${sampleid}.aligned.bam
rm ${temp_dir}/${sampleid}.sorted.bam
rm ${temp_dir}/${sampleid}.marked_duplicates.bam
rm ${temp_dir}/${sampleid}.unmapped.bam
rm ${temp_dir}/${sampleid}.marked_duplicates.metrics
rm ${temp_dir}/${sampleid}.reordered.bam
rm ${temp_dir}/${sampleid}.reordered.bam.bai
rm ${temp_dir}/${sampleid}.interval_list
rm ${temp_dir}/${sampleid}.indels_realigned.bam
rm ${temp_dir}/${sampleid}.indels_realigned.bam.bai
rm ${temp_dir}/${sampleid}.indels_realigned.bai
rm ${temp_dir}/${sampleid}_recal_report.grp
