#!/bin/bash

source /broad/software/scripts/useuse

SampleID=$1
BamPath=$2
ReferenceFile=$3
IntermedFilePath=$4
GvcfPath=$5

### SamToFastq
use Java-1.8
java -Xmx12G -jar /cil/shed/apps/external/picard/current/bin/picard.jar SamToFastq INPUT=${BamPath} FASTQ=${IntermedFilePath}${SampleID}.1.fq SECOND_END_FASTQ=${IntermedFilePath}${SampleID}.2.fq VALIDATION_STRINGENCY=LENIENT

### Align fastqs
use BWA
use Samtools
bwa mem -t 8 -R "@RG\\tID:FLOWCELL_${SampleID}\\tSM:${SampleID}\\tPL:ILLUMINA\\tLB:LIB_$SampleID" ${ReferenceFile} ${IntermedFilePath}${SampleID}.1.fq ${IntermedFilePath}${SampleID}.2.fq | samtools view -bS -> ${IntermedFilePath}${SampleID}.aligned.bam

### Sort Bam
use Java-1.8
java -Xmx8G -jar /cil/shed/apps/external/picard/current/bin/picard.jar SortSam I=${IntermedFilePath}${SampleID}.aligned.bam O=${IntermedFilePath}${SampleID}.sorted.bam SO=coordinate

### Mark Duplicates---generally needs at least -l h_vmem=10g
use Java-1.8
java -Xmx8G -jar /cil/shed/apps/external/picard/current/bin/picard.jar MarkDuplicates I=${IntermedFilePath}${SampleID}.sorted.bam O=${IntermedFilePath}${SampleID}.marked_duplicates.bam M=${IntermedFilePath}${SampleID}.marked_duplicates.metrics

### Reorder BAM
use Java-1.8
java -Xmx8G -jar /cil/shed/apps/external/picard/current/bin/picard.jar ReorderSam I=${IntermedFilePath}${SampleID}.marked_duplicates.bam O=${IntermedFilePath}${SampleID}.reordered.bam R=${ReferenceFile}
use Samtools
samtools index ${IntermedFilePath}${SampleID}.reordered.bam 

### GATK RealignerTargetCreator
use Java-1.8
java -Xmx8G -jar /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T RealignerTargetCreator -nct 1 -nt 24 -R ${ReferenceFile} -I ${IntermedFilePath}${SampleID}.reordered.bam -o ${IntermedFilePath}${SampleID}.interval_list
   
### GATK IndelRealigner
use Java-1.8
use Samtools
java -Xmx4G -jar /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T IndelRealigner -nct 1 -nt 1 -R ${ReferenceFile} -I ${IntermedFilePath}${SampleID}.reordered.bam -targetIntervals ${IntermedFilePath}${SampleID}.interval_list -o ${IntermedFilePath}${SampleID}.indels_realigned.bam
samtools index ${IntermedFilePath}${SampleID}.indels_realigned.bam
    
### GATK BaseRecalibrator
use Java-1.8
java -Xmx4G -jar /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -nt 1 -R ${ReferenceFile} -I ${IntermedFilePath}${SampleID}.indels_realigned.bam -knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/7g8_gb4.combined.final.vcf.gz -knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/hb3_dd2.combined.final.vcf.gz -knownSites /gsap/garage-protistvector/U19_Aim4/Pf3K/3d7_hb3.combined.final.vcf.gz -o ${IntermedFilePath}${SampleID}_recal_report.grp
    
### GATK Print Reads (BaseRecalibrator)
use Java-1.8
use Samtools
java -Xmx4G -jar /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T PrintReads -nct 8 -nt 1 -R ${ReferenceFile} -I ${IntermedFilePath}${SampleID}.indels_realigned.bam -BQSR ${IntermedFilePath}${SampleID}_recal_report.grp -o ${IntermedFilePath}${SampleID}.bqsr.bam
samtools index ${IntermedFilePath}${SampleID}.bqsr.bam

### GATK HaplotypeCaller
use Java-1.8
use Samtools
java -Xmx8G -jar /humgen/gsa-hpprojects/GATK/bin/GenomeAnalysisTK-3.5-0-g36282e4/GenomeAnalysisTK.jar -T HaplotypeCaller -nt 1 -R ${ReferenceFile} --input_file ${IntermedFilePath}${SampleID}.bqsr.bam -ERC GVCF -ploidy 2 --interval_padding 100 -o ${GvcfPath}${SampleID}.g.vcf -variant_index_type LINEAR -variant_index_parameter 128000

### Clean up

rm ${IntermedFilePath}${SampleID}.1.fq
rm ${IntermedFilePath}${SampleID}.2.fq
rm ${IntermedFilePath}${SampleID}.aligned.bam
rm ${IntermedFilePath}${SampleID}.sorted.bam
rm ${IntermedFilePath}${SampleID}.marked_duplicates.bam
rm ${IntermedFilePath}${SampleID}.marked_duplicates.metrics
rm ${IntermedFilePath}${SampleID}.reordered.bam
rm ${IntermedFilePath}${SampleID}.reordered.bam.bai
rm ${IntermedFilePath}${SampleID}.interval_list
rm ${IntermedFilePath}${SampleID}.indels_realigned.bam
rm ${IntermedFilePath}${SampleID}.indels_realigned.bam.bai
rm ${IntermedFilePath}${SampleID}.indels_realigned.bai
rm ${IntermedFilePath}${SampleID}_recal_report.grp
