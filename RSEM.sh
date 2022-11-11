#!/bin/bash

# RSEM quantification pipeline
# usage: from an empty working directory, run
#requires following softwares:
##1) samtools
##2) rsem
##3) R
##4) ucsc_tools
#input: BAM file from STAR
read1=$1
RSEMrefDir=$2
dataType="unstr_PE" # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsRSEM=$3 # number of threads for RSEM
workDIR=$4
cd $workDIR
#name=`basename $read1 .1.trimmed.fastq.gz` #for trimmed
name=`basename $read1 _combined_filtered_R1_val_1.fq.gz` #for val
mkdir $name
cd $name
echo $name

# output: all in the working directory, fixed names
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Quant.genes.results                           # RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                        # RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                     # RSEM diagnostic plots

# executables
RSEM=rsem-calculate-expression        
bedGraphToBigWig=bedGraphToBigWig              

######### RSEM

#### prepare for RSEM: sort transcriptome BAM to ensure the order of the reads, to make RSEM output (not pme) deterministic
trBAMsortRAM=60G
cp /picb/humpopg-bigdata5/ningzhilin/RNA_Seq/mapping_unstrand/"${name}"/"${name}".Aligned.toTranscriptome.out.bam ./raw.bam
cp raw.bam Tr.bam 

case "$dataType" in
str_SE|unstr_SE)
      # single-end data
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | sort -S $trBAMsortRAM -T ./ ) | samtools view -@ $nThreadsRSEM -bS - > $name.Aligned.toTranscriptome.out.bam
      ;;
str_PE|unstr_PE)
      # paired-end data, merge mates into one line before sorting, and un-merge after sorting
      cat <( samtools view -H Tr.bam ) <( samtools view -@ $nThreadsRSEM Tr.bam | awk '{printf "%s", $0 " "; getline; print}' | sort -S $trBAMsortRAM -T ./ | tr ' ' '\n' ) | samtools view -@ $nThreadsRSEM -bS - > $name.Aligned.toTranscriptome.out.bam
      ;;
esac

'rm' Tr.bam


# RSEM parameters: common
RSEMparCommon="--bam --estimate-rspd --no-bam-output --seed 12345"

# RSEM parameters: run-time, number of threads and RAM in MB
RSEMparRun=" -p $nThreadsRSEM "

# RSEM parameters: data type dependent

case "$dataType" in
str_SE)
      #OPTION: stranded single end
      RSEMparType="--forward-prob 0"
      ;;
str_PE)
      #OPTION: stranded paired end
      RSEMparType="--paired-end --forward-prob 0"
      ;;
unstr_SE)
      #OPTION: unstranded single end
      RSEMparType=""
      ;;
unstr_PE)
      #OPTION: unstranded paired end
      RSEMparType="--paired-end"
      ;;
esac


###### RSEM command
echo $RSEM $RSEMparCommon $RSEMparRun $RSEMparType $name.Aligned.toTranscriptome.out.bam $RSEMrefDir $name >& $name.Log.rsem
$RSEM $RSEMparCommon $RSEMparRun $RSEMparType $name.Aligned.toTranscriptome.out.bam $RSEMrefDir $name >& $name.Log.rsem

###### RSEM diagnostic plot creation
# Notes:
# 1. rsem-plot-model requires R (and the Rscript executable)
# 2. This command produces the file $name.pdf, which contains multiple plots
echo rsem-plot-model $name $name.pdf
rsem-plot-model $name $name.pdf

