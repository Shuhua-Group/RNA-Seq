#!/bin/bash
# Written by Lei Tian

# STAR mapping / RSEM quantification pipeline
# usage: from an empty working directory, run
# qsub -V -cwd -l h_vmem=4G -pe shm 12 -j y -N $name ./STAR_RSEM_ENCODE_test.sh (read1) (read2 or "") (STARgenomeDir) (RSEMrefDir) (dataType) (nThreadsSTAR) (nThreadsRSEM)
#Required following softwares:
##1) STAR
##2) samtools
##3) R
##4) ucsc_tools
# input: gzipped fastq file read1 [read2 for paired-end] 
#        STAR genome directory, RSEM reference directory - prepared with STAR_RSEM_prep.sh script
read1=$1 #gzipped fastq file for read1
read2=$2 #gzipped fastq file for read1, use "" if single-end
STARgenomeDir=$3 # STAR reference index
dataType="unstr_PE" # RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=$4 # number of threads for STAR
workDIR=$5
#Change to working directory
cd $workDIR
#name=`basename $read1 .1.fastq.gz` #for untrimmed
#name=`basename $read1 .1.trimmed.fastq.gz` #for trimmed
#name=`basename $read1 _val_1.fq.gz` #for val
#name=`basename $read1 _combined_filtered_R1_val_1.fq.gz` #trimmed
name=`basename $read1`
name=${name%%_*}
mkdir $name
cd $name

# output: all in the working directory, fixed names
# Aligned.sortedByCoord.out.bam                 # alignments, standard sorted BAM, agreed upon formatting
# Log.final.out                                 # mapping statistics to be used for QC, text, STAR formatting
# Signal.{Unique,UniqueMultiple}.strand{+,-}.bw # 4 bigWig files for stranded data
# Signal.{Unique,UniqueMultiple}.unstranded.bw  # 2 bigWig files for unstranded data

# executables
STAR=STAR                             
bedGraphToBigWig=bedGraphToBigWig              


# STAR parameters: common
STARparCommon=" --genomeDir $STARgenomeDir  --readFilesIn $read1 $read2 --outFileNamePrefix "${name}". --outSAMunmapped Within --outFilterType BySJout \
 --outSAMattributes NH HI AS NM MD    --outFilterMultimapNmax 20   --outFilterMismatchNmax 999   \
 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   \
 --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat "

# STAR parameters: run-time, controlled by DCC
#STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad LoadAndKeep  --limitBAMsortRAM 10000000000"
STARparRun=" --runThreadN $nThreadsSTAR --limitBAMsortRAM 10000000000"

# STAR parameters: type of BAM output: quantification or sorted BAM or both
#     OPTION: sorted BAM output
## STARparBAM="--outSAMtype BAM SortedByCoordinate"
#     OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#     OPTION: both
STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"


# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

case "$dataType" in
str_SE|str_PE)
      #OPTION: stranded data
      STARparStrand=""
      STARparWig="--outWigStrand Stranded"
      ;;
      #OPTION: unstranded data
unstr_SE|unstr_PE)
      STARparStrand="--outSAMstrandField intronMotif"
      STARparWig="--outWigStrand Unstranded"
      ;;
esac

# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile commentsENCODElong.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"

## not needed ## --outSAMheaderPG @PG ID:Samtools PN:Samtools CL:"$samtoolsCommand" PP:STAR VN:0.1.18"

# ENCODE metadata BAM comments
#echo -e '@CO\tLIBID:ENCLB175ZZZ
#@CO\tREFID:ENCFF001RGS
#@CO\tANNID:gencode.v19.annotation.gtf.gz
#utSA    
#@CO\tSPIKEID:ENCFF001RTP VN:Ambion-ERCC Mix, Cat no. 445670' > commentsENCODElong.txt

###### STAR command
echo $STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta
$STAR $STARparCommon $STARparRun $STARparBAM $STARparStrand $STARparsMeta

###### bedGraph generation, now decoupled from STAR alignment step
# working subdirectory for this STAR run
mkdir Signal

bamfilename=""${name}".Aligned.sortedByCoord.out.bam"

echo $STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bamfilename --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

$STAR --runMode inputAlignmentsFromBAM --inputBAMfile $bamfilename --outWigType bedGraph $STARparWig --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr

# move the signal files from the subdirectory
mv Signal/Signal*bg .




###### bigWig conversion commands
# exclude spikeins
grep ^chr $STARgenomeDir/chrNameLength.txt > chrNL.txt

case "$dataType" in
str_SE|str_PE)
      # stranded data
      str[1]=-; str[2]=+;
      for istr in 1 2
      do
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str$istr.out.bg > sig.tmp
          $bedGraphToBigWig sig.tmp  chrNL.txt Signal.$imult.strand${str[istr]}.bw
      done
      done
      ;;
unstr_SE|unstr_PE)
      # unstranded data
      for imult in Unique UniqueMultiple
      do
          grep ^chr Signal.$imult.str1.out.bg > sig.tmp
          $bedGraphToBigWig sig.tmp chrNL.txt  Signal.$imult.unstranded.bw
      done
      ;;
esac


