# Expression profiles of east-west highly differentiated genes in Uyghur genomes 

## Introduction <br> 
Here is the pipeline describing the data processing in the manuscript

## RNA-seq data processing 
### 1. Quality assessment 
```shell
fastqc ${fq_dir}/${sampleID}_R1.fastq.gz ${fq_dir}/${sampleID}_R2.fastq.gz --outdir $wd/Quality_Assessment/${sampleID}
```
### 2. trim + QC 
```shell
mkdir -p ${fq_dir}/trim
cd ${fq_dir}/trim
trim_galore -q 20 --trim1 --paired --fastqc ${fq_dir}/${sampleID}_R1.fastq.gz ${fq_dir}/${sampleID}_R2.fastq.gz
```
### 3. mapping 
#### Building the STAR index 
```shell
fasta="$wd/ref/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
gtf="$wd/ref/Homo_sapiens.GRCh37.75.gtf"
mkdir $wd/ref/$STARgenomeDir

STAR --runThreadN 5 --runMode genomeGenerate --genomeDir $STARgenomeDir --genomeFastaFiles $fasta --sjdbGTFfile $gtf --limitGenomeGenerateRAM=50503721344 --limitSjdbInsertNsj 3000000 --limitIObufferSize 6368709120 --sjdbOverhang 99 --outFileNamePrefix $wd/$STARgenomeDir
``` 
#### Mapping with STAR
```shell
mkdir -p $wd/mapping

sh star_mapping.sh ${fq_dir}/trim/${sampleID}_R1.fastq.gz ${fq_dir}/trim/${sampleID}_R2.fastq.gz $wd/ref/$STARgenomeDir  ${num_thread} $wd/mapping
```
### 4. Quantify expression level with RSEM
#### Building the RSEM index
```shell
rsem-prepare-reference $wd/ref/ref_ensemble/ $wd/ref/ref_RSEM/Homo_sapiens.GRCh37.75 --gtf $gtf -p 5
```
#### RSEM running
```shell
sh RSEM.sh ${sampleID} $wd/mapping $wd/ref/ref_RSEM/Homo_sapiens.GRCh37.75 ${num_thread} $wd/RSEM/
```
### 5. Merging and QC of the FPKM matrix
#### Merging the RSEM results to FPKM matrix
```shell
python /picb/humpopg-bigdata5/ningzhilin/RNA_Seq/exp_old/expr_merge.py
```
#### FPKM QC
```shell
Rscript exp_distribu.r
python expr_detect.py
Rscript exp_distribu_filter.r
```
#### PCA analysis of FPKM matrix
```shell
Rscript exp_pca.r
```
#### Normalization of FPKM matrix with PEER
```shell
Rscript exp_peer.r
```

## Functional analysis
### 1. ASE related analysis
Calling of ASE and aseQTL were applied with the Python package ["ASEkit"](https://pypi.org/project/ASEkit/) developed by our team.
### 2. QTL analysis
The format of input data is consistant with examples in R package ["MatrixeQTL"](http://bios.unc.edu/research/genomic_software/Matrix_eQTL/)
```shell
Rscript MatrixEQTL.r $geno $exp $cov $gene_loc $snploc $cis_res $trans_res
```
### 3. Roadmap enrichment
The input are two-columns, TAB-separated, BED-format files
```shell
python2 roadmap.enrich.py --cis ${cis_loci} --bg ${bg_loci} --out ${prefix}
```s
