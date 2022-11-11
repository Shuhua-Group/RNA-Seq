## Load genotype data
library(MatrixEQTL)

snps = SlicedData$new();
snps$fileDelimiter = '\t'; # the TAB character
snps$fileOmitCharacters = 'NA'; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
snps$LoadFile("/picb/humpopg-bigdata5/ningzhilin/RNA_new/plink/genotype_QC.txt");

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = '\t'; # the TAB character
gene$fileOmitCharacters = 'NA'; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
gene$LoadFile("/picb/humpopg-bigdata5/ningzhilin/RNA_new/exp/GE_FPKM_PEERed_s.txt");

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = 2000; # read file in one piece
cvrt$LoadFile("/picb/humpopg-bigdata5/ningzhilin/RNA_new/eQTL/UYG90/covariates_s.txt");


gene.loc=read.table('/picb/humpopg-bigdata5/ningzhilin/RNA_new/exp/geneloc.txt',header = 1, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", as.is = 1)


snp.loc=read.table('/picb/humpopg-bigdata5/ningzhilin/RNA_new/plink/snpsloc.txt',header = 1, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", as.is = 1)
me_ln= Matrix_eQTL_main(snps=snps,gene=gene,cvrt=cvrt,output_file_name='trans.1M.txt',pvOutputThreshold=1e-6,useModel=modelLINEAR,pvOutputThreshold.cis=1,output_file_name.cis='cis_fpkm.1M.txt',snpspos=snp.loc,genepos=gene.loc,cisDist=1e6)
