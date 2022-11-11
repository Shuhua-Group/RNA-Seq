# Written by Xinjiang Tan

library(clusterProfiler)
library(msigdf)
library(dplyr)
library(org.Hs.eg.db)
vignette("msigdf")

args <- commandArgs(T)

cis_file <- args[1]
bg_file <- args[2]
out_name <- args[3]

cis_genes <- read.table(cis_file,stringsAsFactors = F)
cis_gene.df <- bitr(cis_genes[,1],fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)

bg_genes <- read.table(bg_file,stringsAsFactors = F)
bg_gene.df <- bitr(bg_genes[,1],fromType = "ENSEMBL",toType = c("SYMBOL"),OrgDb = org.Hs.eg.db)

c2 <- msigdf.human %>% filter(category_code=='c2')
c2 <- data.frame(geneset=c2[,3],symbol=c2[,4])

res <- enricher(cis_gene.df$SYMBOL,universe=bg_gene.df$SYMBOL,TERM2GENE = c2)

pdf(out_name,width = 16,height = 9)
dotplot(res,showCategory=50)
dev.off()
