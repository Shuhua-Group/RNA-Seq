# Written by Zhilin Ning

#exp <- read.table("expr_merge.txt",head=TRUE)
exp <- read.table("expr_merge_unstrand.txt",head=TRUE)
gene <- exp[,1]
exp_s <- subset(exp,select=-gene_id)
rownames(exp_s) <- gene
sample <- colnames(exp_s)
head(exp_s)
pc_result = prcomp(t(exp_s))
compo1 <- pc_result$x[,1]
compo2 <- pc_result$x[,2]

info <- read.table("ind_list_HU_info.txt")
pop <- as.vector(info$V2)
sex <- as.vector(info$V3)

color <-c(1)
for (i in pop)
if (i=="HAN"){color <- c(color,"red")} else{color <- c(color,"blue")}
cl <- color[-1]
pattern <- c(100)
for (i in sex)
if (i=="Male"){pattern <- c(pattern,0)} else{pattern <- c(pattern,1)}
ptn <- pattern[-1]

#pdf("exp_pca.pdf")
pdf("exp_pca_unstrand.pdf")
opar <- par(no.readonly=TRUE)
par(mfrow = c(1, 1))
plot(compo1,compo2,col=cl,pch=ptn)
text(compo1,compo2,sample,cex=0.4,pos=2,col="black")
legend("topleft",inset=.05,legend=c("HAN","UYG"),pch=c(20,20),col=c("red","blue"))
par(opar)
dev.off()
