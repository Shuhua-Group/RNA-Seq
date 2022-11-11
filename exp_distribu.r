# Written by Zhilin Ning

exp <- read.table("expr_merge.txt",head=TRUE)
sample <- colnames(exp)
gene <- exp[,1]
exp_s <- subset(exp,select=-gene_id)
rownames(exp_s) <- gene
head(exp_s)
s1 <- exp_s[,1]
#usually add 1 then log2
s1_2 <- log(s1+1,2)
pdf("sample_1_exp_distribution.pdf")
plot(density(s1_2))
xfit <- seq(min(s1_2),max(s1_2),length=15)
yfit <- dnorm(xfit, mean=mean(s1_2),sd=sd(s1_2))
lines(xfit,yfit,col="red")
dev.off()

