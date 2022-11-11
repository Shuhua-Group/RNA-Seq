# Written by Zhilin Ning

exp <- read.table("expr_detect.txt",head=TRUE)
sample <- colnames(exp)
gene <- exp[,1]
exp_s <- subset(exp,select=-gene_id)
rownames(exp_s) <- gene
head(exp_s)
s1 <- exp_s[,1]

pdf("sample_1_exp_distribution_filtered.pdf")
plot(density(s1))
xfit <- seq(min(s1),max(s1),length=15)
yfit <- dnorm(xfit, mean=mean(s1),sd=sd(s1))
lines(xfit,yfit,col="red")
dev.off()

