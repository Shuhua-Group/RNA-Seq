# Written by Zhilin Ning

exp <- read.table("expr_detect.txt",head=TRUE)
sample <- colnames(exp)
sample <- sample[-1]
gene <- exp[,1]
exp_s <- subset(exp,select=-gene_id)
rownames(exp_s) <- gene
head(exp_s)
exp_ss <- exp_s

info <- read.table("ind_list_HU_info.txt")
pop <- as.vector(info$V2)
sex <- as.vector(info$V3)
age <- as.vector(info$V4)

color <-c(1)
for (i in pop)
if (i=="HAN"){color <- c(color,"red")} else{color <- c(color,"blue")}
cl <- color[-1]
pattern <- c(100)
for (i in sex)
if (i=="Male"){pattern <- c(pattern,0)} else{pattern <- c(pattern,1)}
ptn <- pattern[-1]

library(peer)
peer_function <- function(x,hidden){
    model = PEER()
    PEER_setPhenoMean(model,x)
    covs = cbind(pop, sex, age)
    rownames(covs)=rownames(x)
    colnames(covs)=c('pop','sex', 'age')
    PEER_setCovariates(model, as.matrix(covs))
    PEER_setNk(model,hidden)
    PEER_getNk(model)
    PEER_update(model)
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)
    results <- list(residuals=residuals,factors=factors,weights=weights,precision=precision)
    return(results)
}

re_3 <- peer_function(t(exp_ss),3)
result <- re_3$residuals
colnames(result) <- gene
rownames(result) <- sample
#pc_3 <- prcomp(re_3$residuals)
#compo1_3 <- pc_3$x[,1]
#compo2_3 <- pc_3$x[,2]

##re_4 <- peer_function(t(exp_ss),4)
##pc_4 <- prcomp(re_4$residuals)
##compo1_4 <- pc_4$x[,1]
##compo2_4 <- pc_4$x[,2]

##re_5 <- peer_function(t(exp_ss),5)
##pc_5 <- prcomp(re_5$residuals)
##compo1_5 <- pc_5$x[,1]
##compo2_5 <- pc_5$x[,2]

#pdf("after_peer_pca_h3.pdf")
#opar <- par(no.readonly=TRUE)
#par(mfrow = c(1, 1))
#plot(compo1_3,compo2_3,main="hidden factor is 3",col=cl,pch=ptn)
#legend("topleft",inset=.05,legend=c("HAN","UYG"),pch=c(20,20),col=c("red","blue"))
##text(compo1_3,compo2_3,sample,cex=0.5,pos=4,col="black")
##plot(compo1_4,compo2_4,main="hidden factor is 4",col=cl,pch=ptn)
##legend("topleft",inset=.05,legend=c("HAN","UYG"),pch=c(20,20),col=c("red","blue"))
##text(compo1_4,compo2_4,sample,cex=0.5,pos=4,col="black")
##plot(compo1_5,compo2_5,main="hidden factor is 5",col=cl,pch=ptn)
##legend("topleft",inset=.05,legend=c("HAN","UYG"),pch=c(20,20),col=c("red","blue"))
##text(compo1_5,compo2_5,sample,cex=0.5,pos=4,col="black")
#par(opar)
#dev.off()

write.table(result, file="expr_after_peer.txt")
