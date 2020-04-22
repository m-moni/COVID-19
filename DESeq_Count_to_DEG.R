library(DESeq2)
library(DESeq)
# see vignette for suggestions on generating
# count tables from RNA-Seq data
data<-read.delim("~/USyd/NuroDiseases/Alzheimers/T2D/GSE107489_gene_counts.txt", header =TRUE, row.names=1)
head(data)
row.names(data) <- make.names(data[,1],TRUE)
head(data)
data<-as.matrix(data[,-1])
#data<-data[,-1]



#cnts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10)
cnts<-as.matrix(data)


#cond <- factor(rep(1:2, each=5))
#cond <- factor(c(2, rep(2:1, each=5)))
cond <- factor(rep(1:2, c(10,22)) )



# object construction
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)

###Filtering
keep <- rowSums(counts(cnts)) >= 10  ##no of row usually
dds <- dds[keep,]
# standard analysis
dds2 <- DESeq(dds)
res <- results(dds2)

# an alternate analysis: likelihood ratio test
ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
resLRT <- results(ddsLRT)
write.csv(res,"E-MTAB-8871_res.csv" )
write.csv(resLRT,"E-MTAB-8871_resLRT.csv" )


##############
#EDGER
library(edgeR)
x<-data
> y <- DGEList(counts=x)
A grouping factor can be added at the same time:
  group <- c(1,1,1, 2,2, 2)
y <- DGEList(counts=x, group=group)

######
> x <- read.delim("TableOfCounts.txt",row.names="Symbol")
x<-data
> group <- factor(c(1,1,2,2))
> y <- DGEList(counts=x,group=group)
> y <- calcNormFactors(y)
> design <- model.matrix(~group)
> y <- estimateDisp(y,design)
To perform quasi-likelihood F-tests:
  > fit <- glmQLFit(y,design)
> qlf <- glmQLFTest(fit,coef=2)
result1<-topTags(qlf, n=Inf)
To perform likelihood ratio tests:
  > fit <- glmFit(y,design)
> lrt <- glmLRT(fit,coef=2)
> topTags(lrt)
result2<-topTags(lrt, n=Inf)

