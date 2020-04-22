
#############################GDS Download
#GDS download and processing test
source("http://bioconductor.org/biocLite.R")
biocLite()
library(GEOquery)
##________________
data<- "GSE57195"

gse<-"GSE57195"
gds <- getGEO(data)
#Convert to ExpressionSet
eset <- GDS2eSet(gds,do.log2=FALSE)
data<-data.frame(fData(eset),exprs(eset))
Gene.symbol
write.csv(,file=paste(gse,".csv",sep="",collapse = NULL),row.names=FALSE)
paste(gse,".txt",sep="",collapse = NULL)
####################


gse<-"GSE51693_2"

data<- "GDS5249"

gds <- getGEO(data)
data<-Table(gds)
data<-data[!is.na(data$IDENTIFIER),]

write.csv(data,file=paste(gse,".csv",sep="",collapse = NULL),row.names=FALSE)
####################GSE
gse<-"GSE33075"

gset <-getGEO(gse, GSEMatrix =TRUE, AnnotGPL=TRUE)
idx<-1
gset <- gset[[idx]]
#ex <- exprs(gset)


feset <- fData(gset)
teset <- exprs(gset)
	
cnames <- c( c("ID_REF", "IDENTIFIER"), colnames(teset))

data<-data.frame(feset$"ID", feset$"Gene symbol",teset)
colnames(data)<-cnames
#attach(data)
data<-data[!is.na(data$IDENTIFIER),]
write.csv(data,file=paste(gse,".csv",sep="",collapse = NULL),row.names=FALSE)

#############################GSE Download
#GDS download and processing test
source("http://bioconductor.org/biocLite.R")
biocLite()
library(GEOquery)
gset <-getGEO("GSE57195", GSEMatrix =FALSE, AnnotGPL=TRUE)

if (length(gset) > 1) idx <- grep("GPL16570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#Convert to ExpressionSet
 ex <- exprs(gset)
 qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
 LogC <- (qx[5] > 100) ||
   (qx[6]-qx[1] > 50 && qx[2] > 0) ||
   (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
 if (LogC) { ex[which(ex <= 0)] <- NaN
 exprs(gset) <- log2(ex) }
write.csv(data.frame(GeneSymbol=fData(gset)$"Gene symbol",exprs(gset)),file="GSE57195.csv",row.names=FALSE)


#ma%%%%%%%%Differentially Expression Analysis
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Aug 7 22:20:59 EDT 2018

#############################################Differential expression analysis with limma###################
#   Differential expression analysis with limma
#Examples

# Simulate gene expression data for 100 probes and 6 microarrays
# Microarray are in two groups
# First two probes are differentially expressed in second group
# Std deviations vary between genes with prior df=4
sd <- 0.3*sqrt(4/rchisq(100,df=4))
y <- matrix(rnorm(100*6,sd=sd),100,6)
rownames(y) <- paste("Gene",1:100)
y[1:2,4:6] <- y[1:2,4:6] + 2
design <- cbind(Grp1=1,Grp2vs1=c(0,0,0,1,1,1))
options(digits=3)
#####Data
data<-GSE00000
rownames(data) = make.names(data[,1], unique=TRUE)
data<-data[,-1]
design <- cbind(Grp1=1,Grp2vs1=c(1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0))


# Ordinary fit
fit <- lmFit(data1[,1:17],design)
fit <- eBayes(fit)
topTable(fit,coef=2)
dim(fit)
colnames(fit)
rownames(fit)[1:10]
names(fit)

# Fold-change thresholding
fit2 <- treat(fit,lfc=0.1)
topTreat(fit2,coef=2)

# Volcano plot
volcanoplot(fit,coef=2,highlight=2)

# Mean-difference plot
plotMD(fit,column=2)

# Q-Q plot of moderated t-statistics
qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)

# Various ways of writing results to file
## Not run: write.fit(fit,file="exampleresults.txt")
## Not run: write.table(fit,file="exampleresults2.txt")

# Fit with correlated arrays
# Suppose each pair of arrays is a block
block <- c(1,1,2,2,3,3)
dupcor <- duplicateCorrelation(y,design,block=block)
dupcor$consensus.correlation
fit3 <- lmFit(y,design,block=block,correlation=dupcor$consensus)

# Fit with duplicate probes
# Suppose two side-by-side duplicates of each gene
rownames(y) <- paste("Gene",rep(1:50,each=2))
dupcor <- duplicateCorrelation(y,design,ndups=2)
dupcor$consensus.correlation
fit4 <- lmFit(y,design,ndups=2,correlation=dupcor$consensus)
dim(fit4)
fit4 <- eBayes(fit4)
topTable(fit4,coef=2)
