setwd("~/Documents/Well_Come_Trust_Project/Rotterdam/GEFOS/Result")

snpsList<-read.delim("SNP_Gene_10+16Gene_1MB_up_down_stream_WCT_sorted.txt", header=TRUE, sep="\t", strip.white=TRUE)
head(snpsList)
dim(snpsList)




Gefos<-read.delim("Homologus.txt", header=TRUE, sep="\t", strip.white=TRUE)
MergesnpsList <- merge(snpsList,Gefos,by="Region")

head(MergesnpsList)
dim(MergesnpsList)
write.table(MergesnpsList, "MergesnpsList.txt", sep="\t", row.names = F)



#Gefos<-read.delim("AOGC_only_chr.ALL.assoc", header=TRUE, sep="",, strip.white=TRUE)
#head(Gefos)
#dim(Gefos)



Gefos<-read.delim("GEFOS2.FNBMD.MEN.GC.txt", header=TRUE, sep="\t",, strip.white=TRUE)
Gefos<-read.delim("GEFOS2.FNBMD.WOMEN.GC.txt", header=TRUE, sep="")
Gefos<-read.delim("GEFOS2.FNBMD.POOLED.GC.txt", header=TRUE, sep="")
head(Gefos)
dim(Gefos)
colnames(Gefos)[1]<-"SNP"
#class(Gefos)

#snpsList<-MergesnpsList
snpsList<-Gene_250KB
head(snpsList)
##InnerJoin
Merge <- merge(snpsList,Gefos,by="SNP")

head(Merge)

dim(Merge)
write.table(Merge, "GEFOS2.FNBMD.MEN.GC_10+16Gene_250KB_up_down_stream_WCT_sorted.txt", sep="\t", row.names = F)

require(graphics)
x<-p.adjust(p=Gefos$P.value, method = "bonferroni", n = length(Gefos$P.value))
head(x)
x
.............
#############
#PName<-read.delim("mmuPathIdPathName.txt", header=TRUE, sep="\t",, strip.white=TRUE)
#Merge <- merge(PName,data,by="pathwayId..mmu00010.")

#head(Merge)






tab5rows <- read.table("snp142.txt", sep="\t", nrows = 10, ncolumns = 10)
classes <- sapply(tab5rows, class)
tabAll <- read.table("snp142.txt", , colClasses = classes)
df <- read.table(pipe("cut -f1,2,3,4,5 snp142.txt"))

head(df)
write.table(df, "Allsnp152.txt", sep="\t", row.names = F)
