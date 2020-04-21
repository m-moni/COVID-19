setwd("F:\\Green University\\Personal\\Research\\Covid-19\\Covid-19 and lung cancer survival")
#################################################################
##############################working with clinical data#########
#################################################################
dato<-read.csv("lusc_tcga_clinical_data.csv",header = TRUE,stringsAsFactors = FALSE)

library(dplyr)
datos<-select(dato,Sample.ID,Race.Category,Diagnosis.Age,Tumor.Tissue.Site,American.Joint.Committee.on.Cancer.Tumor.Stage.Code,Patient.ID,Patient.Primary.Tumor.Site,Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,Censor.Status,Overall.survival.in.days)
#datos<-select(dato,Sample.ID,Race.Category,Diagnosis.Age,Tumor.Tissue.Site,American.Joint.Committee.on.Cancer.Tumor.Stage.Code,Patient.ID,Patient.Primary.Tumor.Site,Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,Censor.Status,Overall.survival.in.days)
#colnames(dato)[c(82,83)]

#ranme the column
#datos<-dplyr::rename(datos,Patient_ID=Patient.ID,anatomic_site=Tumor.Disease.Anatomic.Site,histologic_grade=Neoplasm.Histologic.Grade,Rectime=Overall.survival.in.days,age=Diagnosis.Age,tumour_site=Primary.Tumor.Site,cancer_stage=Neoplasm.American.Joint.Committee.on.Cancer.Clinical.Group.Stage,race=Race.Category)
datos<-dplyr::rename(datos,Patient_ID=Patient.ID,anatomic_site=Tumor.Tissue.Site,histologic_grade=American.Joint.Committee.on.Cancer.Tumor.Stage.Code,Rectime=Overall.survival.in.days,age=Diagnosis.Age,tumour_site=Patient.Primary.Tumor.Site,cancer_stage=Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,race=Race.Category)
#replace the "" in race category
datos$race=as.character(lapply(datos$race,function(x){gsub("^$","Others",x)}))
#replace the "" in cancer_site with NA
datos$tumour_site=as.character(lapply(datos$tumour_site,function(x){gsub("^$",NA,x)}))
#replace "" in cancer_stage with NA
#datos$cancer_stage=as.character(lapply(datos$cancer_stage,function(x){gsub("^$",NA,x)}))
datos$anatomic_site=as.character(lapply(datos$anatomic_site,function(x){gsub("^$",NA,x)}))
datos$histologic_grade=as.character(lapply(datos$histologic_grade,function(x){gsub("^$",NA,x)}))
#datos$cancer_stage=as.character(lapply(datos$,function(x){gsub("^$",NA,x)}))

library(plyr)
count(datos,'cancer_stage')
count(datos,'race')
count(datos,'tumour_site')
count(datos,'anatomic_site')
count(datos,'histologic_grade')


#################################################################
##############################working with mRNA data#########
#################################################################

GSE147507_Result <- read.csv("~/Bangladesh/Corona/GSE147507_Result_sur.csv")
list<-GSE147507_Result[,1]
name(list)<-Hugo_Symbol


#read rna expression data
expr_mrA<-read.csv("data_mRNA_median_Zscores.csv",header=FALSE,stringsAsFactors = FALSE)
expr_mr<-read.csv("data_mRNA_median_Zscores.csv",header=TRUE,stringsAsFactors = FALSE)
colnames(expr_mr)<-expr_mrA[1,]

dim(expr_mr)
class(expr_mr)

expr_mr2<-expr_mr[expr_mr[,1] %in% GSE147507_Result[,1],]
rownames(expr_mr2)<-expr_mr2[,1]
expr_mr2<-expr_mr2[,-1]
head(expr_mr2)
class(expr_mr2)
dim(expr_mr2)
AVG<-colMeans(expr_mr2)
expr_mr2<-rbind(expr_mr2, AVG)

#function for labelling each expression value
altered_test<-function(x){{
  if (x>=1.5){d="Cancer"
    }else d="Control"
   
    }
  d
 }


#altered_test(10)


#apply the function over all the colulmn to convert altered unaltered

applyfunc<-function(df,f){
  ds<-matrix(0,nrow = nrow(df),ncol=ncol(df))
  colnames(ds)<-colnames(df)
  for (i in seq(1:ncol(df))){
    ds[,i]<-(sapply(df[,i],f))
  }
  ds<-as.data.frame(ds)
}
gene_status<-applyfunc(expr_mr2,altered_test)

row.names(gene_status)<-row.names(expr_mr2)
#remove the 01 from patient iD
write.csv(t(gene_status), "gene_status.csv")
remove_01<-function(x){
  x<-unlist(strsplit(x,split=""))
  x<-paste(x[0:(length(x)-3)],collapse = "")
  x
}

#remove_01("TCGA-04-1332-01")

gene_status$Patient_ID<-as.character(gene_status$Patient_ID)
gene_status$Patient_ID=unlist(lapply(gene_status$Patient_ID,remove_01))
#View(gene_status)
####################################################################
#########################Merge the tables###########################
###################################################################
#gene_status$Patient_ID=as.character(gene_status$Patient_ID)
combined<-datos%>%inner_join(gene_status)
#View(combined)  
#View(datos)
#View(gene_status)

#relevel the genes as normal as reference factor
applyrevel<-function(combined){
  
  col_names<-colnames(combined)[11:ncol(combined)]
  for(i in col_names){
    #combined$i<-as.factor(combined$i)
    combined[,i]<-as.factor(combined[,i])
    combined[,i]<-relevel(combined[,i],ref="Control")
    #combined$i<-relevel(fctocombined$i,ref="Normal")
  }
  combined
}
combined<-applyrevel(combined)
######################################################
################ Univariate analysis ################
#####################################################

library(survival)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
applycox<-function(combined){
  models<-list()
  col_names<-colnames(combined)[11:ncol(combined)]
  for(i in col_names){
    
    fit<-coxph(kmsurvo~factor(combined[,i]),data=combined)
    tss<-summary(fit)
    coefs<-c(tss$coefficients[c(1,2,5)])
    models[[i]]=coefs
  }
  
  final_mode<-as.data.frame(models) 
  final_model=t(final_mode)
  colnames(final_model)<-c("coef","exp.coef","p")
  
  as.data.frame(final_model)
}
fs<-applycox(combined)

#class(fs)
fs<-fs%>%mutate(gene=rownames(.))
fs
######################################################
################ Multivariate Analysis################
#####################################################

fitt<-coxph(kmsurvo~.,data=combined[,11:ncol(combined)])
fitt
######################################################
##Clinical and rna Expression survival analysis#######
#####################################################

#fit_age<-coxph(kmsurvo~age,data=combined)
#lapply(combined,getClass)

#there is some "" in Race ,replace that with "others"
# library(methods)
# d<-c("fruit","Saiful")
# d<-as(d,'list')
# lapply(d, function(x){gsub("f","F",x)})

#combined$Race=as.character(lapply(combined$Race,function(x){gsub("^$","Others",x)}))
#View(combined)

#unique(datos[,7])

combined_reor<-combined[,c(6,9,10,2:5,7,8,11:ncol(combined))]
#View(combined_reor)
#kmsurvo1<-Surv(combined_reor$Rectime,combined_reor$Censor.Status)
#fitt_grand1<-coxph(kmsurvo1~.-1,data=combined_reor[,4:ncol(combined_reor)])
combined_reor$race=factor(combined_reor$race)
fitt_grand<-coxph(kmsurvo~.-1,data=combined_reor[,4:ncol(combined_reor)])

fitt_grand

