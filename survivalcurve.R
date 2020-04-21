library(survival)
library(survminer)
kmsurvo<-Surv(combined$Rectime,combined$Censor.Status)
sfit <- survfit(Surv(combined$Rectime,combined$Censor.Status)~PTPRT, data=combined)
ggsurvplot(sfit, pval=TRUE, 
           legend.labs=c("Normal", "Altered"), legend=c(.75,.75),  
           title="BAMBI")
ggsurvplot(sfit, legend = c(0.2, 0.2))
plot(sfit)

ggsurv<-ggsurvplot(sfit, size = 1.5, 
           pval=0.022,pval.size=8, #pval.coord=c(0,.05),#
           legend.title="", legend.labs=c("Normal", "Altered"),
           xlab = "Time in days",
           font.x = c(18, "bold.italic", "darkred"),
           font.y = c(18, "bold.italic", "darkred"),
           legend=c(.75,.75),  
           title="PTPRT", main = "Survival curve",
           font.main = c(22, "bold", "darkblue"),
           font.tickslab = c(16, "plain", "darkgreen"))

ggsurv$plot <- ggsurv$plot + 
  theme(legend.text = element_text(size = 18, color = "black", face = "bold"))
ggsurv

