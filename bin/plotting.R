lambda<-function(pvalues){
chisq <- qchisq(1-pvalues,1)
lambda=median(chisq,na.rm=TRUE)/qchisq(0.5,1)
return(lambda)
}


#Plotting manhattan plots
library(data.table)
ANA_C1_v1_filt<-fread(SAIGEOUTPUT)

library(qqman)

colnames(ANA_C1_v1_filt)[which(colnames(ANA_C1_v1_filt)=="POS")]="BP"
colnames(ANA_C1_v1_filt)[which(colnames(ANA_C1_v1_filt)=="p.value")]="P"

png("ANA_C1_v1_filt_manhattan.png",width=12,height=6,units='in',res=300)

p1<-manhattan(ANA_C1_v1_filt, main = "ANA_C1_v1", cex = 0.6,
    cex.axis = 0.6, col = c("blue4", "orange3"), suggestiveline = FALSE)

dev.off()

#qqplot
library(GWASTools)

png("ANA_C1_v1_filt_qqplotCI.png")

qqPlot(ANA_C1_v1_filt$P, main=paste0("ANA_C1_v1",",","lambda=",round(lambda(ANA_C1_v1_filt$P),2)))

dev.off()

#check whether we can just save the R objects and also give the option for exporting through a future airlock
