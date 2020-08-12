lambda<-function(pvalues){
chisq <- qchisq(1-pvalues,1)
lambda=median(chisq,na.rm=TRUE)/qchisq(0.5,1)
return(lambda)
}

#Loading required packages

library(data.table)
library(GWASTools)
library(qqman)


#Plotting manhattan plots
args = commandArgs(trailingOnly=TRUE)
SAIGE_outputfile=args[1]
analysis_tag=args[2]

analysis<-fread(SAIGE_outputfile)

colnames(analysis)[which(colnames(analysis)=="POS")]="BP"
colnames(analysis)[which(colnames(analysis)=="p.value")]="P"

png(paste0(analysis_tag,"_manhattan.png"),width=12,height=6,units='in',res=300)
p1<-manhattan(analysis, main = analysis_tag, cex = 0.6, cex.axis = 0.6, col = c("blue4", "orange3"), suggestiveline = FALSE)
dev.off()

#qqplot

png(paste0(analysis_tag,"_qqplotCI.png"))
qqPlot(analysis$P, main=paste0(analysis_tag,",","lambda=",round(lambda(analysis$P),2)))
dev.off()

#saveRDS(analysis, file = paste0(analysis_tag,".Rdata"))
#check whether we can just save the R objects and also give the option for exporting through a future airlock
