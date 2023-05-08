
install.packages("CMplot")
library(CMplot)
setwd("D:\\Rdemo\\CMplot")

data<-read.csv('dataY.csv')
LR<-data[,10]
datap7<-pchisq(LR,6,lower.tail = F)
LR<-data[,7]
datap4<-pchisq(LR,6,lower.tail = F)
datag<-data[,c(1:3)]
data<-cbind(datag,datap1,datap2,datap3,datap4,datap5,datap6,datap7)

write.csv(data, file = "G:/??R?ĵط?/P??.csv",row.names = FALSE, col.names = FALSE)

datapp<-data[,c(1:6)]
datatpp<-data[,c(1,2,3,7,8)]
dataxpp<-data[,c(1,2,3,9,10)]


library(RColorBrewer)

CMplot(dataxpp,plot.type="c",
       r=0.4,col=c(c(brewer.pal(5,"Pastel1"))),
       chr.labels=paste("Chr",c(1:5),sep=""),
       cir.chr.h=1.5,
       amplify=TRUE,threshold.lty=c(1,2),
       threshold.col=c("red","blue"),signal.line=1,
       signal.col=c("red","green"),
       chr.den.col=c("darkgreen","yellow","red"),
       bin.size=1e6,outward=FALSE,file="tiff",
       memo="quan",dpi=300,file.output=TRUE,verbose=TRUE)



