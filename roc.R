##ROCcurve그리기###

getwd()
setwd("/home/gangtl95/data/ROC/")

data1 <- read.table("./FAC_NWCOAD_seg.txt",header=T,sep="\t",stringsAsFactors = F)
data2 <- read.table("./NEX_NWCOAD_seg.txt",header=T,sep="\t",stringsAsFactors = F)
data3 <- read.table("./SEQ_NWCOAD_seg.txt",header=T,sep="\t",stringsAsFactors = F)
#"4N-A93T-01A","A6-2677-01A","A6-6652-01A","AA-3655-01A","AA-3848-01A","AA-3854-01A","CK-6746-01A","CM-5862-01A","QG-A5YX-01A","SS-A7HO-01A"

library(ggplot2)
library(ROCR)

dim(data1)
dim(data2)
dim(data3)

table(data1$OL.Ratio,data1$GS.Result)
table(data2$OL.Ratio,data2$GS.Result)
table(data3$OL.Ratio,data3$GS.Result)

tiff(filename=paste("ROCcurve.","SEQ_Ploidy2",".tiff",sep=""),width=3072,height=2304,res=300)
pred <- prediction(data1$OL.Ratio,data1$GS.Result)
perf <- performance(pred,"tpr","fpr")  #(수행 data, measure, 2nd measure)ROC 곡선의 경우 measure = tpr , x.measure= fpr
par(bg="white",mai=c(1.2,1.5,1,1),pch=0)
plot(perf,col="seagreen3",lwd=3.5,
     print.cutoffs.at=seq(0,1,by=0.3),text.adj=NA,xlim=c(0,1),xaxis.cex.axis=1.5,
     yaxis.cex.axis=1.5,xaxis.lwd=2,yaxis.lwd=2, cex.lab=2)  #colorize 는 설정한 cutoff 에 따라 곡선의 색이 변함



par(pch=1)
pred2 <- prediction(data2$OL.Ratio,data2$GS.Result)
perf2 <- performance(pred2,"tpr","fpr")
plot(perf2,add=TRUE,col="goldenrod1",lwd=3.5,
     print.cutoffs.at=seq(0,1,by=0.3),text.adj=NA)



par(pch=2)
pred3 <- prediction(data3$OL.Ratio,data3$GS.Result)
perf3 <- performance(pred3,"tpr","fpr")
plot(perf3,add=TRUE,col="darkorchid1",lwd=3.5,
     print.cutoffs.at=seq(0,1,by=0.3),text.adj=NA)
abline(0,1,lty=2)
legend(0.67,0.28,c('FACETS','NEXUS','SEQUENZA'),col=c("seagreen3", "goldenrod1", "darkorchid1"),lwd=6,cex = 1.4)
invisible(dev.off())

auc.perf=performance(pred,measure="auc")
fac=auc.perf@y.values
fac 
auc.perf2=performance(pred2,measure="auc")
nex=auc.perf2@y.values
nex 
auc.perf3=performance(pred3,measure="auc")
seq=auc.perf3@y.values
seq

auc.all.info=data.frame(fac[[1]],nex[[1]],seq[[1]],row.names = "AUC")
auc.all.info
colnames(auc.all.info)=c("FACETS","NEXUS","SEQUENZA")
auc.all.info



##Specificity, sensitivity 구하기


data1[data1$OL.Ratio==1]
# install.packages("remotes")
# remotes::install_github("cardiomoon/multipleROC")
#https://bookdown.org/cardiomoon/roc/intro.html#tab-5
#devtools::install_github("cardiomoon/webr") >> 안됨
library(multipleROC)
#library(webr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
x=multipleROC(GS.Result~OL.Ratio,data=data3)  #x=OL.Ratio, y=GS.Result

calSens=function(x,y){
  newx=sort(c(unique(x),max(x,na.rm=TRUE)+1))
  completeTable=function(res){
    if(nrow(res)==1){
      res1=matrix(c(0,0),nrow=1)
      temp=setdiff(c("TRUE","FALSE"),attr(res,"dimnames")[[1]][1])
      if(temp=="FALSE") res=rbind(res1,res)
      else res=rbind(res,res1)
      res
    }
    res
  }
  
  getSens=function(cut){
    res=table(x>=cut,y)
    res=completeTable(res)
    sens=res[2,2]/sum(res[,2])
    spec=res[1,1]/sum(res[,1])
    ppv=res[2,2]/sum(res[2,])
    npv=res[1,1]/sum(res[1,])
    data.frame(x=cut,sens=sens,spec=spec,fpr=1-spec,ppv=ppv,npv=npv,sum=sens+spec)
  }
  map_dfr(newx,function(cut){getSens(cut)})
}

result=calSens(data3$OL.Ratio,data3$GS.Result)
result

longdf <- result %>% pivot_longer(cols=sens:spec,names_to = "rate")
ggplot(data=longdf,aes(x=x,y=value,color=rate))+geom_line()
result[which.max(result$sum),]
