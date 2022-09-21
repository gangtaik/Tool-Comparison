# mp=.libPaths()
# mp=c(mp,"/home/gangtl95/R/x86_64-pc-linux-gnu-library/4.0")

if(!require(dplyr)){install.packages("dplyrr");library(dplyr)}
if(!require(CopyNumberPlots)){BiocManager::install("CopyNumberPlots");library(CopyNumberPlots)}
if(!require(GenomicRanges)){BiocManager::install("GenomicRanges");library(GenomicRanges)}
install.packages("devtools")
if(!require(svglite)){devtools::install_github("r-lib/svglite");library(svglite)}

getwd() 
setwd("/home/gangtl95/karyoplot/SVG/")

# centromere position
# Go to the Table Browser: http://genome.ucsc.edu/cgi-bin/hgTables
# Choose the Mapping and Sequencing group
# Select the "Chromosome Band (Ideogram)" track
# Select filter, and enter "acen" in the gieStain field
# Press "submit" and then "get output"
# Each chromosome will have two entries which overlap.
# https://groups.google.com/a/soe.ucsc.edu/g/genome/c/SaR2y4UNrWg/m/XsGdMI3AazgJ
# filter: cytoBandIdeo.gieStain like '%cen'
#chrom	chromStart	chromEnd	name	gieStain

#hg38 cent
cent38=read.table("./data/hg38cenlen.txt",header=F,sep="\t",stringsAsFactors = F)
colnames(cent38)=c("chrom","chromStart","chromEnd","name","gieStain")
cent38=cent38 %>% filter(chrom !="chrX"& chrom !="chrY"& chrom !="chrM")
cent.38=data.frame(unique(cent38$chrom))
colnames(cent.38)="chrom"
cent.38=cent.38 %>% mutate(centStart=1) %>% mutate(centEnd=1)
for (i in 1:length(cent.38$chrom)){
  cent.38[i,2]=cent38[2*i-1,2]
  cent.38[i,3]=cent38[2*i,3]
}

wh_len.38=read.table("./data/hg38_chromlen.txt",sep="\t",stringsAsFactors = F,
                     header=T);wh_len.38
wh_len.38=wh_len.38 %>% mutate(Start=1,End=Length)
wh_len.38=wh_len.38[,c(1,3,4,2)]
wh_len.38=wh_len.38 %>% filter(Chr!="chrX" & Chr!="chrY" & Chr!="chrMT");wh_len.38
pt.rg=GRanges(wh_len.38)
#hg38로 수정#cent.19 <-> cent.38 #wh_len.19 <-> wh_len.38
# "GRCh38" <-> "GRCh37" # "hg38" <-> "hg19"
#"4N-A93T-01A","A6-2677-01A","A6-6652-01A","AA-3655-01A","AA-3848-01A","AA-3854-01A","CK-6746-01A","CM-5862-01A","QG-A5YX-01A","SS-A7HO-01A"
for (code in c("4N-A93T-01A","A6-2677-01A","A6-6652-01A","AA-3655-01A","AA-3848-01A","AA-3854-01A","CK-6746-01A","CM-5862-01A","QG-A5YX-01A","SS-A7HO-01A")){
  code="QG-A5YX-01A"
  grcgh="GRCh38"
  gb="hg38"
  smpl_n=paste("TCGA",code,sep="-")
  fil_n=paste("COAD",code,sep="-") #file name
  a=paste("SVG.Wh_",smpl_n,"_",gb,".svg",sep="")
  chr="chr5"
  print(code)
  print(smpl_n) 
  print(fil_n)
  print(a)
  
  
  #빈 dataframe에 채울 supple data 만들기##
  spp=matrix(c(fil_n,"tool",chr,1,1,1,1),nrow=1)
  spp=as.data.frame(spp,stringsAsFactors = F)
  colnames(spp)=c("Sample", "Tool","Chr", "ST", "END","Length", "loh")
  spp$ST=as.numeric(spp$ST)
  spp$END=as.numeric(spp$END)
  spp$Length=as.numeric(spp$Length)
  spp$loh=as.numeric(spp$loh)
  
  ##lrr&baf data 만들기##
  lrr=read.table(paste("./data/",gb,"/",fil_n,"_",gb,"_lrr.txt",sep=""),
                 header=T,stringsAsFactors = F,sep="\t")
  lrr=lrr[,1:4]
  colnames(lrr)=c("Chr", "ST", "END", "Log.R.Ratio")
  lrr$Log.R.Ratio=as.numeric(lrr$Log.R.Ratio)
  # lrr=lrr %>% filter(lrr$ST>=sT)
  # lrr=lrr %>% filter(lrr$Chr==chr)
  # lrr=lrr %>% filter(lrr$END<=eNd)
  probes.data <- toGRanges(lrr)
  probes.data
  names(mcols(probes.data)) <- "lrr"
  seqlevelsStyle(probes.data) <- "UCSC"
  tail(probes.data)
  
  baf=read.table(paste("./data/",gb,"/",fil_n,"_",gb,"_baf.txt",sep=""),
                 header=T,stringsAsFactors=F,sep="\t")
  baf=baf[,1:4]
  colnames(baf)=c("Chr", "ST", "END", "B.Allele.Freq")
  # baf=baf %>% filter(baf$Chr==chr)
  # baf=baf %>% filter(baf$ST>=sT)
  # baf=baf %>% filter(baf$END <= eNd)
  baf$B.Allele.Freq=as.numeric(baf$B.Allele.Freq)
  snps.data <- toGRanges(baf)
  snps.data
  names(mcols(snps.data)) <- "baf"
  seqlevelsStyle(snps.data) <- "UCSC"
  tail(snps.data)
  
  
  
  #https://www.bioconductor.org/packages/devel/bioc/vignettes/karyoploteR/inst/doc/karyoploteR.html
  #plot parameter info
  plot.params.4 <- list(leftmargin = 0.05, rightmargin = 0.01, 
                        topmargin = 30, bottommargin = -20, ideogramheight = 10, 
                        ideogramlateralmargin = 0.003, data1height = 200, 
                        data1inmargin = 10, data1outmargin = 0, data1min = 0, 
                        data1max = 1, data2height = 0, data2inmargin = 0, 
                        data2outmargin = 0, data2min = 0, data2max = 1)
  chr_lis=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22")
  chr_num=gsub("chr","",chr_lis)
  
  #svglite(filename="1.svg",width=44,height=8) #inches
  #tiff(filename=a,width=1700,height=300,units = "mm", res=300)
  pdf(file ="test.pdf" )
  kp <- plotKaryotype(genome=gb,plot.type = 4,ideogram.plotter = NULL,plot.params = plot.params.4,
                      chromosomes=chr_lis,labels.plotter=NULL)
  kpDataBackground(kp, r0=0.73, r1=1.12,data.panel=2,color=c("#A5E5FF","white"))
  plotLRR(kp,probes.data,labels=NA,ymin=-2,ymax =2,r0=0.73, r1=1.12,add.axis = F,lwd=3)
  lmax=2
  lmin=-2
  f1=function(var){
    ag1=(1.12-0.73)/(lmax-lmin)
    ag1*(var-lmax)+1.12
  }
  
  #kpAbline(kp, v=c(cent.38[N,2],cent.38[N,3]), r0=0.675, r1=1,col="tomato",lty=2,lwd=3) #centromere pos
  kpAbline(kp,h=c(f1(-1.1),f1(-0.15),f1(0),f1(0.1),f1(0.7)),col="blue",lty=1,lwd=3) #lrr thresholod line
  #kpAbline(kp,h=c(0.73,1.12),,col="black",lty=1,lwd=0.5)
  kpAxis(kp, r0=0.73, r1=1.12,ymin=-2,ymax =2, numticks=7, tick.pos=c(2,0.7,0.1,0,-0.15,-1.1,-2),cex=1,labels=NA,lwd=3)
  
  #labels=c("2","0.7","0.1\t","0","-0.15\t","-1.1","-2")
  
  #kpAddChromosomeNames(kp,chr.names = chr_num)
  kpDataBackground(kp, r0=0.33, r1=0.7,color=c("#A5E5FF","white")) #grey90, #A5E5FF
  plotBAF(kp, snps.data, labels=NA, r0=0.33, r1=0.7,add.axis = F)
  bmax=1
  bmin=0
  f2=function(var){
    ag2=(0.7-0.33)/(bmax-bmin)
    ag2*(var-bmax)+0.7
  }
  kpAbline(kp,h=c(f2(0.85),f2(0.15)),col="blue",lty=1,lwd=3) #baf thresholod line 
  kpAxis(kp, ymin=0, ymax=1, r0=0.33, r1=0.7, numticks=5, tick.pos=c(1,0.85,0.5,0.15,0),cex=1,labels=NA,lwd=3)
  #labels=c("1","0.85","0.5","0.15","0")
  
  
  
  
  
  
  ##Bottome Panel for STANDRAD and Tool Results## 
  kpDataBackground(kp, r0=0.02, r1=0.30,color=c("#A5E5FF","white"))
  kpAbline(kp,h=c(0.253,0.197,0.141,0.085),col="black",lty=1,lwd=3.5) #baf thresholod line 
  
  ##STADARD Results##
  loh_a=read.table(paste("../../len/hg38_ar/",fil_n,"_ar.txt",sep=""),
                   sep="\t",stringsAsFactors = F,header=T);loh_a
  loh_a=loh_a %>% mutate(loh=1) #뒤에 loh 갯수에 대한 field 추가 (default)
  loh_a
  if(is.na(loh_a[1,1])==TRUE){
    loh_a=spp
  }
  loh_a=loh_a[,c(3:5,7)]
  colnames(loh_a)=c("Chr", "ST", "END", "loh");loh_a
  
  for ( i in 1:length(loh_a$Chr)){
    xs=loh_a$ST[i]
    xe=loh_a$END[i]
    kpRect(kp, chr=loh_a[i,1],x0=xs,x1=xe,y0=0.253, y1=0.273,col="tomato")
  }
  
  ##FACETS Results##
  loh_f=read.table(paste("../../len/hg38_fa/",fil_n,"_fa.txt",sep=""),
                   sep="\t",stringsAsFactors = F,header=T);loh_f
  loh_f=loh_f %>% mutate(loh=1)
  if(is.na(loh_f[1,1])==TRUE){
    loh_f=spp
  }
  loh_f=loh_f[,c(3:5,7)]
  colnames(loh_f)=c("Chr", "ST", "END", "loh")
  
  for ( i in 1:length(loh_f$Chr)){
    xs=loh_f$ST[i]
    xe=loh_f$END[i]
    kpRect(kp, chr=loh_f[i,1],x0=xs,x1=xe,y0=0.197, y1=0.227,col="tomato")
  }
  ##NEXUS Results##
  loh_n=read.table(paste("../../len/hg38_nx/",fil_n,"_nx.txt",sep=""),
                   sep="\t",stringsAsFactors = F,header=T);loh_n
  loh_n=loh_n %>% mutate(loh=1)
  if(is.na(loh_n[1,1])==TRUE){
    loh_n=spp
  }
  loh_n=loh_n[,c(3:5,7)]
  colnames(loh_n)=c("Chr", "ST", "END", "loh")
  
  for ( i in 1:length(loh_n$Chr)){
    xs=loh_n$ST[i]
    xe=loh_n$END[i]
    kpRect(kp, chr=loh_n[i,1],x0=xs,x1=xe,y0=0.141, y1=0.171,col="tomato")
  }
  
  ##SEQUENZA Results##
  loh_s=read.table(paste("../../len/hg38_sq/",fil_n,"_sq.txt",sep=""),
                   sep="\t",stringsAsFactors = F,header=T);loh_s
  loh_s=loh_s %>% mutate(loh=1)
  if(is.na(loh_s[1,1])==TRUE){
    loh_s=spp
  }
  loh_s=loh_s[,c(3:5,7)];loh_s
  colnames(loh_s)=c("Chr", "ST", "END", "loh")
  
  for ( i in 1:length(loh_s$Chr)){
    xs=loh_s$ST[i]
    xe=loh_s$END[i]
    kpRect(kp, chr=loh_s[i,1],x0=xs,x1=xe,y0=0.085, y1=0.11,col="tomato")
  }
  
  invisible(dev.off())
}