# mp=.libPaths()
# mp=c(mp,"/home/gangtl95/R/x86_64-pc-linux-gnu-library/4.0")

#Single Chromosome 
if(!require(dplyr)){install.packages("dplyrr");library(dplyr)}
if(!require(CopyNumberPlots)){BiocManager::install("CopyNumberPlots");library(CopyNumberPlots)}
if(!require(GenomicRanges)){BiocManager::install("GenomicRanges");library(GenomicRanges)}

getwd()
setwd("/home/gangtl95/karyoplot/")

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
cent.38

#hg19 cent
# cent19=read.table("./hg19centlen.txt",header=F,sep="\t",stringsAsFactors = F)
# colnames(cent19)=c("chrom","chromStart","chromEnd","name","gieStain")
# cent19=cent19 %>% filter(chrom !="chrX"& chrom !="chrY"& chrom !="chrM")
# cent.19=as.data.frame(unique(cent19$chrom))
# colnames(cent.19)="chrom"
# cent.19=cent.19 %>% mutate(centStart=1) %>% mutate(centEnd=1)
# 
# 
# for (i in 1:length(cent.19$chrom)){
#   cent.19[i,2]=cent19[2*i-1,2]
#   cent.19[i,3]=cent19[2*i,3]
# }
# 
# cent.19

#hg38 length
wh_len.38=read.table("./data/hg38_chromlen.txt",sep="\t",stringsAsFactors = F,
                     header=T);wh_len.38
wh_len.38=wh_len.38 %>% mutate(Start=1,End=Length)
wh_len.38=wh_len.38[,c(1,3,4,2)]
wh_len.38=wh_len.38 %>% filter(Chr!="chrX" & Chr!="chrY" & Chr!="chrMT");wh_len.38

#hg19 length
# wh_len.19=read.table("./hg19_chromlen.txt",sep="\t",stringsAsFactors = F,
#                 header=F)
# colnames(wh_len.19)=c("chrom","chromStart","chromEnd","name")
# wh_len.19=as.data.frame(wh_len.19)
# wh_len.19=wh_len.19 %>% mutate(Start=chromStart+1,End=chromEnd,Length=chromEnd)
# wh_len.19=wh_len.19[,c(1,5,6,7)]
# wh_len.19=wh_len.19 %>% filter(chrom!="chrX" & chrom!="chrY" & chrom!="chrMT");wh_len.19

#초기 설정값 
#지수형태의 숫자를 full로 나열 
#반대로 하고 싶은 경우 -100
#"4N-A93T-01A","A6-2677-01A","A6-6652-01A","AA-3655-01A","AA-3848-01A","AA-3854-01A","CK-6746-01A","CM-5862-01A","QG-A5YX-01A","SS-A7HO-01A"

for (N in 1:22){
  for ( code in c("4N-A93T-01A","A6-2677-01A","A6-6652-01A","AA-3655-01A","AA-3848-01A","AA-3854-01A","CK-6746-01A","CM-5862-01A","QG-A5YX-01A","SS-A7HO-01A")){
    
    #hg38로 수정#cent.19 <-> cent.38 #wh_len.19 <-> wh_len.38
    # "GRCh38" <-> "GRCh37" # "hg38" <-> "hg19"
    # code="SS-A7HO-01A"
    # N=1
    grcgh="GRCh38"
    gb="hg38"
    smpl_n=paste("TCGA",code,sep="-")
    fil_n=paste("COAD",code,sep="-") #file name
    chr=paste("chr",N,sep="")
    sT=wh_len.38[N,2];sT #chromosome start site
    eNd=wh_len.38[N,3];eNd #chromosome end site
    smain=paste(chr,":",sT,"-",eNd,sep="") #subtitle
    a=paste("NW",smpl_n,"_",gb,"_",chr,".tiff",sep="")
    print(N)
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
    
    #hg38_bam
    #paste("./",gb,"/",fil_n,"_",gb,"_lrr.txt",sep="")
    lrr=read.table(paste("./data/",gb,"/",fil_n,"_",gb,"_lrr.txt",sep=""),
                   header=T,stringsAsFactors = F,sep="\t")
    lrr=lrr[,1:4]
    head(lrr)
    colnames(lrr)=c("Chr", "ST", "END", "Log.R.Ratio")
    lrr$Log.R.Ratio=as.numeric(lrr$Log.R.Ratio)
    lrr=lrr %>% filter(lrr$ST>=sT)
    lrr=lrr %>% filter(lrr$Chr==chr)
    lrr=lrr %>% filter(lrr$END<=eNd)
    probes.data <- toGRanges(lrr)
    probes.data
    names(mcols(probes.data)) <- "lrr"
    seqlevelsStyle(probes.data) <- "UCSC"
    tail(probes.data)
    
    baf=read.table(paste("./data/",gb,"/",fil_n,"_",gb,"_baf.txt",sep=""),
                   header=T,stringsAsFactors=F,sep="\t")
    baf=baf[,1:4]
    colnames(baf)=c("Chr", "ST", "END", "B.Allele.Freq")
    baf=baf %>% filter(baf$Chr==chr)
    baf=baf %>% filter(baf$ST>=sT)
    baf=baf %>% filter(baf$END <= eNd)
    baf$B.Allele.Freq=as.numeric(baf$B.Allele.Freq)
    snps.data <- toGRanges(baf)
    snps.data
    names(mcols(snps.data)) <- "baf"
    seqlevelsStyle(snps.data) <- "UCSC"
    tail(snps.data)
    
    
    
    ##STADARD Results##
    loh_a=read.table(paste("./data/len/hg38_ar/",fil_n,"_ar.txt",sep=""),
                     sep="\t",stringsAsFactors = F,header=T);loh_a
    loh_a=loh_a %>% filter(loh_a$Chr!="chr0")  %>% filter(loh_a$Chr==chr)
    if(is.na(loh_a[1,1])==TRUE){
      loh_a=spp
    }
    loh_a=loh_a %>% filter(loh_a$Chr==chr) %>% mutate(loh=1) #뒤에 loh 갯수에 대한 field 추가 (default)
    loh_a=loh_a[,c(3:5,7)]
    colnames(loh_a)=c("Chr", "ST", "END", "loh");loh_a
    loh_a.data=toGRanges(loh_a);loh_a.data
    
    ##FACETS Results##
    loh_f=read.table(paste("./data/len/hg38_fa/",fil_n,"_fa.txt",sep=""),
                     sep="\t",stringsAsFactors = F,header=T);loh_f
    loh_f=loh_f %>% filter(loh_f$Chr!="chr0") %>% filter(loh_f$Chr==chr)
    if(is.na(loh_f[1,1])==TRUE){
      loh_f=spp
    }
    loh_f=loh_f %>% filter(loh_f$Chr==chr) %>% mutate(loh=1)
    loh_f=loh_f[,c(3:5,7)]
    colnames(loh_f)=c("Chr", "ST", "END", "loh")
    loh_f.data=toGRanges(loh_f);loh_f.data
    
    ##NEXUS Results##
    loh_n=read.table(paste("./data/len/hg38_nx/",fil_n,"_nx.txt",sep=""),
                     sep="\t",stringsAsFactors = F,header=T);loh_n
    loh_n=loh_n %>% filter(loh_n$Chr!="chr0") %>% filter(loh_n$Chr==chr)
    if(is.na(loh_n[1,1])==TRUE){
      loh_n=spp
    }
    loh_n=loh_n %>% filter(loh_n$Chr==chr) %>% mutate(loh=1)
    loh_n=loh_n[,c(3:5,7)]
    colnames(loh_n)=c("Chr", "ST", "END", "loh")
    loh_n.data=toGRanges(loh_n);loh_n.data
    
    ##SEQUENZA Results##
    loh_s=read.table(paste("./data/len/hg38_sq/",fil_n,"_sq.txt",sep=""),
                     sep="\t",stringsAsFactors = F,header=T);loh_s
    loh_s=loh_s %>% filter(loh_s$Chr!="chr0") %>% filter(loh_s$Chr==chr)
    if(is.na(loh_s[1,1])==TRUE){
      loh_s=spp
    }
    loh_s=loh_s %>% filter(loh_s$Chr==chr) %>% mutate(loh=1)
    loh_s=loh_s[,c(3:5,7)];loh_s
    colnames(loh_s)=c("Chr", "ST", "END", "loh")
    loh_s.data=toGRanges(loh_s);loh_s.data
    
    
    ##Plot Main title, subtitle labeling##
    #plotKaryotype
    #getDefaultPlotParams
    #https://bernatgel.github.io/karyoploter_tutorial//Tutorial/PlotParams/PlotParams.html
    plot.params.4 <- list(leftmargin = 0.1, rightmargin = 0.01, 
                          topmargin = 30, bottommargin = -20, ideogramheight = 10, 
                          ideogramlateralmargin = 0.003, data1height = 200, 
                          data1inmargin = 10, data1outmargin = 0, data1min = 0, 
                          data1max = 1, data2height = 0, data2inmargin = 0, 
                          data2outmargin = 0, data2min = 0, data2max = 1) #adjusting plot parameters plot.type=4
    
    #kpAddMainTitle
    #내 방식대로 수정 (kpAddMainTitle_m)
    kpAddMainTitle_m=function (karyoplot, main = NULL, ...) 
    {
      if (!is.null(main)) {
        karyoplot$beginKpPlot()
        on.exit(karyoplot$endKpPlot())
        bb <- getMainTitleBoundingBox(karyoplot)
        x <- (bb$x0 + bb$x1)/2
        y <- (bb$y0 + bb$y1+10)/2
        graphics::text(x = x, y = y, labels = main, ...)
      }
      invisible(karyoplot)
    }
    #내 방식대로 수정 (kpAddSubTitle)
    kpAddSubTitle=function (karyoplot, smain = NULL, ...) 
    {
      if (!is.null(smain)) {
        karyoplot$beginKpPlot()
        on.exit(karyoplot$endKpPlot())
        bb <- getMainTitleBoundingBox(karyoplot)
        x <- (bb$x0 + bb$x1)/2
        y <- (bb$y0 + bb$y1-14)/2
        graphics::text(x = x, y = y, labels = smain, ...)
      }
      invisible(karyoplot)
    }
    
    tiff(filename=a,width=3072,height=2304,res=300)
    ##Plot Title (main, sub)
    kp=plotKaryotype(genome=gb,chromosomes = chr, plot.params = plot.params.4)
    #kpAddMainTitle_m(kp, main=smpl_n, cex=1.5)
    #kpAddSubTitle(kp, smain=paste(smain," (",grcgh,"/",gb,")",sep=""), cex=1.2)
    kpAddMainTitle_m(kp, main=smain, cex=2)
    
    ##Plot BAF #second panel
    kpDataBackground(kp, r0=0.325, r1=0.675)
    #kpAbline(kp, h=0.675, col="darkgray")
    kpAbline(kp, v=c(cent.38[N,2],cent.38[N,3]), r0=0.325, r1=0.675,col="tomato",lty=2,lwd=3)  #centromere pos
    plotBAF(kp, snps.data, labels=NA, r0=0.325, r1=0.675,add.axis = F)
    bmax=1
    bmin=0
    f2=function(var){
      ag2=(0.675-0.325)/(bmax-bmin)
      ag2*(var-bmax)+0.675
    }
    kpAbline(kp,h=c(f2(0.85),f2(0.15)),col="blue",lty=1,lwd=2) #baf thresholod line 
    kpAxis(kp, ymin=0, ymax=1, r0=0.325, r1=0.675, numticks=5, tick.pos=c(1,0.85,0.5,0.15,0),
           labels=NA,cex=0.6,lwd=2) #tick.pos에서의 위치는 축을 기준으로 한다. 
    #c("1","0.85","0.5","0.15","0")
    #kpAddLabels(kp, labels ="BAF", srt=0, pos=2, label.margin = 0.05, cex=1.3, r0=0.6,r1=0.625)
    
    ##Plot LRR #first panel
    kpDataBackground(kp, r0=0.7, r1=1.05,data.panel=2)
    plotLRR(kp,probes.data,labels=NA,ymin=-2,ymax =2,r0=0.7, r1=1.05,add.axis = F)
    lmax=2
    lmin=-2
    f1=function(var){
      ag1=(1.05-0.7)/(lmax-lmin);ag1
      ag1*(var-lmax)+1.05
    }
    
    kpAbline(kp, v=c(cent.38[N,2],cent.38[N,3]), r0=0.7, r1=1.05,col="tomato",lty=2,lwd=3) #centromere pos
    kpAbline(kp,h=c(f1(-1.1),f1(-0.15),f1(0.1),f1(0.7)),col="blue",lty=1,lwd=2) #lrr thresholod line
    kpAxis(kp, r0=0.7, r1=1.05,ymin=-2,ymax =2, numticks=7, tick.pos=c(2,0.7,0.1,0,-0.15,-1.1,-2),
           labels=NA,cex=0.6,lwd=2) #labels=c("2","0.7 (High Gain)","0.1 (Gain)   ","0","-0.15 (Loss)   ","-1.1 (Big Loss)","-2")
    #kpAddLabels(kp, labels = "LRR", srt=0, pos=2, label.margin = 0.05, cex=1.3, r0=0.99,r1=1)
    ##Bottome Panel STANDARD## 
    kpDataBackground(kp, r0=0.1, r1=0.3)
    
    for ( i in 1:length(loh_a$Chr)){
      xs=loh_a$ST[i]
      xe=loh_a$END[i]
      if(xs==1 && xe==1){
        kpRect(kp, chr = chr,x0=0,x1=0,y0=0.25, y1=0.27,col="tomato",lwd=0)
      }else{
        kpRect(kp, chr = chr,x0=xs,x1=xe,y0=0.25, y1=0.27,col="tomato",lwd=1)
      }
    } #Results Bar expression
    for ( i in 1:(length(loh_a$Chr)-1)){
      xs=loh_a$ST[i+1]
      xe=loh_a$END[i]
      print (xs)
      print (xe)
      dif=xs-xe
      print (dif)
      if (length(loh_a$Chr)==1){
        print("Only one segments")
      }else if(dif>0&dif <6.5e5){
        print("Multiple segments")
        kpPoints(kp,chr=chr,x=(xe+xs)/2,y=0.24,pch=17,col="black",cex=1.1)
      }else{
        print("Interval is long")
      }
    }
    kpAddLabels(kp, labels="STANDARD", r0=0.25, r1=0.27,cex=1) #Labeling
    
    ##Bottome Panel FACETS##
    for ( i in 1:length(loh_f$Chr)){
      xs=loh_f$ST[i]
      xe=loh_f$END[i]
      if(xs==1 && xe==1){
        kpRect(kp, chr = chr,x0=0,x1=0,y0=0.21, y1=0.23,col="tomato",lwd=0)
      }else{
        kpRect(kp, chr = chr,x0=xs,x1=xe,y0=0.21, y1=0.23,col="tomato",lwd=1)
      }
    }
    for ( i in 1:(length(loh_f$Chr)-1)){
      xs=loh_f$ST[i+1]
      xe=loh_f$END[i]
      print (xs)
      print (xe)
      dif=xs-xe
      print (dif)
      if (length(loh_f$Chr)==1){
        print("Only one segments")
      }else if(dif>0&dif <6.5e5){
        print("Multiple segments")
        kpPoints(kp,chr=chr,x=(xe+xs)/2,y=0.2,pch=17,col="black",cex=1.1)
      }else{
        print("Interval is long")
      }
    }
    kpAddLabels(kp, labels="FACETS", r0=0.21, r1=0.23,cex=1)
    
    
    ##Bottome Panel NEXUS##
    for ( i in 1:length(loh_n$Chr)){
      xs=loh_n$ST[i]
      xe=loh_n$END[i]
      if(xs==1 && xe==1){
        kpRect(kp, chr = chr,x0=0,x1=0,y0=0.17, y1=0.19,col="tomato",lwd=0)
      }else{
        kpRect(kp, chr = chr,x0=xs,x1=xe,y0=0.17, y1=0.19,col="tomato",lwd=1)
      }
    }
    for ( i in 1:(length(loh_n$Chr)-1)){
      xs=loh_n$ST[i+1]
      xe=loh_n$END[i]
      print (xs)
      print (xe)
      dif=xs-xe
      print (dif)
      if (length(loh_n$Chr)==1){
        print("Only one segments")
      }else if(dif>0&dif <6.5e5){
        print("Multiple segments")
        kpPoints(kp,chr=chr,x=(xe+xs)/2,y=0.16,pch=17,col="black",cex=1.1)
      }else{
        print("Interval is long")
      }
    }
    kpAddLabels(kp, labels="NEXUS", r0=0.17, r1=0.19, cex=1)
    
    ##Bottome Panel SEQUENZA##
    for ( i in 1:length(loh_s$Chr)){
      xs=loh_s$ST[i]
      xe=loh_s$END[i]
      if(xs==1 && xe==1){
        kpRect(kp, chr = chr,x0=0,x1=0,y0=0.13, y1=0.15,col="tomato",lwd=0)
      }else{
        kpRect(kp, chr = chr,x0=xs,x1=xe,y0=0.13, y1=0.15,col="tomato",lwd=1)
      }
    }
    for ( i in 1:(length(loh_s$Chr)-1)){
      xs=loh_s$ST[i+1]
      xe=loh_s$END[i]
      print (xs)
      print (xe)
      dif=xs-xe
      print (dif)
      if (length(loh_s$Chr)==1){
        print("Only one segments")
      }else if(dif>0&dif <6.5e5){
        print("Multiple segments")
        kpPoints(kp,chr=chr,x=(xe+xs)/2,y=0.12,pch=17,col="black",cex=1.1)
      }else{
        print("Interval is long")
      }
    }
    kpAddLabels(kp, labels="SEQUENZA", r0=0.13, r1=0.15, cex=1)
    invisible(dev.off())
    wh_len.38
    print(paste(a," is finished",sep=""))
  }
}
