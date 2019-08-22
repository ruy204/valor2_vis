#Visualization_R

library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(Gviz)
library(gridExtra)

#################
#   Functions   #
#################

#1. check specific region on chromosome

chromosome_check<-function(chr,start_position,chunck_length,annotate="population",types=unique(bedpe$type)){
  
  bedpe<-read.delim("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",header = F)
  colnames(bedpe)<-c("chrom1","start1","end1","chrom2","start2","end2",
                     "type","score","number1","number2","number3",
                     "samples","Description")
  bedpe[,c(2,3,5,6,8,9,10,11)]<-sapply(bedpe[,c(2,3,5,6,8,9,10,11)],function(x){as.numeric(as.character(x))})
  bedpe$chrom1<-factor(bedpe$chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$chrom2<-factor(bedpe$chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe<-bedpe %>% dplyr::filter(type %in% types)
  
  end_position=start_position + chunck_length #start_position*0.1  #16300000
  
  #1. create genome axis track and ideogram track
  
  gtrack<-GenomeAxisTrack()
  # plotTracks(gtrack,from=0,to=250000000)
  gen="hg19"
  names(gen)=chr
  itrack<-IdeogramTrack(genome = gen,chromosome=chr)
  
  #2. create gene track (requires network!)
  
  biomTrack <- BiomartGeneRegionTrack(genome = "hg19",
                                      chromosome = chr, start = 0, end = 250000000,
                                      name = "ENSEMBL",transcriptAnnotation="symbol", background.title="brown")
  
  #3. create data track
  
  lim<-c(start_position,end_position)
  subdata = bedpe %>% dplyr::filter(chrom1==chr, start1 %in% seq(start_position,end_position,1))
  dat<-as.numeric(subdata$type)
  coords <- sort(c(subdata$start1,lim[2]))
  
  dtrack<-DataTrack(data=dat,start=coords[-length(coords)], end=coords[-1], chromosome=chr,
                    genome=gen, name="SV",cex=1.5)
  
  aTrack<-AnnotationTrack(start=subdata$start1,width = abs(subdata$end1-subdata$start1),chromosome = chr,
                          group=subdata$samples,genome="hg19",name="SV",
                          stacking = "squish")
  
  if(annotate=="population"){
    identifier(aTrack)=as.character(subdata$Description)
  }else{
    identifier(aTrack)=as.character(subdata$samples)
  }
  
  
  plotTracks(list(itrack,gtrack,biomTrack,aTrack), col.line = NULL, col = NULL,
             from=start_position,to=end_position,extend.left = 0.1,extend.right = 0.1,
             background.panel="#FFFEDB", background.title="darkblue",groupAnnotation="group")
}


























