#Visualization_R

library(dplyr)
library(tidyr)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(GenomicRanges)
library(Gviz)
library(gridExtra)
library(Sushi)
library(BioCircos)

# color extraction function
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

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

#2. check intra-chromosome events for specific chromosome/region

intra_check<-function(chromosome,start,end,plot_type="loops"){
  
  #read-in bedpe data
  bedpe<-read.delim("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",header = F)
  colnames(bedpe)<-c("chrom1","start1","end1","chrom2","start2","end2",
                     "type","score","number1","number2","number3",
                     "samples","Description")
  bedpe[,c(2,3,5,6,8,9,10,11)]<-sapply(bedpe[,c(2,3,5,6,8,9,10,11)],function(x){as.numeric(as.character(x))})
  bedpe$chrom1<-factor(bedpe$chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$chrom2<-factor(bedpe$chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$length1<-bedpe$end1-bedpe$start1
  bedpe$length2<-bedpe$end2-bedpe$start2
  
  #chromosome preparation
  df1=as.data.frame(cbind(paste("chr",c(1:22,"X","Y"),sep=""),rep(1,24),
                          c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
                            135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,
                            48129895,51304566,155270560,59373566)))
  colnames(df1)<-c("Chrom1","start","end")
  df1$Chrom1<-factor(df1$Chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  df1$end<-as.numeric(as.character(df1$end))
  
  chra = chromosome
  sub<-bedpe %>% dplyr::filter(chrom1==chra & chrom2==chra) %>% dplyr::arrange(start1)
  # sub<-bedpe %>% dplyr::filter(start1 %in% c(142685331,142690334))
  pbpe = plotBedpe(sub,chrom = chra,chromstart = max(1,start),
                   chromend=min(dplyr::filter(df1,Chrom1==chra)$end,end),
                   heights = sub$score,bty='n',plottype=plot_type,offset=1,flip=F,
                   colorby=as.numeric(as.factor(sub$type)),colorbycol=hcl.colors,
                   border="black",lwd = 5)
  labelgenome(chra, 1, dplyr::filter(df1,Chrom1==chra)$end,side=1,scipen=20,n=30,scale="Mb",line=.18,chromline=.5,scaleline=0.5)
  legend("topright",inset =0.01,legend=levels(as.factor(sub$type)),col=hcl.colors(length(levels(as.factor(sub$type)))),pch=19,bty='n',text.font=2)
  # legend("topright",inset =0.01,legend=levels(as.factor(sub$samples)),col=f("Set1")[1:length(levels(as.factor(sub$samples)))],pch=19,bty='n',text.font=2)
  axis(side=2,las=2,tcl=.2)
  
}

intra_check("chr22",15000000,30000000)

#3. inter-chromosome activities

inter_check<-function(events){
  
  #read-in bedpe data
  bedpe<-read.delim("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",header = F)
  colnames(bedpe)<-c("chrom1","start1","end1","chrom2","start2","end2",
                     "type","score","number1","number2","number3",
                     "samples","Description")
  bedpe[,c(2,3,5,6,8,9,10,11)]<-sapply(bedpe[,c(2,3,5,6,8,9,10,11)],function(x){as.numeric(as.character(x))})
  bedpe$chrom1<-factor(bedpe$chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$chrom2<-factor(bedpe$chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$length1<-bedpe$end1-bedpe$start1
  bedpe$length2<-bedpe$end2-bedpe$start2
  
  sub<-bedpe %>% dplyr::filter(type %in% events)
  links_chromosomes_1 = gsub("chr","",sub$chrom1)
  links_chromosomes_2 = gsub("chr","",sub$chrom2)
  links_pos_1 = sub$start1
  links_pos_2 = sub$start2
  links_labels = NA
  
  tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 1,
                                       borderSize = 0, fillColors = "#EEFFEE")  
  
  tracklist = tracklist + BioCircosLinkTrack('myLinkTrack', links_chromosomes_1, links_pos_1,
                                             links_pos_1 , links_chromosomes_2, 
                                             links_pos_2, links_pos_2 ,
                                             maxRadius = 1, labels = links_labels)
  
  BioCircos(tracklist, genomeFillColor = "Paired",
            chrPad = 0.02, displayGenomeBorder = F, yChr =  T,
            genomeTicksDisplay = FALSE,  
            genomeLabelTextSize = "10pt", genomeLabelDy = 0,zoom = T,
            TEXTModuleDragEvent = TRUE)
}

inter_check("inverted-reciprocal")

#4. number of samples for inter-chromosome activities

inter_count<-function(events,resolution=1){
  
  #read-in bedpe data
  bedpe<-read.delim("C:/PhD/Rotations/Rotation_1/data/SV2/bedpe/combine/jointdf.bedpe",header = F)
  colnames(bedpe)<-c("Chrom1","start1","end1","Chrom2","start2","end2",
                     "type","score","number1","number2","number3",
                     "samples","Description")
  bedpe[,c(2,3,5,6,8,9,10,11)]<-sapply(bedpe[,c(2,3,5,6,8,9,10,11)],function(x){as.numeric(as.character(x))})
  bedpe$Chrom1<-factor(bedpe$Chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$Chrom2<-factor(bedpe$Chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  bedpe$length1<-bedpe$end1-bedpe$start1
  bedpe$length2<-bedpe$end2-bedpe$start2
  
  bedpe2<-bedpe
  bedpe2[,c(2,5)]<-sapply(bedpe2[,c(2,5)],function(x){round(x/resolution)})
  sv_hotspot1<-bedpe2 %>% group_by(Chrom1,start1,Chrom2,start2,type) %>% 
    summarise(n_sample = n_distinct(samples)) %>% ungroup() %>% 
    dplyr::arrange(Chrom1,start1)
  a<-sv_hotspot1 %>% dplyr::filter(type %in% events)
  a<-a %>% group_by(Chrom1,Chrom2) %>% summarise(sum_sample=sum(n_sample))
  a$Chrom1<-factor(a$Chrom1,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  a$Chrom2<-factor(a$Chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
  a1<-a %>% spread(Chrom2,sum_sample)
  a2<-as.matrix(a1[,2:ncol(a1)])
  rownames(a2)<-a1$Chrom1
  length(union(rownames(a2),colnames(a2))) #17 chromosomes involved
  allchr = union(rownames(a2),colnames(a2))
  orderlist = c(paste("chr",sort(as.numeric(gsub("chr","",allchr[which(!allchr %in% c("chrX","chrY"))]))),sep=""),
                allchr[which(allchr %in% c("chrX","chrY"))])
  
  f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
  grid.col = c(f("Set1"),f("Paired")[1:(length(orderlist)-length(f("Set1")))]) #extract 9 colors for each strain 
  names(grid.col)=orderlist
  chordDiagram(a2,directional = TRUE,order = orderlist,transparency = 0.6, grid.col = grid.col)
  
}

inter_count("inverted-reciprocal")








