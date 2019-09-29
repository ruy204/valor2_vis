#Hotspot_visualization

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
source("C:/PhD/Rotations/Rotation_1/scripts/Visualization_functions.R")

##################
#   Data Input   #
##################

#1. read reciprocal hotspot (from locations)

fromdf<-read.csv("C:/PhD/Rotations/Rotation_1/plots/valor2/all/tables/reciprocal_inv_from_hotspot.csv")
todf<-read.csv("C:/PhD/Rotations/Rotation_1/plots/valor2/all/tables/reciprocal_inv_to_hotspot.csv")

#####################
#   Visualization   #
#####################

#1. visualize reciprocal hotspot (from locations)

ethrank = c('Finnish','SOUTHERN HAN CHINESE','PUERTO RICAN','UTAH/MORMON','YORUBA/Nigeria','HAN CHINESE/China',
            'JAPANESE/Japan','USA/MEXICAN','USA/AFRICAN','ITALY/TOSCANI','Caucasian')
fromdf$Description<-factor(fromdf$Description,levels=ethrank)
fromdf$Chrom2<-factor(fromdf$Chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))

#chromosome preparation
df1=as.data.frame(cbind(paste("chr",c(1:22,"X","Y"),sep=""),rep(1,24),
                        c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,
                          135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,
                          48129895,51304566,155270560,59373566)))
colnames(df1)<-c("Chrom2","start","end")
df1$Chrom2<-factor(df1$Chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
df1$end<-as.numeric(as.character(df1$end))

ggplot(data=df1,aes(Chrom2,end))+geom_bar(stat="identity",width = 0.2,alpha=0.5)+ylim(0,max(df1$end))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(data=fromdf,aes(Chrom2,start2,color=type),alpha=0.5,size=3)+facet_wrap(~Description,nrow=2)+
  scale_color_brewer(palette = "Set1")+labs(title="Target locations of (inv-)reciprocal hotspot on Chromosome 22 location 44Mb")

#2. visualize reciprocal hotspot (to locations)

todf$Description<-factor(todf$Description,levels=ethrank)
todf$Chrom2<-factor(todf$Chrom2,levels=paste("chr",c(1:22,"X","Y"),sep=""))
ggplot(data=df1,aes(Chrom2,end))+geom_bar(stat="identity",width = 0.2,alpha=0.5)+ylim(0,max(df1$end))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_point(data=todf,aes(Chrom1,start1,color=type),alpha=0.5,size=3)+facet_wrap(~Description,nrow=2)+
  scale_color_brewer(palette = "Set1")+labs(title="Original locations of (inv-)reciprocal hotspot on Chromosome 22 location 44Mb")

#############
#   Merge   #
#############

todf2<-todf %>% group_by(type,Chrom1,from_location) %>% mutate(n_event=n()) %>% ungroup()

ggplot(data=df1,aes(Chrom2,end))+geom_bar(stat="identity",width = 0.2,alpha=0.5)+ylim(0,max(df1$end))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_jitter(data=todf2,aes(Chrom1,start1,color=n_event),alpha=0.8,size=3)+facet_wrap(~type,nrow=1)+
  labs(title="Original locations of (inv-)reciprocal hotspot on Chromosome 22 location 44Mb")+
  scale_color_gradient(low="pink",high="darkblue")
















