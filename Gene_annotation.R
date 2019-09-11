#Gene_annotation_Pathway

library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)
library(KEGGREST)
library(org.Hs.eg.db)

# color extraction function
f <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)

########################
#   Read in bed file   #
########################

bedfile<-read.csv("C:/PhD/Rotations/Rotation_1/writing_valor2/table/inv_del_gene.csv")
ethdf<-bedfile %>% group_by(samples) %>% summarise(ethnicity=unique(ethnicity)) %>% ungroup()

#########################################
#   Convert ID for pathway enrichment   #
#########################################

listMarts()
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
searchFilters(mart=ensembl,pattern="ensembl.*id")
searchFilters(mart=ensembl,pattern="symbol")
idtable<-getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id_version",
                              "hgnc_symbol","uniprot_gn_symbol"),mart=ensembl)
colnames(idtable)<-c("ensembl_gene_id","transcription","gene_symbol","uniprot_symbol")
bedfile<-bedfile %>% left_join(idtable,by="transcription")
allgenes<-unique(idtable$gene_symbol)[which(unique(idtable$gene_symbol)!="") & which(!is.na(unique(idtable$gene_symbol)))]

##########################
#   Pathway enrichment   #
##########################

#a) retrieve all pathways

keggdf<-as.data.frame(keggList("pathway","hsa"))
keggdf$pathwayidx<-rownames(keggdf)
colnames(keggdf)<-c("pathway","idx")

#b) get all pathway gene lists

kegglist<-list()
for (i in 1:nrow(keggdf)){
  print(i)
  nm<-gsub("path:","",keggdf$idx[i])
  gnlist<-keggGet(nm)[[1]]$GENE
  if(length(gnlist)>=2){
    # oddvals <- seq(1, length(gnlist), by=2) # extract only odd indexed gene names (only "Rvxxxx")
    evenvals <- seq(2, length(gnlist), by=2)
  }else{
    gns=gnlist
  }
  gns<-as.character(sapply(gnlist[evenvals],function(x){strsplit(x,";")[[1]][1]}))
  kegglist[[i]]<-gns
}
names(kegglist)<-keggdf$idx
saveRDS(kegglist,"C:/PhD/Rotations/KEGG_HSA_genelist.rds")

#c) hypergeometric test

hypergeofun<-function(sample_name){
  
  sub<-bedfile %>% dplyr::filter(samples==sample_name)
  genes<-unique(sub$gene_symbol)
  genes<-genes[which(genes!="")]
  genes<-genes[!is.na(genes)]
  
  plist = lapply(kegglist,function(k){
    gnset = k
    p=phyper(q=length(unique(intersect(gnset,genes)))-1,
             m=length(unique(gnset)),
             n=length(unique(allgenes[which(!allgenes %in% unique(gnset))])),
             k=length(unique(genes)),lower.tail = F, log.p = F)
  })
  pdf = as.data.frame(do.call("rbind",plist))
  pdf$pathway<-rownames(pdf)
  colnames(pdf)[1]<-"pvalue"
  pdf$padj<-p.adjust(pdf$pvalue,method="fdr")
  pdf<-pdf %>% dplyr::filter(padj<=0.1)
  return(pdf)
  
}
sample_enrichment<-lapply(unique(bedfile$samples),hypergeofun)
names(sample_enrichment)<-as.character(unique(bedfile$samples))

sublist<-sample_enrichment[names(which(sapply(sample_enrichment,nrow)>0))]
sigdf<-as.data.frame(do.call("rbind",sublist))
sigdf$samples<-sapply(rownames(sigdf),function(k){strsplit(k,".",fixed = TRUE)[[1]][1]})
sigdf<-sigdf %>% left_join(ethdf,by="samples")
table(sigdf$ethnicity)
table(sigdf$pathway)
sigsamples<-unique(sigdf$samples)

bedfile %>% dplyr::filter(samples %in% sigsamples) #all hitting UGT2B4























