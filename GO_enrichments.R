#Enrichments in GO terms
BiocManager::install("clusterProfiler")
BiocManager::install("topGO")
library("clusterProfiler") 
library("topGO")
library("dplyr")

setwd("/Volumes/ADATA_cam/Transcriptomas fosfito/genes_unicos _compartidos/")


outpathcount = ("/Volumes/ADATA_cam/Transcriptomas fosfito/genes_unicos _compartidos/")

Id_Go<- readMappings ("Tatroviride_GO_V3.txt",sep="\t", IDsep = ";")

##Biological Process##

GoF_BP<-annFUN.gene2GO(whichOnto = "BP", feasibleGenes = NULL, gene2GO = Id_Go)


termGo_BP<-data.frame()

for(i in names(GoF_BP)){              
  y=GoF_BP[[i]]                       
  ln=length(GoF_BP[[i]])              
  x=data.frame(GO=rep(i,ln),y)        
  termGo_BP=rbind(termGo_BP,x)       
}

tGo_BP<-buildGOmap(termGo_BP)
term_BP<-tGo_BP$GO
addterm_BP<-go2term(term_BP)
colnames(addterm_BP)<- c("GO", "Term")

Gomap_BP<-left_join(tGo_BP, addterm_BP, by = "GO") 

term2gene_bp<-Gomap_BP[,c("GO", "Gene")]
term2name_bp<-Gomap_BP[,c("GO", "Term")]  

##Celullar Component##

GoF_CC<-annFUN.gene2GO(whichOnto = "CC", feasibleGenes = NULL, gene2GO = Id_Go)


termGo_CC<-data.frame()

for(i1 in names(GoF_CC)){              
  y1=GoF_CC[[i1]]                       
  ln1=length(GoF_CC[[i1]])              
  x1=data.frame(GO=rep(i1,ln1),y1)        
  termGo_CC=rbind(termGo_CC,x1)        
}


tGo_CC<-buildGOmap(termGo_CC)
term_CC<-tGo_CC$GO
addterm_CC<-go2term(term_CC)
colnames(addterm_CC)<- c("GO", "Term")
#b<-unique(tGo_CC$GO) #Corroborar cu?ntos GOs hay
Gomap_CC<-left_join(tGo_CC, addterm_CC, by = "GO") 

term2gene_cc<-Gomap_CC[,c("GO", "Gene")]
term2name_cc<-Gomap_CC[,c("GO", "Term")]

## Molecular Function ##

GoF_MF<-annFUN.gene2GO(whichOnto = "MF", feasibleGenes = NULL, gene2GO = Id_Go)


termGo_MF<-data.frame()

for(i2 in names(GoF_MF)){              
  y2=GoF_MF[[i2]]                       
  ln2=length(GoF_MF[[i2]])              
  x2=data.frame(GO=rep(i2,ln2),y2)        
  termGo_MF=rbind(termGo_MF,x2)        
}


tGo_MF<-buildGOmap(termGo_MF)
term_MF<-tGo_MF$GO
addterm_MF<-go2term(term_MF)
colnames(addterm_MF)<- c("GO", "Term")
#b<-unique(tGo_MF$GO) #Corroborar cu?ntos GOs hay
Gomap_MF<-left_join(tGo_MF, addterm_MF, by = "GO") 
term2gene_mf<-Gomap_MF[,c("GO", "Gene")]
term2name_mf<-Gomap_MF[,c("GO", "Term")]

#########################
list.files()
arch = list.files(pattern="go")
arch

list_modulo <- list()


for (i in arch) {
  name = sub("", "", i)
  print(name)
  data = read.table(i, header = TRUE, sep="\t",  quote= "", row.names=1, comment.char="")
  list_modulo[[name]] = data
}

list_enrch_BP<-list()
list_enrch_CC<-list()
list_enrch_MF<-list()
list_FC<-list()



for (i in arch){
  names_file<-sub("","",i)
  print(names_file)
  data.dge<-list_modulo[[names_file]]
  #print(head(data.dge))
  dge_fc<-data.dge$net.colors
  names(dge_fc)<-rownames(data.dge)
  dge_fc<-sort(dge_fc, decreasing = TRUE)
  list_FC[[names_file]]<-dge_fc
  dge_id<-rownames(data.dge)
  dge.erch_BP<-enricher(dge_id, TERM2GENE=term2gene_bp, TERM2NAME=term2name_bp, pAdjustMethod = "fdr", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1)
  dge.erch_CC<-enricher(dge_id, TERM2GENE=term2gene_cc, TERM2NAME=term2name_cc, pAdjustMethod = "fdr", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1)
  dge.erch_MF<-enricher(dge_id, TERM2GENE=term2gene_mf, TERM2NAME=term2name_mf, pAdjustMethod = "fdr", pvalueCutoff = 0.5, qvalueCutoff = 1, minGSSize = 1)
  list_enrch_BP[[names_file]]<-dge.erch_BP
  list_enrch_CC[[names_file]]<-dge.erch_CC
  list_enrch_MF[[names_file]]<-dge.erch_MF
  write.table(dge.erch_BP, file=paste(outpathcount,"BP_",names_file, sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
  write.table(dge.erch_CC, file=paste(outpathcount,"CC_",names_file, sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
  write.table(dge.erch_MF, file=paste(outpathcount,"MF_",names_file, sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
}

