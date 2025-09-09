#To generate the accounts table we use "featureCounts"

BiocManager::install("Rsubread")
library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(sva)
library(pamr)
library(Biobase)
library(GenomicRanges)
library(Rsubread)
library(DESeq2)
?featureCounts

# Cuntificación usando featureCounts
control_1mMPi_R1<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_41.sam",
                                annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                                isGTFAnnotationFile = TRUE,
                                GTF.featureType = "mRNA",
                                GTF.attrType = "ID",
                                countMultiMappingReads = TRUE,
                                fraction = TRUE,
                                allowMultiOverlap = TRUE,
                                isPairedEnd=TRUE)


control_1mMPi_R2<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_42.sam",
                                annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                                isGTFAnnotationFile = TRUE,
                                GTF.featureType = "mRNA",
                                GTF.attrType = "ID",
                                countMultiMappingReads = TRUE,
                                fraction = TRUE,
                                allowMultiOverlap = TRUE,
                                isPairedEnd=TRUE)


control_1mMPi_R3<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_43.sam",
                                annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                                isGTFAnnotationFile = TRUE,
                                GTF.featureType = "mRNA",
                                GTF.attrType = "ID",
                                countMultiMappingReads = TRUE,
                                fraction = TRUE,
                                allowMultiOverlap = TRUE,
                                isPairedEnd=TRUE)


T01mMPhi_R1<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/E250050987_clean_L01_13.bam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)




T01mMPhi_R2<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_45.sam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)

T01mMPhi_R3<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_46.sam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)

T05mMPhi_R1<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_47.sam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)

T05mMPhi_R2<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_48.sam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)


T05mMPhi_R3<-featureCounts("/Volumes/ADATA_cam/Transcriptomas fosfito/archivos_sam/V350255651_clean_L02_57.sam",
                           annot.ext = "/Volumes/ADATA_cam/Genoma Trichoderma atroviride version 3 final/Functional_annotation/Tatroviride_V3.gff3",
                           isGTFAnnotationFile = TRUE,
                           GTF.featureType = "mRNA",
                           GTF.attrType = "ID",
                           countMultiMappingReads = TRUE,
                           fraction = TRUE,
                           allowMultiOverlap = TRUE,
                           isPairedEnd=TRUE)

#We create the accounts table
matriz_cuentas<-cbind(control_1mMPi_R1$counts,control_1mMPi_R2$counts,control_1mMPi_R3$counts,
                      T01mMPhi_R1$counts,T01mMPhi_R2$counts,T01mMPhi_R3$counts,
                      T05mMPhi_R1$counts,T05mMPhi_R2$counts,T05mMPhi_R3$counts)
#We download the table of accounts in .csv and .txt format
write.csv(matriz_cuentas, file = "/Volumes/ADATA_cam/Transcriptomas fosfito/matriz_cuentas.csv")
write.table(matriz_cuentas, file = "/Volumes/ADATA_cam/Transcriptomas fosfito/matriz_cuentas_nuevo.txt")

#we create an output directory
outpathcount2 = "/Volumes/ADATA_cam/Transcriptomas fosfito/"
dir.create(outpathcount2, showWarnings=FALSE)

#we create the input directory
setwd("/Volumes/ADATA_cam/Transcriptomas fosfito/")

#Install the packages we are going to use
BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("tidyverse")
BiocManager::install("ggplot2")
BiocManager::install("ggfortify")
#Let's call the libraries of each of our packages
library(ggfortify)
library (ggplot2)
library(rafalib)
library(limma)
library(edgeR)


#let's read the table
counts = read.table("matriz_cuentas_nuevo.txt") 


#This command will eliminate those genes that have
#less than 2 reads and the sum of all the rows gives less than 5 reads.


counts = counts[rowSums(cpm(counts) >= 2) >=5,]

dim(counts


#This is to declare the groups we use in the DGELis object


grp = c("control_1mMPi","control_1mMPi","control_1mMPi",
        "T01mMPhi","T01mMPhi","T01mMPhi",
        "T05mMPhi","T05mMPhi","T05mMPhi")
#This command is to make the comparison matrix for the contrasts 
#taking into account the groups that we just formed
dge = DGEList(counts=counts, group=grp)
dge
#This will show us the size of the libraries and samples
dge$samples$lib.size
dge$samples
#We can graph the size of our libraries
barplot(dge$samples$lib.size, names=colnames(dge),las=2, col= c(2,2,2,3,3,3,4,4,4))

library(ggfortify)
pca_res <- prcomp(dge, scale. = TRUE)
autoplot(pca_res)

mds <- plotMDS(dge,
               xlim = c(-2, 2),
               xlab = "Dim1",
               ylim = c(-2, 2),
               ylab = "Dim2",
               cex.lab = 1.5,
               cex = 1.5,
               col = c(rep("red", 3), rep("black", 3), rep("blue", 3)),
               pch = 16,
               labels = NULL)



#Let's convert our numbers to log2 to graph the normalization.
logcounts= cpm(dge,log=TRUE)
#This command will graph the normalization of the libraries using boxplots.
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2)

dge = calcNormFactors(dge)
#Let's plot again to observe the changes after normalizing by MDS and boxplot
mds = plotMDS(dge,xlim=c(-2,2),xlab = "Dim1", ylim=c(-2,2), 
              ylab = "Dim2",cex.lab= 1.5, cex =1.5, #pch = 19,
              col=c(rep("red",3), rep("black",3),rep("blue",3),
                    cex.lab = 1, cex=1.5 ))#pch = 16, cex = 3))

logcounts= cpm(dge,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2)
#We can view the normalization values ​​that the program gives us.
dge$samples

#Performing differential expression analysis
#We will estimate the dispersion of our libraries so that we can use the value 
#in the differential expression analysis.
dge = estimateCommonDisp(dge)
dge$common.dispersion
dge = estimateGLMTagwiseDisp(dge)
dge$tagwise.dispersion

#With this command, we will graph the dispersion of the total libraries, 
#after the GMLTagwise analysis.
plotBCV(dge)
abline(h=c(-0.3,0.3), col="red")

dge$samples
levels(dge$samples$group)
#we compare different groups

design = model.matrix(~0+group, data=dge$samples) 
colnames(design) = levels(dge$samples$group)
design
#We make the comparisons and DEGs
fit = glmQLFit(dge, design)
my.contrasts = makeContrasts(T01mMPhivscontrol_1mMPi=T01mMPhi-control_1mMPi, T05mMPhivscontrol_1mMPi=T05mMPhi-control_1mMPi,  levels=design)


qlf.T01mMPhivscontrol_1mMPi = glmQLFTest(fit, contrast=my.contrasts[,"T01mMPhivscontrol_1mMPi"])
topTags(qlf.T01mMPhivscontrol_1mMPi)

qlf.T05mMPhivscontrol_1mMPi = glmQLFTest(fit, contrast=my.contrasts[,"T05mMPhivscontrol_1mMPi"])
topTags(qlf.T05mMPhivscontrol_1mMPi)

deTabT01mMPhivscontrol_1mMPi = deTab = topTags(qlf.T01mMPhivscontrol_1mMPi, n=Inf)$table
deTabT01mMPhivscontrol_1mMPi[c(2),]
deGenesT01mMPhi.control_1mMPi = rownames(deTabT01mMPhivscontrol_1mMPi)[deTabT01mMPhivscontrol_1mMPi$PValue < 0.05 & abs(deTabT01mMPhivscontrol_1mMPi$logFC) > 1 & abs(deTabT01mMPhivscontrol_1mMPi$FDR < 0.05)]
length(deGenesT01mMPhi.control_1mMPi)


plotSmear(dge[,c(4,5,6,1,2,3)], de.tags=deGenesT01mMPhi.control_1mMPi, cex=0.5)

plotSmear(dge[, c(4, 5, 6, 1, 2, 3)], 
          de.tags = deGenesT01mMPhi.control_1mMPi, 
          cex = 0.5,
          main = "Smear plot: T01mMPhi vs Control",
          xlab = "Average log CPM",
          ylab = "log-Fold Change")


write.table(deTabT01mMPhivscontrol_1mMPi[deGenesT01mMPhi.control_1mMPi,], file=paste(outpathcount2, "T01mMPhivscontrol_1mMPi_nuevo.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.csv(matriz_cuentas, file = "/Volumes/ADATA_cam/Transcriptomas fosfito/matriz_cuentas.csv")


deTabT05mMPhivscontrol_1mMPi = deTab = topTags(qlf.T05mMPhivscontrol_1mMPi, n=Inf)$table
deTabT05mMPhivscontrol_1mMPi[c(2),]
deGenesT05mMPhi.control_1mMPi = rownames(deTabT05mMPhivscontrol_1mMPi)[deTabT05mMPhivscontrol_1mMPi$PValue < 0.05 & abs(deTabT05mMPhivscontrol_1mMPi$logFC) > 1 & abs(deTabT05mMPhivscontrol_1mMPi$FDR < 0.05)]
length(deGenesT05mMPhi.control_1mMPi)
plotSmear(dge[,c(7,8,9,1,2,3)], de.tags=deGenesT05mMPhi.control_1mMPi, cex=0.7)
write.table(deTabT05mMPhivscontrol_1mMPi[deGenesT05mMPhi.control_1mMPi,], file=paste(outpathcount2, "T05mMPhivscontrol_1mMPi_nuevo.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")



