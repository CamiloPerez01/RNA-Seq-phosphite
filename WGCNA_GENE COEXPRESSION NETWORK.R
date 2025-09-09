#Gene coexpression network (WGCNA)
setwd("/Volumes/ADATA_cam/redes_fosfito/")
outpathcount2 = ("/Volumes/ADATA_cam/redes_fosfito/")

red = read.table("matrix_accounts.txt")

library("limma")
library("edgeR")
library("statmod")
library("DESeq2")

red


counts_red = red[rowSums(cpm(red) >= 10) >=2,]
bibliosmico = as.matrix(counts_red)

bibliosmico

bibliosmicolnc <- round(bibliosmico)

vsd <- vst(bibliosmicolnc, blind=FALSE)
vsd

BiocManager :: install ("WGCNA")
library(WGCNA)

options(stringsAsFactors = FALSE)

dim(vsd)
names(vsd)

WTvsPatol_cpm <- (vsd)
head(WTvsPatol_cpm)

datExpr= t(WTvsPatol_cpm[,])
names(datExpr)= row.names(WTvsPatol_cpm)
rownames(datExpr)=names(WTvsPatol_cpm)
dim(datExpr)
head(datExpr)

gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust
sampleTree = hclust(dist(datExpr), method = "average")

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

abline(h = 90, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 80, minSize = 1)
table(clust)

#Clust 1 contains the samples we want to keep -------
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2, main = "Sample dendrogram", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

powers = c(c(1:10), seq(from = 11, to=20, by=1))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

cor <- WGCNA::cor



#-----------------------------------------

net = blockwiseModules(datExpr, 
                       #similarity matrix options 
                       corType = "pearson",
                       #opciones de matriz de adyacencia
                       power = 14,
                       networkType = "signed",
                       #opciones TOM
                       TOMType = "signed",
                       #Opciones de identificación de modulos
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "red_fosfito",
                       verbose = 3)
table(net$colors)
head(net)


sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
mergedColors

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


mergedColors = labels2colors(net$colors)
mergedColors
table(mergedColors)

ls (net)

library(factoextra)

install.packages("dplyr")
library(dplyr)

x<-as.data.frame(net$colors)
x

x0<-x%>%
  filter(net$colors=="0")

x1<-x%>%
  filter(net$colors=="1")

x2<-x%>%
  filter(net$colors=="2")

x3<-x%>%
  filter(net$colors=="3")

x4<-x%>%
  filter(net$colors=="4")

x5<-x%>%
  filter(net$colors=="5")

x6<-x%>%
  filter(net$colors=="6")

x7<-x%>%
  filter(net$colors=="7")

x8<-x%>%
  filter(net$colors=="8")

x9<-x%>%
  filter(net$colors=="9")
x10<-x%>%
  filter(net$colors=="10")
x11<-x%>%
  filter(net$colors=="11")

x12<-x%>%
  filter(net$colors=="12")
x13<-x%>%
  filter(net$colors=="13")
x14<-x%>%
  filter(net$colors=="14")
x15<-x%>%
  filter(net$colors=="15")
x16<-x%>%
  filter(net$colors=="16")
x17<-x%>%
  filter(net$colors=="17")

write.table(x, file=paste(outpathcount, "modulos_totales_de_la_red_power18_pearson_lncRNAs.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x0, file=paste(outpathcount, "module_0_grey_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x1, file=paste(outpathcount, "module_1_turquoise_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x2, file=paste(outpathcount, "module_2_blue_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x3, file=paste(outpathcount, "module_3_brown_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x4, file=paste(outpathcount, "module_4_yellow_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x5, file=paste(outpathcount, "module_5_green_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x6, file=paste(outpathcount, "module_6_red_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x7, file=paste(outpathcount, "module_7_black_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x8, file=paste(outpathcount, "module_8_pink_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x9, file=paste(outpathcount, "module_9_magenta_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")


write.table(x10, file=paste(outpathcount, "module_10_purple_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x11, file=paste(outpathcount, "module_11_greenyellow_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x12, file=paste(outpathcount, "module_12_tan_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x13, file=paste(outpathcount, "module_13_salmon_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x14, file=paste(outpathcount, "module_14_cyan_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x15, file=paste(outpathcount, "module_15_midnightblue_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x16, file=paste(outpathcount, "module_16_lightcyan_pearson_power14.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")
write.table(x17, file=paste(outpathcount, "module_17_grey60_pearson_power18.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")

MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")



dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 14);


MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes


TOM = TOMsimilarityFromExpr(datExpr, power = 14);

probes = colnames(datExpr)
names(datExpr)
probes

cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput-edges0.35_red_pow14_pearson", ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes0.35_red_pow14_pearson", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.35,
                               nodeNames = probes,
                               nodeAttr = moduleColors);


colors <- net$colors

topHubs <- chooseTopHubInEachModule(datExpr, 
                                    colors,
                                    omitColors = "grey",
                                    power = 14,
                                    type = "signed")

gen_hub
print(topHubs)

