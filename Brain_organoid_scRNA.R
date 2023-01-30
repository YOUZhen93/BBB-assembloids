### 
11/14/2022 7:21:03 PM


## Doublet filtration
Using_Scrublet.py

generate double_res.tsv for each batch








# workdir
/public/home/gonglianggroup/zheny/temp/Brain_Organoid


options(stringsAsFactors=F)

library(Seurat)
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(harmony)


## step 1.

filelist
[1] "filtered_feature_bc_matrix_CCM_One"  
[2] "filtered_feature_bc_matrix_CCM_Sec"  
[3] "filtered_feature_bc_matrix_CCM_Third"
[4] "filtered_feature_bc_matrix_UHC_One"  
[5] "filtered_feature_bc_matrix_UHC_Sec"  
[6] "filtered_feature_bc_matrix_UHC_Third"


creatSeuO <- function(i) {
	pbmc = Read10X(i, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
	seurat_object = CreateSeuratObject(counts = pbmc, project = "organoid")
	doublet = read.table(paste0(i, "/double_res.tsv"), header = T, sep = "\t")
	seurat_object$doublet_scores = doublet$doublet_scores0.06
	seurat_object$predicted_doublets = doublet[,7]
	seurat_object = subset(seurat_object, subset = predicted_doublets == "False")
	ID = gsub("filtered_feature_bc_matrix_", "", i); seurat_object$batch = ID
	seurat_object$group = unlist(strsplit(ID, "_"))[1]
	print(paste0("Cell remaining after doublets remove: ", length(colnames(seurat_object))))
	saveRDS(seurat_object, file = paste0("./", i, "/raw_SeuObj.RDS"))
}

sapply(filelist, creatSeuO)


## Step 2.

setwd("/public/home/gonglianggroup/zheny/temp/Brain_Organoid")
filelist = c('filtered_feature_bc_matrix_CCM_One','filtered_feature_bc_matrix_CCM_Sec','filtered_feature_bc_matrix_CCM_Third',
	'filtered_feature_bc_matrix_UHC_One','filtered_feature_bc_matrix_UHC_Sec','filtered_feature_bc_matrix_UHC_Third')

firstMove <- function(i) {
	print(paste("Loading Seurat Object from", i, "......"))
	pbmc = readRDS(paste0(i, "/raw_SeuObj.RDS"))
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)

	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

	pbmc <- SCTransform(pbmc, verbose = FALSE)
	pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)

	pbmc <- RunPCA(pbmc, verbose = FALSE)
	pdf("elbow.pdf")
	ElbowPlot(pbmc)
	dev.off()
	print("saving Seurat Object .....")
	saveRDS(pbmc, file = paste0("./", i, "/PCA_SeuObj.RDS"))

}

sapply(filelist, firstMove)


## Step 3.

CCM group:
CCM = grep("CCM", filelist, value=T)

ONE = readRDS("./filtered_feature_bc_matrix_CCM_One/PCA_SeuObj.RDS")
SEC = readRDS("./filtered_feature_bc_matrix_CCM_Sec/PCA_SeuObj.RDS")
THI = readRDS("./filtered_feature_bc_matrix_CCM_Third/PCA_SeuObj.RDS")

pbmc.big <- merge(ONE, y = c(SEC, THI), add.cell.ids = c("One", "Sec", "Third"), project = "CCM")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc.big <- SCTransform(pbmc.big, verbose = FALSE)
pbmc.big <- CellCycleScoring(pbmc.big, s.features = s.genes, g2m.features = g2m.genes)

pbmc.big <- RunPCA(pbmc.big, verbose = FALSE)

pbmc.big = RunHarmony(pbmc.big, group.by.vars = "batch")
pdf("CCM_merged_elbow.pdf")
ElbowPlot(pbmc.big)
dev.off()

pbmc.big <- RunUMAP(pbmc.big, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:15, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE, resolution = 1)
markers = FindAllMarkers(pbmc.big, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
saveRDS(pbmc.big, file = "CCM_merged_SeuObj.RDS")
saveRDS(markers, file = "CCM_merged_raw_markers.RDS")




UHC group:
UHC = grep("UHC", filelist, value=T)

ONE = readRDS("./filtered_feature_bc_matrix_UHC_One/PCA_SeuObj.RDS")
SEC = readRDS("./filtered_feature_bc_matrix_UHC_Sec/PCA_SeuObj.RDS")
THI = readRDS("./filtered_feature_bc_matrix_UHC_Third/PCA_SeuObj.RDS")

pbmc.big <- merge(ONE, y = c(SEC, THI), add.cell.ids = c("One", "Sec", "Third"), project = "UHC")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc.big <- SCTransform(pbmc.big, verbose = FALSE)
pbmc.big <- CellCycleScoring(pbmc.big, s.features = s.genes, g2m.features = g2m.genes)

pbmc.big <- RunPCA(pbmc.big, verbose = FALSE)

pbmc.big = RunHarmony(pbmc.big, group.by.vars = "batch")
pdf("UHC_merged_elbow.pdf")
ElbowPlot(pbmc.big)
dev.off()

pbmc.big <- RunUMAP(pbmc.big, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:15, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE, resolution = 1)
markers = FindAllMarkers(pbmc.big, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
saveRDS(pbmc.big, file = "UHC_merged_SeuObj.RDS")
saveRDS(markers, file = "UHC_merged_raw_markers.RDS")






## Step 4.
# Futher QC

pbmc = readRDS("CCM_merged_SeuObj.RDS")
qc1 = as.data.frame(pbmc@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nCount_RNA), list(mean=mean, median=median)))
qc1 = qc1[order(qc1$mean), ]; thres = quantile(pbmc$nCount_RNA, 0.1)
delc1 = qc1$seurat_clusters[which(qc1$mean < thres)]
pbmc = subset(pbmc, idents = as.character(delc1), invert = T)
## QC genes expressed

'''
qc2 = as.data.frame(pbmc@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nFeature_RNA), list(mean=mean, median=median)))
qc2 = qc2[order(qc2$mean), ]; thres = quantile(pbmc$nFeature_RNA, 0.1)
delc2 = qc2$seurat_clusters[which(qc2$mean < thres)]
'''

uhc = readRDS("UHC_merged_SeuObj.RDS")
qc1 = as.data.frame(uhc@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nCount_RNA), list(mean=mean, median=median)))
qc1 = qc1[order(qc1$mean), ]; thres = quantile(uhc$nCount_RNA, 0.1)
delc1 = qc1$seurat_clusters[which(qc1$mean < thres)]
uhc = subset(uhc, idents = as.character(4), invert = T)




pbmc.big <- merge(pbmc, y = uhc, add.cell.ids = c("CCM", "UHC"), project = "Organoid")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc.big <- SCTransform(pbmc.big, verbose = FALSE)
pbmc.big <- CellCycleScoring(pbmc.big, s.features = s.genes, g2m.features = g2m.genes)

pbmc.big <- RunPCA(pbmc.big, verbose = FALSE)

pbmc.big = RunHarmony(pbmc.big, group.by.vars = "batch")
pdf("TwoGrps_merged_elbow.pdf")
ElbowPlot(pbmc.big)
dev.off()

pbmc.big <- RunUMAP(pbmc.big, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:15, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE, resolution = 1)
markers = FindAllMarkers(pbmc.big, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
saveRDS(pbmc.big, file = "TwoGrps_merged_SeuObj.RDS")
saveRDS(markers, file = "TwoGrps_merged_raw_markers.RDS")





pbmc = readRDS("TwoGrps_merged_SeuObj.RDS")
ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
pdf("TwoGrps_merged_DimPlot.pdf", w=10, h=8)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = colorRampPalette(brewer.pal(8, "Set1"))(6), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("red", "navy"), group.by = 'group') + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()

markers = readRDS("TwoGrps_merged_raw_markers.RDS")







ct = c(
	6,11,
13, 5,
21,
0,
17, 20)

ctn = c(
	rep("Proliferative Astrocytes", 2), 
	rep("Astrocytes 1", 2), ## FABP7 PTN PTPRZ1 SLC1A3
	"Astrocytes 2",  # VEGFA BNIP3 SLC16A3 SLC1A3 IGFBP2 PGK1 FAM162   cell migration/protein binding
	"Astrocytes 3",  # CXCL14
	rep("G2M phase Astrocytes", 2)) # CCNB2 BIRC5 CENPF SLC1A3 
	





ct = c(
	22, 
	26,
	23,
	3,
	2,
	7,18,
	15)




ctn = c(
	rep("GABAergic neurons 1",1), # GAD1 GAD2 DLX2 DLX1
	rep("Proliferative GABAergic neurons 2",1), # GAD1 GAD2 DLX2 DLX1
	rep('Glutamatergic neurons 1', 1), # GABBR2  GRIN2B BHLHE22
	rep('Glutamatergic neurons 2', 1), #  GRIN2B NEUROD6 NEUROD2 BHLHE22 Precursors of Glut
	rep('Glutamatergic neurons 3', 1), # NEUROD6 NEUROD2 BHLHE22 SLA Precursors of Glut
	rep('Glutamatergic neurons 4', 2), # NHLH1 EOMES IGFBPL1  UBC_like
	rep('Glutamatergic neurons 5', 1)) # DDIT4 PGK1 STMN2 IGFBP2 MALAT1  cell migration/protein binding




ct = c(
	19, 24,
	8,
	1,
	10,
	12,
	4,
	9)
	


ctn = c(
	rep("Brain vascular endothelials",2), ## SLC2A1 ENO1, PGK1, IGFBP2   cell migration/protein binding
	"Proliferative cells", # TOP2A MKI67 CENPF
	rep("Vascular leptomeningeal cells 1(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM  TAC1
	rep("Vascular leptomeningeal cells 2(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM TSHZ2 IGF2 TGFBI
	rep("Tendon cells", 1), # THBS4  SCXA
	rep('Pericytes 1', 1), # PDGFRB RGS5
	rep('Pericytes 2', 1)) # PDGFRB BST2 MEF2C REN


	
ct = c(
	14,
	25,
	27,
	16)
	


ctn = c(
	"Vascular endothelials 1", # CD34 CD93 CLDN5
	"Vascular endothelials 2", # KDR CLDN5 IL32 CD93
	"Vascular endothelials 3",  # PRND CLDN5  ESM1
	"Smooth muscle cells (SMCs)" # ACTA2 TAGLN
	)) 






ct = c(
	6,11,
13, 5,
21,
0,
17, 20,
22, 
	26,
	23,
	3,
	2,
	7,18,
	15,
	19, 24,
	8,
	1,
	10,
	12,
	4,
	9,
	14,
	25,
	27,
	16)

ctn = c(
	rep("Proliferative Astrocytes", 2), # SLC1A3 HIST1H4C
	rep("Astrocytes 1", 2), ## FABP7 PTN PTPRZ1 SLC1A3
	"Astrocytes 2",  # VEGFA BNIP3 SLC16A3 SLC1A3 IGFBP2 PGK1 FAM162   cell migration/protein binding
	"Astrocytes 3",  # CXCL14
	rep("G2M phase Astrocytes", 2),
	rep("GABAergic neurons 1",1), # GAD1 GAD2 DLX2 DLX1
	rep("Proliferative GABAergic neurons 2",1), # GAD1 GAD2 DLX2 DLX1
	rep('Glutamatergic neurons 1', 1), # GABBR2  GRIN2B BHLHE22
	rep('Glutamatergic neurons 2', 1), #  GRIN2B NEUROD6 NEUROD2 BHLHE22 Precursors of Glut
	rep('Glutamatergic neurons 3', 1), # NEUROD6 NEUROD2 BHLHE22 SLA Precursors of Glut
	rep('Glutamatergic neurons 4', 2), # NHLH1 EOMES IGFBPL1  UBC_like
	rep('Glutamatergic neurons 5', 1),
	rep("Brain vascular endothelials",2), ## SLC2A1 ENO1, PGK1, IGFBP2   cell migration/protein binding
	"Proliferative cells", # TOP2A MKI67 CENPF
	rep("Vascular leptomeningeal cells 1(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM  TAC1
	rep("Vascular leptomeningeal cells 2(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM TSHZ2 IGF2 TGFBI
	rep("Tendon cells", 1), # THBS4  SCXA
	rep('Pericytes 1', 1), # PDGFRB RGS5
	rep('Pericytes 2', 1),
	"Vascular endothelials 1", # CD34 CD93 CLDN5
	"Vascular endothelials 2", # KDR CLDN5 IL32 CD93
	"Vascular endothelials 3",  # PRND CLDN5  ESM1
	"Smooth muscle cells (SMCs)") # CCNB2 BIRC5 CENPF SLC1A3 
	



names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc))]
names(cellT) = names(Idents(pbmc))
Idents(pbmc) = factor(cellT)
pbmc$Sub_celltype = as.character(Idents(pbmc))
pbmc$Global_celltype = gsub("\\ \\d", "", pbmc$Sub_celltype)



pbmc$Sub_celltype = factor(pbmc$Sub_celltype, levels = c("Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "GABAergic neurons 1", 
	'Glutamatergic neurons 1', 'Glutamatergic neurons 2', 'Glutamatergic neurons 3', 'Glutamatergic neurons 4', 'Glutamatergic neurons 5',
	"Vascular leptomeningeal cells 1(VLMCs)", "Vascular leptomeningeal cells 2(VLMCs)", 'Pericytes 1', 'Pericytes 2', "Tendon cells",
	"Brain vascular endothelials", "Vascular endothelials 1", "Vascular endothelials 2", "Vascular endothelials 3",
	"Smooth muscle cells (SMCs)", "Proliferative Astrocytes", "G2M phase Astrocytes", "Proliferative GABAergic neurons 2", "Proliferative cells"))

pbmc$Global_celltype = factor(pbmc$Global_celltype, levels = c("Astrocytes", "GABAergic neurons", 'Glutamatergic neurons',
	"Vascular leptomeningeal cells(VLMCs)", 'Pericytes', "Tendon cells", "Brain vascular endothelials", "Vascular endothelials",
	"Smooth muscle cells (SMCs)", "Proliferative Astrocytes", "G2M phase Astrocytes", "Proliferative GABAergic neurons", 
	"Proliferative cells"))

shortNa = c("Ast 1", "Ast 2", "Ast 3", "GABAs 1", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5", "VLMCs 1", "VLMCs 2", "Peri 1", 
	"Peri 2", "TCs", "Brain VEs", "VEs 1", "VEs 2", "VEs 3", "SMCs", "Pro Asts", "G2M Asts", "Pro GABAs", "PCs")
names(shortNa) = levels(pbmc$Sub_celltype)

pbmc$short_SubCT = shortNa[match(pbmc$Sub_celltype, names(shortNa))]
pbmc$short_SubCT = factor(pbmc$short_SubCT, levels=shortNa)
identical(as.numeric(table(pbmc$short_SubCT)), as.numeric(table(pbmc$Sub_celltype)))

shortNa = c("Ast", "GABAs", "GluNs", "VLMCs", "Peri", "TCs", "Brain VEs", "VEs", "SMCs", "Pro Asts", "G2M Asts", "Pro GABAs", "PCs")
names(shortNa) = levels(pbmc$Global_celltype)
pbmc$short_GlobalCT = shortNa[match(pbmc$Global_celltype, names(shortNa))]
pbmc$short_GlobalCT = factor(pbmc$short_GlobalCT, levels=shortNa)
identical(as.numeric(table(pbmc$short_GlobalCT)), as.numeric(table(pbmc$Global_celltype)))



ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
pdf("TwoGrps_merged_DimPlot_withCT_Sub.pdf", w=17, h=8)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T, group.by = 'Sub_celltype') +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()

pdf("TwoGrps_merged_DimPlot_withCT_Glo.pdf", w=13, h=8)
DimPlot(pbmc, label = TRUE, cols = colorRampPalette(brewer.pal(8, "Set1"))(13), label.size = 6, repel = T, group.by = 'Global_celltype') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()


pdf("Feature_CellCycleScore.pdf", w = 16, h=7)
FeaturePlot(pbmc, features = c("S.Score", "G2M.Score"), pt.size = 0, cols = c("grey", "#db8495", "#b31031"), ncol = 2, raster = TRUE)
dev.off()


saveRDS(pbmc, file = "TwoGrps_merged_SeuObj.RDS")


## some stats
CCM_One   CCM_Sec CCM_Third   UHC_One   UHC_Sec UHC_Third 
     7403      8352      6690      8191      6640      9017 


CCM   UHC 
22445 23848 




## cell composition in two groups
df = pbmc@meta.data
df = df[,c('short_SubCT', 'group')]
plotdf = as.data.frame(df %>% count(short_SubCT, group, name="Count"))
plotdf = rbind(plotdf, data.frame(short_SubCT="VEs 3", group="CCM", Count=0))
plotdf$Count = plotdf$Count/sum(plotdf$Count)
plotdf$Count = round(plotdf$Count * 100, digit = 3)
#plotdf = plotdf %>% group_by(Shortname) %>% summarise(Count = sum(Count))
#plotdf = as.data.frame(plotdf[order(plotdf$Count, decreasing = T), ])
pdf("cell_composition_Sub.pdf", w = 10, h = 6)
ggplot(plotdf, aes(x = short_SubCT, y = Count, fill = group)) + geom_bar(stat="identity", width = 0.5, position=position_dodge()) + scale_fill_manual(values = c("red","navy")) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black'), axis.text.x=element_text(angle=90, hjust=1, vjust=1)) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))
dev.off()


# Global
df = pbmc@meta.data
df = df[,c('short_GlobalCT', 'group')]
plotdf = as.data.frame(df %>% count(short_GlobalCT, group, name="Count"))
plotdf$Count = plotdf$Count/sum(plotdf$Count)
plotdf$Count = round(plotdf$Count * 100, digit = 3)
#plotdf = plotdf %>% group_by(short_GlobalCT) %>% summarise(Count = sum(Count))
#plotdf = as.data.frame(plotdf[order(plotdf$Count, decreasing = T), ])
pdf("cell_composition_Global.pdf", w = 10, h = 6)
ggplot(plotdf, aes(x = short_GlobalCT, y = Count, fill = group)) + geom_bar(stat="identity", width = 0.5, position=position_dodge()) + scale_fill_manual(values = c("red","navy")) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black'), axis.text.x=element_text(angle=90, hjust=1, vjust=1)) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))
dev.off()





## Step 5.

pbmc = readRDS("TwoGrps_merged_SeuObj.RDS")
Idents(pbmc) = pbmc$short_SubCT
markers = FindAllMarkers(pbmc, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
saveRDS(markers, file = "TwoGrps_merged_markers_SubCT.RDS")

mm2 = c("FABP7","PTN","VEGFA","SLC1A3","BNIP3","CXCL14","GAD1","GAD2","DLX2","GABBR2",
	"GRIN2B","NEUROD6","SLA","STMN2","EOMES","PDGFRB","LUM","COL3A1","CLDN5",
	"TAC1","C7","ACTA2","THBS4","RGS5","IGF1","CD34","MECOM", "SLC2A1","KDR","HIST1H4C","IGFBP2","PGK1")	

markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
top10 <- newmarkers %>% arrange(p_val_adj) %>%
    group_by(cluster) %>% 
    slice(1:5)


top10$cluster = factor(top10$cluster, levels = levels(Idents(pbmc)))
top10 = top10[order(top10$cluster), ];
pdf("Heatmap.pdf", w = 12, h = 9)
DoHeatmap(pbmc, features = top10$gene, raster=T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))
dev.off()



-------------------- make stackVlnplot

modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-5, 0, -5, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-3, 0, -3, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, color="black"), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mymarker = c('SLC1A3', 'HES5', 'VEGFA', 'BNIP3', 'CXCL14', 'GAD2', 'DLX1', 'GABBR2', 'GRIN2B', 'BHLHE22', 'NEUROD2', 
	'EOMES', 'DDIT4', 'LUM', 'COL3A1', 'TSHZ2', 'PDGFRB', 'RGS5', 'REN', 'THBS4', 'SLC2A1', 'CLDN5', 'CD34', 'KDR', 'PRND',
	 'ACTA2', 'TAGLN', 'HIST1H1B', 'PTTG1', 'TOP2A')
pdf("stackedVlnPlot.pdf", w= 15, h=40)
StackedVlnPlot(obj = pbmc, features = mymarker, cols = mycolors)
dev.off()

