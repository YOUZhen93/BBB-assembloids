library(Seurat)
library(harmony)
library(data.table)
library(Matrix)

options(stringsAsFactors=F)





 
#####   ALS brain Organoid scRNA-Seq 20220208  resolution 0.5
ct = c(
	7,
	9,
	0,
	3,
	2,
	11,
	6,
	1,
	4,
	5,
	10,
	8
	)



ctn = c(
	"proliferative  Progenitors", ## MKI67, TOP2A, CENPF
	"proliferative  Progenitors", ## MKI67, NUSAP1, CCNB1
	"Neural progenitors 1",  ## SOX2, NES, HES1 CD9 CALB1 ZIC1
	"Neural progenitors 2",  ## SOX2, NES, HES1 PTPRZ1 SLC1A3 FABP7
	"GABAergic neurons",  ## SLC17A7, TBR1, SLC17A6, GAD2
	"Glutamatergic neurons",  ## GABBR2, GRIN2B, SLC17A7, SLC17A6, GABBR1
	"Proliferative cells", ## UBE2C TOP2A CENPF MKI67
	rep("Vascular leptomeningeal cells (VLMCs)", 1), ## LUM, COL3A1, COL1A1
	"Vascular endothelials 1", ## IGFBP3 ALDOA, ENO1
	"Perivascular adipocyte", ## IGF1, MEOX2, MECOM
	"Tendon cells",  ##  TNMD, THBS4
	"Vascular endothelials 2" ## VEGFA, ALDOA, ENO1
	)

names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc))]
names(cellT) = names(Idents(pbmc))
Idents(pbmc) = factor(cellT)
pbmc$celltype = as.character(Idents(pbmc))


#####   ALS brain Organoid scRNA-Seq 20220208  resolution 1
# 3/9/2022 6:03:08 PM
ct = c(
	14,8,
	2,6,16,
	4,11,
	9,
	7, 
	17,
	15,
	18, 13,
	1, 5,
	0,
	3,
	10,
	12
	)



ctn = c(
	rep("Proliferative  astrocytes", 2), ## MKI67, TOP2A, CENPF
	rep("Astrocytes", 3),  ## SLC1A3 and NOTCH1
	rep("Neural progenitors", 2),  ## SOX2, NES, HES1 CD9 CALB1 ZIC1
	"GABAergic neurons",  ## GAD2, GAD1, DLX2 DLX5 DLX1 
	"Glutamatergic neurons 1",  ## GRIN2B, SLC17A6, GABBR2, TBR1, STMN2, DCX
	"Glutamatergic neurons 2",  ## SLC17A6, TBR1, NEUROD6, IGFBPL1
	"Glutamatergic neurons 3",  ## GABBR2, SLC17A6, GRIN2B, SLC17A7, TBR1, LMO3, ARPP21
	rep("Proliferative cells", 2), ## UBE2C TOP2A CENPF MKI67
	rep("Vascular leptomeningeal cells (VLMCs)", 2), ## LUM, COL3A1, COL1A1, IGF2, PLAT
	"Brain vascular endothelials", ## SLC2A1
	"Vascular endothelials 2", ## VEGFA, ALDOA, ENO1
	"Perivascular adipocyte", ## IGF1, MEOX2, MECOM
	"Tendon cells"  ##  TNMD, THBS4, SCXA
	)

names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc))]
names(cellT) = names(Idents(pbmc))
Idents(pbmc) = factor(cellT)
pbmc$celltype = as.character(Idents(pbmc))


library(tidyverse)

markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10




Idents(pbmc) = factor(Idents(pbmc), levels = c("Neural progenitors", "Astrocytes", "Proliferative astrocytes", 'GABAergic neurons', 
	'Glutamatergic neurons 1', 'Glutamatergic neurons 2', 'Glutamatergic neurons 3', 'Vascular leptomeningeal cells (VLMCs)', 'Brain vascular endothelials', 'Vascular endothelials 2', 
	'Perivascular adipocyte','Tendon cells', 'Proliferative cells'))

# short celltype
Idents(pbmc) = factor(Idents(pbmc), levels = c("NPs", "Ast", "PAst", 'GABAs', 
	'GluNs1', 'GluNs2', 'GluNs3', 'VLMCs', 'Brain VEs', 'VEs2', 'PAs','TCs', 'PCs'))


top10$cluster = factor(top10$cluster, levels = c("Neural progenitors", "Astrocytes", "Proliferative astrocytes", 'GABAergic neurons', 
	'Glutamatergic neurons 1', 'Glutamatergic neurons 2', 'Glutamatergic neurons 3', 'Vascular leptomeningeal cells (VLMCs)', 'Brain vascular endothelials', 'Vascular endothelials 2', 
	'Perivascular adipocyte','Tendon cells', 'Proliferative cells'))
top10 = top10[order(top10$cluster), ]

DoHeatmap(pbmc, features = top10$gene) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))


plot1 = FeaturePlot(pbmc, features = c("GAD1", "GFAP", "S100B", "ALDH1L1", "PDGFRB", "ACTA2", "SLC2A1", "CLDN5", "TJP1", "OCLN", "AIF1"),cols = c("grey", "red"), ncol = 4, pt.size = 1, raster = T, combine=F) 
plot1 <- lapply(plot1, FUN = function(x) x + theme(axis.text = element_text(size=16,color='black'), axis.title = element_text(size=18), plot.title = element_text(size = 18,face='bold')))
CombinePlots(plots = plot1)


ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T, pt.size = 1) + NoLegend() + 
    theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))





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

mymarker = c('SLC1A3', 'NOTCH1', 'SOX2', 'NES', 'TOP2A', 'CENPF','MKI67', 'GAD1', 'GAD2', 'DLX2', 'SLC17A6', 'SLC17A7', 'TBR1', 'GABBR2', 'GRIN2B', 'GABBR1', 
	'LUM', 'COL3A1', 'COL1A1', 'ALDOA', 'SLC2A1', 'FN1', 'RSPO3', 'IGF1', 'MEOX2', 'MECOM', 'TNMD', 'THBS4')
StackedVlnPlot(obj = pbmc, features = mymarker, cols = mycolors)





#### Second batch of brain organoids scRNA-Seq
#### 3/11/2022 9:29:09 AM

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(DoubletFinder)
options(stringsAsFactors = F)

setwd('F:/Brain_organoid_scRNA/Sec-filtered_feature_bc_matrix/filtered_feature_bc_matrix')
x = Read10X(".", unique.features = TRUE, strip.suffix = FALSE)
pbmc <- CreateSeuratObject(counts = x, project = "brain_organoids", min.cells = 3, min.features = 200)
pbmc$batch = "batch2"

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)

# rerun mark
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc)

# pbmc = RunHarmony(pbmc, group.by.vars = "samples")
pbmc <- RunUMAP(pbmc, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 1)

ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ncols)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))

## QC 
## QC UMI
qc1 = as.data.frame(pbmc@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nCount_RNA), list(mean=mean, median=median)))
qc1 = qc1[order(qc1$mean), ]; thres = quantile(pbmc$nCount_RNA, 0.1)
delc1 = qc1$seurat_clusters[which(qc1$mean < thres)]
## QC genes expressed
qc2 = as.data.frame(pbmc@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nFeature_RNA), list(mean=mean, median=median)))
qc2 = qc2[order(qc2$mean), ]; thres = quantile(pbmc$nFeature_RNA, 0.1)
delc2 = qc2$seurat_clusters[which(qc2$mean < thres)]

## remove cluster 14
pbmc = pbmc[,which(pbmc$seurat_clusters != "14")]
# rerun from rerun mark
pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
ElbowPlot(pbmc)
pbmc <- RunUMAP(pbmc, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)

### DoubletFinder part cannot be run in my own PC
sweep.res.pbmc <- paramSweep_v3(pbmc, PCs = 1:15, sct = TRUE) # optimize params
sweep.stats_kidney <- summarizeSweep(sweep.res.pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_kidney)
barplot(bcmvn_pbmc$BCmetric, names.arg = bcmvn_pbmc$pK, las=2)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.075*nrow(pbmc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
# ********************************************************


pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 1)
ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ncols)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))

markers = FindAllMarkers(pbmc, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]



### batch remove
pbmc1 = readRDS("../../filtered_feature_bc_matrix/brain_organoids_SeuratObj.RDS")
pbmc$batch = "batch2"
pbmc1$batch = "batch1"

pbmc.combined <- merge(pbmc1, y = pbmc, add.cell.ids = c("batch1", "batch2"), project = "brain_organoids")

pbmc.combined <- SCTransform(pbmc.combined, verbose = FALSE)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

pbmc.combined = RunHarmony(pbmc.combined, group.by.vars = "batch")
ElbowPlot(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:20, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:20, verbose = FALSE)
pbmc.combined <- FindClusters(pbmc.combined, verbose = FALSE, resolution = 1)
markers = FindAllMarkers(pbmc.combined, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
save(pbmc.combined, markers, file = "combined_organoids_res.Rdata")






# 3/13/2022 1:58:46 PM

load("combined_organoids_res.Rdata")

DimPlot(pbmc.combined, label = TRUE, cols = c("red", "navy"), label.size = 6, repel = T, group.by = "batch") + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))

ncols = length(unique(Idents(pbmc.combined)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.combined, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))


# QC
qc1 = as.data.frame(pbmc.combined@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nCount_RNA), list(mean=mean, median=median)))
qc1 = qc1[order(qc1$mean), ]; thres = quantile(pbmc.combined$nCount_RNA, 0.1)
delc1 = qc1$seurat_clusters[which(qc1$mean < thres)]
## QC genes expressed
qc2 = as.data.frame(pbmc.combined@meta.data %>% group_by(seurat_clusters) %>% summarise_at(vars(nFeature_RNA), list(mean=mean, median=median)))
qc2 = qc2[order(qc2$mean), ]; thres = quantile(pbmc.combined$nFeature_RNA, 0.1)
delc2 = qc2$seurat_clusters[which(qc2$mean < thres)]



ct = c(
	9,
	11,
	18,
	7,
	2,5,
	21,
	23,
	15,
	19,
	0,
	13,
	16,
	1,10,
	3,
	17,
	6,
	8,12,
	4,14,
	22,
	20
	)



ctn = c(
	"Proliferative  cells", ## MKI67, TOP2A, CENPF
	"Proliferative  astrocytes", ## MKI67, TOP2A, CENPF
	"Proliferative  astrocytes", ## MKI67, TOP2A, CENPF
	"Brain vascular endothelials", ## SLC2A1
	rep("Astrocytes 1", 2), ## FABP7 PTN PTPRZ1 SLC1A3 NOTCH1
	"Neural progenitors 2", # PAX3 CNPY1 EN2
	"Endothelial cells", # CLDN5 ESAM CDH5
	"Astrocytes 2", 	# CALB1 LPL PLP1 	NTRK2
	"Astrocytes 3",  # VEGFA BNIP3 SLC16A3 SLC1A3
	"Neural progenitors 1", # C1QTNF4 CD9 HES1 SOX2
	"Astrocytes 4", # CXCL14 NOG SLC1A3
	"GABAergic neurons", # GAD1 GAD2 DLX2 DLX1
	rep('Glutamatergic neurons 1', 2), # SLC17A6 GABBR2 GRIN2B STMN2 DCX
	'Glutamatergic neurons 2', # NEUROD6 LMO3 BHLHE22 NEUROD2   Precursors of Glut
 	'Glutamatergic neurons 3', # NHLH1 EOMES IGFBPL1
 	"Tendon cells",  ##  TNMD, THBS4, SCXA
 	rep("Perivascular adipocyte", 2), ## IGF1, MEOX2, MECOM
 	rep("Vascular leptomeningeal cells 1(VLMCs)", 2), ## LUM, COL3A1, COL1A1, IGF2, PLAT
 	"Vascular leptomeningeal cells 2(VLMCs)", # GUCY1A3 C7 CELF3 
 	"Vascular endothelials 3" ## EIF1, ALDOA, ENO1 FTL
)



names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc.combined))]
names(cellT) = names(Idents(pbmc.combined))
Idents(pbmc.combined) = factor(cellT)
pbmc.combined$celltype = as.character(Idents(pbmc.combined))
pbmc.combined$celltype2 = pbmc.combined$celltype
pbmc.combined$celltype2 = gsub("\\ \\d", "", pbmc.combined$celltype)
shortN = unique(pbmc.combined$celltype);
names(shortN) = c("PAs", "Brain VEs", "Ast 2", "VLMCs 1", "Ast 1", "PAst", "GluNs 1", "TCs", "VEs 3", "GluNs 2", "GABAs", "PCs", "Ast 4", "GluNs 3", "NPs 1", 
	"Ast 3", "VLMCs 2", "NPs 2", "ECs")
pbmc.combined$short_celltype = names(shortN)[match(pbmc.combined$celltype, shortN)]

ncols = length(unique(Idents(pbmc.combined)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.combined, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))



markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10




#Idents(pbmc) = factor(Idents(pbmc), levels = c("Neural progenitors", "Astrocytes", "Proliferative astrocytes", 'GABAergic neurons', 
#	'Glutamatergic neurons 1', 'Glutamatergic neurons 2', 'Glutamatergic neurons 3', 'Vascular leptomeningeal cells (VLMCs)', 'Brain vascular endothelials', 'Vascular endothelials 2', 
#	'Perivascular adipocyte','Tendon cells', 'Proliferative cells'))

# short celltype
Idents(pbmc.combined) = factor(Idents(pbmc.combined), levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))


top10$cluster = factor(top10$cluster, levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))
top10 = top10[order(top10$cluster), ]

DoHeatmap(pbmc.combined, features = top10$gene, raster = T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))


#plot1 = FeaturePlot(pbmc, features = c("GAD1", "GFAP", "S100B", "ALDH1L1", "PDGFRB", "ACTA2", "SLC2A1", "CLDN5", "TJP1", "OCLN", "AIF1"),cols = c("grey", "red"), ncol = 4, pt.size = 1, raster = T, combine=F) 
#plot1 <- lapply(plot1, FUN = function(x) x + theme(axis.text = element_text(size=16,color='black'), axis.title = element_text(size=18), plot.title = element_text(size = 18,face='bold')))
#CombinePlots(plots = plot1)
saveRDS(pbmc.combined, file = "combined_organoids_celltype.Rdata")




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

ncols = length(unique(Idents(pbmc.combined)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)

mymarker = c('SLC1A3', 'NOTCH1', 'SOX2', 'NES', 'TOP2A', 'CENPF','MKI67', 'GAD1', 'GAD2', 'DLX2', 'SLC17A6', 'SLC17A7', 'TBR1', 'GABBR2', 'GRIN2B', 'GABBR1', 
	'LUM', 'COL3A1', 'COL1A1', 'ALDOA', 'SLC2A1', 'FN1', 'RSPO3', 'IGF1', 'MEOX2', 'MECOM', 'TNMD', 'THBS4', 'CLDN5', 'ESAM')
StackedVlnPlot(obj = pbmc.combined, features = mymarker[1:15], cols = mycolors)
StackedVlnPlot(obj = pbmc.combined, features = mymarker[16:30], cols = mycolors)


### cell type proportions
df = pbmc.combined@meta.data
df = df[,c('short_celltype', 'batch')]
plotdf = as.data.frame(df %>% count(short_celltype, batch, name="Count"))

plotdf$Count[which(plotdf$batch == "batch1")] = plotdf$Count[which(plotdf$batch == "batch1")]/9342
plotdf$Count[which(plotdf$batch == "batch2")] = plotdf$Count[which(plotdf$batch == "batch2")]/8550
plotdf$Count = round(plotdf$Count * 100, digit = 3)

plotdf$short_celltype = factor(plotdf$short_celltype, levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))

ggplot(plotdf, aes(x = batch, y = Count, fill = short_celltype)) + geom_bar(stat="identity", width = 0.5) + scale_fill_manual(values = mycolors) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))



plotdf = as.data.frame(df %>% count(short_celltype, batch, name="Count"))
plotdf$Count = plotdf$Count/17892
plotdf$Count = round(plotdf$Count * 100, digit = 3)
plotdf1 = plotdf %>% group_by(short_celltype) %>% summarise(sum(Count))
plotdf1 = plotdf1[order(plotdf1$`sum(Count)`, decreasing = T), ]
plotdf$short_celltype = factor(plotdf$short_celltype, levels = as.character(plotdf1$short_celltype))
ggplot(plotdf, aes(x = short_celltype, y = Count, fill = batch)) + geom_bar(stat="identity", width = 0.5) + scale_fill_manual(values = c('red', 'navy')) + 
    theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
                            axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black'), 
                            axis.text.x = element_text(size=18, color="black", angle=45, hjust=1, vjust=1)) + xlab("") + ylab("Percentage (%)") + 
    scale_y_continuous(expand=c(0,0))


## level 2
ncols = 11
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)

shortN = unique(pbmc.combined$celltype2)
names(shortN) = c("PAs", "Brain VEs", "Ast", "VLMCs", "PAst", "GluNs", "TCs", "VEs", "GABAs", "PCs", "NPs", "ECs")
pbmc.combined$short_celltype2 = names(shortN)[match(pbmc.combined$celltype2, shortN)]
df = pbmc.combined@meta.data
df = df[,c('short_celltype2', 'batch')]
plotdf = as.data.frame(df %>% count(short_celltype2, batch, name="Count"))

plotdf$Count = plotdf$Count/17892
plotdf$Count = round(plotdf$Count * 100, digit = 3)
plotdf = plotdf %>% group_by(short_celltype2) %>% summarise(Count = sum(Count))
plotdf = as.data.frame(plotdf[order(plotdf$Count, decreasing = T), ])

plotdf$short_celltype2 = factor(plotdf$short_celltype2, levels = plotdf$short_celltype2)

ggplot(plotdf, aes(x = batch, y = Count, fill = short_celltype2)) + geom_bar(stat="identity", width = 0.5) + scale_fill_manual(values = mycolors) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))



plotdf1 = plotdf %>% group_by(short_celltype2) %>% summarise(sum(Count))
plotdf1 = plotdf1[order(plotdf1$`sum(Count)`, decreasing = T), ]
plotdf$short_celltype2 = factor(plotdf$short_celltype2, levels = as.character(plotdf1$short_celltype2))
ggplot(plotdf, aes(x = short_celltype2, y = Count, fill = batch)) + geom_bar(stat="identity", width = 0.5) + scale_fill_manual(values = c('red', 'navy')) + 
    theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
                            axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black'), 
                            axis.text.x = element_text(size=18, color="black", angle=45, hjust=1, vjust=1)) + xlab("") + ylab("Percentage (%)") + 
    scale_y_continuous(expand=c(0,0))

plotdf = as.data.frame(table(pbmc.combined$short_celltype))
plotdf$Freq = 100 * plotdf$Freq/sum(plotdf$Freq)
plotdf = plotdf[order(plotdf$Freq, decreasing = T), ]; plotdf$Freq = round(plotdf$Freq, digit = 2)
plotdf$Var1 = factor(plotdf$Var1, levels = plotdf$Var1)
ncols = length(unique(plotdf$Var1))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
ggplot(plotdf, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat="identity", width = 0.7) + geom_text(aes(label = Freq), size=6, color="black", vjust = 0) + 
scale_fill_manual(values = mycolors) + theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text.x = element_text(size=18, color="black", angle=45, hjust=1, vjust=1),
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0)) + theme(legend.position = 'none')

plotdf = as.data.frame(table(pbmc.combined$short_celltype2))
plotdf$Freq = 100 * plotdf$Freq/sum(plotdf$Freq)
plotdf = plotdf[order(plotdf$Freq, decreasing = T), ]; plotdf$Freq = round(plotdf$Freq, digit = 2)
plotdf$Var1 = factor(plotdf$Var1, levels = plotdf$Var1)
ncols = length(unique(plotdf$Var1))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
ggplot(plotdf, aes(x = Var1, y = Freq, fill = Var1)) + geom_bar(stat="identity", width = 0.7) + geom_text(aes(label = Freq), size=6, color="black", vjust = 0) + 
scale_fill_manual(values = mycolors) + theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text.x = element_text(size=18, color="black", angle=45, hjust=1, vjust=1),
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0)) + theme(legend.position = 'none')



### cell cell communication
### analysis
# 3/10/2022 11:41:59 AM
# walk through
require(pbmcapply)
pbmc = readRDS("combined_organoids_celltype.Rdata")
markers = read.table("celltype_enriched_markers.xls", header=T, sep="\t")
Idents(pbmc) = pbmc$short_celltype
Idents(pbmc) = factor(Idents(pbmc), levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))
shortN = unique(Idents(pbmc))
names(shortN) = unique(pbmc$short_celltype)
markers$cluster = names(shortN)[match(markers$cluster, shortN)]

# normalized
countTable = pbmc@assays$RNA@counts
LR = read.csv("C:/Users/youzh/Downloads/mouse_human_lrdb_version2.csv")
countTable = as.matrix(countTable[rownames(countTable) %in% unique(c(LR$Ligand.Human, LR$Receptor.Human)), ])
cTMM = TMMnormalization(countTable)
celltypes = levels(unique(Idents(pbmc)))


celltype1 = ""; celltype2 = ""

condition = as.character(Idents(pbmc));



#  make all combination df
workdf = expand.grid(celltypes, celltypes)
names(workdf) = c("Sender", "Receiver")

totaldf = data.frame()
CCTObj = list()
for (i in seq(nrow(workdf))) {
    print(paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> "))
    tempdf = workmyassoff(cTMM, workdf$Sender[i], workdf$Receiver[i], 0.95, condition, 1)
    tdf1 = tempdf@result_df; tdf1$Sender = workdf$Sender[i]; tdf1$Receiver = workdf$Receiver[i];
    totaldf = rbind(totaldf, tdf1)
    CCTObj[[paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> ")]] = c(tempdf@paramL, tempdf@paramR)
}


### outputs: CCTObj list params of all ligands and receptors
###          tempdf CellCrossTalker object
###          totaldf final KL divergences data frame of all celltypes

##  Check results
ranks = seq(200)
resL <- unlist(pbmcapply::pbmclapply(ranks, PMFofDGBD, 
    param=unlist(CCTObj[['NPs >>> GABAs']]["CTSD"]), mc.cores=1))
resR <- unlist(pbmcapply::pbmclapply(ranks, PMFofDGBD, 
    param=unlist(CCTObj[['NPs >>> GABAs']]["ATP6AP2"]), mc.cores=1))
plot(c(resL, resR), pch = 19, col = factor(rep(c('a', 'b'), each=200)))

resL <- unlist(pbmcapply::pbmclapply(ranks, PMFofDGBD, 
    param=unlist(CCTObj[['NPs >>> GABAs']]["HSP90AB1"]), mc.cores=1))
resR <- unlist(pbmcapply::pbmclapply(ranks, PMFofDGBD, 
    param=unlist(CCTObj[['NPs >>> GABAs']]["NRIP1"]), mc.cores=1))
plot(c(resL, resR), pch = 19, col = factor(rep(c('a', 'b'), each=200)))









library(tidyverse)
library(reshape2)
library(circlize)
library(ComplexHeatmap)

totaldf$cell_cell = paste(totaldf$Sender, totaldf$Receiver, sep = ">>>")
totaldf$LR = paste(totaldf$Ligand, totaldf$Receptor, sep = "---")
totaldf = arrange(totaldf, cell_cell, LRscore)


##  KL divergence 0 occur
if (min(totaldf$KLscore) == 0) {
	Mt = floor(log10(min(totaldf$KLscore[which(totaldf$KLscore != 0)])))
	Mt = Mt - 2
	totaldf$KLscore[which(totaldf$KLscore == 0)] = 10^Mt
}
totaldf$LRscore = -log10(totaldf$KLscore) ## The smaller the KL divergence the greater connections
totaldf$LRscore = totaldf$LRscore + abs(min(totaldf$LRscore)) + 0.1
totaldf = arrange(totaldf, cell_cell, desc(LRscore))

save(totaldf, CCTObj, file = "CCC_CellCrossTalker2.Rdata")





x = as.data.frame(totaldf %>% group_by(cell_cell) %>% top_n(5, wt = LRscore))
x$Sender = as.character(x$Sender); x$Receiver = as.character(x$Receiver)
x = x[-grep("PAst", x$cell_cell), ]; x = x[-grep("PCs", x$cell_cell), ]
x = x[-which(x$Sender %in% c("Brain VEs", "VEs 3", "VLMCs 1", "VLMCs 2", "TCs", "PAs") & x$Receiver %in% c("Brain VEs", "VEs 3", "VLMCs 1", "VLMCs 2", "TCs", "PAs")), ]

x = x[-which(x$Sender == x$Receiver), ]

x1 = x[which(x$Sender %in% c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3") & x$Receiver %in% c("VLMCs 1", "VLMCs 2", "Brain VEs", "TCs", "PAs", "VEs 3", "ECs")), ]
x1 = x[which(x$Sender %in% c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3") & x$Receiver %in% c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3")), ]
x1 = x[which(x$Sender %in% c("VLMCs 1", "VLMCs 2", "Brain VEs", "TCs", "PAs", "VEs 3", "ECs") & x$Receiver %in% c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3")), ]
x1 = x[which(x$Sender %in% c("VLMCs 1", "VLMCs 2", "Brain VEs", "TCs", "PAs", "VEs 3", "ECs") & x$Receiver %in% c("VLMCs 1", "VLMCs 2", "Brain VEs", "TCs", "PAs", "VEs 3", "ECs")), ]


ggplot(data = x1, mapping = aes_string(x = "cell_cell", 
                                             y = "LR")) + 
  geom_point(mapping = aes_string(size="LRscore",color="LRscore"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=median(x1$LRscore),limits=c(min(x1$LRscore),max(x1$LRscore)))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(range = c(2,10))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
# 11 X 8 

totaldf$LRscore = totaldf$LRscore + abs(min(totaldf$LRscore))
CCC_table = aggregate(totaldf$LRscore, by = list(totaldf$cell_cell), FUN=sum)
names(CCC_table) = c("cell_cell", "Aggr_score")
CCC_table = CCC_table[order(CCC_table$Aggr_score, decreasing=T), ]
# CCC_table$Aggr_score = CCC_table$Aggr_score + abs(min(CCC_table$Aggr_score)) + 1

x = CCC_table[-grep("PCs", CCC_table$cell_cell), ]
x = x[-grep("PAst", x$cell_cell), ]
x$Sender = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
x$Receiver = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
mtx1 = acast(x, Sender~Receiver, value.var="Aggr_score") # convert to nonsymmetric matrix
fords = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'ECs')
reindex = match(fords, colnames(mtx1))
mtx1 = mtx1[reindex, reindex]

col_fun = colorRamp2(c(min(mtx1)+median(mtx1)/sd(mtx1), median(mtx1), max(mtx1)-median(mtx1)/sd(mtx1)), c("#0000ff", "#fffafa", "#ff0000"))  
Heatmap(mtx1, name = "CCC",col=col_fun,column_names_side = "top", rect_gp = gpar(col = "black", lwd = 2),
            
            column_names_gp = gpar(fontsize = 12,fontface="bold"),
            column_names_rot = -45,row_names_side = "left", row_names_gp = gpar(fontsize = 12,fontface="bold"),
            
            cluster_rows = FALSE,cluster_columns = FALSE,   heatmap_legend_param = list(
              title = "Score",grid_width=unit(0.6,"cm")
            ))








###  Wed May 11 13:55:22 2022
###  Three batches combined

pbmc.combined <- SCTransform(object = pbmc.combined, verbose = FALSE)

NE = pbmc[,Idents(pbmc) %in% c(15, 3, 8, 10, 9, 7, 21)]
NE <- FindNeighbors(NE, dims = 1:15, verbose = FALSE)
NE <- FindClusters(NE, verbose = FALSE, resolution = 0.5)
markers_NE = FindAllMarkers(NE, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)


ct = c(
0,1,3,
9,
2,
4,
5,6,
7,
8
)

ctn = c(
	rep("Astrocytes 1", 3), ## FABP7 PTN PTPRZ1 SLC1A3
	"Astrocytes 3",  # VEGFA BNIP3 SLC16A3 SLC1A3
	"Astrocytes 2", 	# CALB1 LPL PLP1 	NTRK2
	"Astrocytes 4",  # CXCL14
	rep("Neural progenitors 1", 2), # SIX3 SIX6 ZIC4 
	"Neural progenitors 2", # EN2 PAX3 CNPY1
	"Proliferative  astrocytes") # KIAA0101 MCM4 CENPK

names(ctn) = ct
cellT = ctn[as.character(Idents(NE))]
names(cellT) = names(Idents(NE))
Idents(NE) = factor(cellT)
NE$celltype2 = as.character(Idents(NE))
NE$celltype1 = gsub("\\ \\d", "", NE$celltype2)




EC = pbmc[,Idents(pbmc) %in% c(7,15,13,21,0,1,17,14)]; 

EC <- FindNeighbors(EC, dims = 1:15, verbose = FALSE)
EC <- FindClusters(EC, verbose = FALSE, resolution = 0.8)
markers_EC = FindAllMarkers(EC, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers_EC = markers_EC[which(markers_EC$p_val_adj < 0.01), ]

ct = c(
	0,5,
	2,
	8,
	12,
	11,
	1,9,
	4,
	7,
	6,
	3,
	10
	)
	


ctn = c(
	rep('Pericytes', 2), # PDGFRB
	rep("Vascular leptomeningeal cells 1(VLMCs)", 1), ## LUM, COL3A1, COL1A1, PLAT
 	rep("Vascular leptomeningeal cells 2(VLMCs)", 1), # TAC1 C7  Fibromyocyte
 	rep("Vascular leptomeningeal cells 3(VLMCs)", 1), # PLAU COLEC10 FOXF1 
	"Smooth muscle cells (SMCs)", # ACTA2
 	rep("Brain vascular endothelials",2), ## SLC2A1
	"Vascular endothelials 1", # NPTX2 MECOM
 	rep("Tendon cells", ),  ##  TNMD, THBS4, SCXA
 	rep("Venous Endothelials 1", 2), # TSHZ2, CPE
 	"Venous Endothelials 2",   #AUTS2, LMO4 expressed in prefrontal, and frontal/occipital cortex respectively
 	"Venous Endothelials 3")   #ENO1, PGK1, IGFBP2 



names(ctn) = ct
cellT = ctn[as.character(Idents(EC))]
names(cellT) = names(Idents(EC))
Idents(EC) = factor(cellT)
EC$celltype2 = as.character(Idents(EC))
EC$celltype1 = gsub("\\ \\d", "", EC$celltype2)






GN = pbmc[,Idents(pbmc) %in% c(1,14,5,11,19)]; 

GN <- FindNeighbors(GN, dims = 1:15, verbose = FALSE)
GN <- FindClusters(GN, verbose = FALSE, resolution = 0.8)
markers_GN = FindAllMarkers(GN, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers_GN = markers_GN[which(markers_GN$p_val_adj < 0.01), ]

ct = c(
	5,
	8,
	0,7,10,
	3,4,
	2,9,6,
	1
)



ctn = c(
	rep("GABAergic neurons",1), # GAD1 GAD2 DLX2 DLX1
	rep('Glutamatergic neurons 1', 1), # GABBR2 GRIN2B  
	rep('Glutamatergic neurons 2', 3), # NEUROD6  BHLHE22 NEUROD2   Precursors of Glut
	rep('Glutamatergic neurons 3', 2), # NHLH1 EOMES IGFBPL1  UBC_like
	rep('Glutamatergic neurons 4', 3), # STMN2, STMN1, GAP43
	rep('Glutamatergic neurons 5', 1)) # EMP3, TRMT112, CDKN1A



names(ctn) = ct
cellT = ctn[as.character(Idents(GN))]
names(cellT) = names(Idents(GN))
Idents(GN) = factor(cellT)
GN$celltype2 = as.character(Idents(GN))
GN$celltype1 = gsub("\\ \\d", "", GN$celltype2)


meta = rbind(NE@meta.data[,c("celltype2","celltype1")], GN@meta.data[,c("celltype2","celltype1")], EC@meta.data[,c("celltype2","celltype1")])
pbmc$celltype1 = meta$celltype1[match(colnames(pbmc), rownames(meta))]
pbmc$celltype2 = meta$celltype2[match(colnames(pbmc), rownames(meta))]
pbmc$celltype1[which(pbmc$seurat_clusters %in% c(13,18))] <- pbmc$celltype2[which(pbmc$seurat_clusters %in% c(13,18))] <- 'Proliferative  astrocytes'
pbmc$celltype1[which(pbmc$seurat_clusters %in% c(12,20))] <- pbmc$celltype2[which(pbmc$seurat_clusters %in% c(12,20))] <- 'Proliferative  cells'








# remove 25 26 clusters
library(ggpubr)
UCs = c("Artery_UCell","Capillary_UCell","Venule_UCell","Venous_UCell","Pericyte_UCell","Smooth_Muscle_UCell",          
"Fibromyocyte_UCell","Perivascular_Fibroblast_UCell")
EC <- AddModuleScore_UCell(EC, features = signatures)
FeaturePlot(EC, "Venous_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Venule_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Pericyte_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Capillary_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Smooth_Muscle_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Fibromyocyte_UCell", pt.size = 0, cols = c("blue", "grey", "red"))
FeaturePlot(EC, "Perivascular_Fibroblast_UCell", pt.size = 0, cols = c("blue", "grey", "red"))

res = list()
for (i in c("Artery_UCell","Capillary_UCell","Venule_UCell","Venous_UCell","Pericyte_UCell","Smooth_Muscle_UCell",          
"Fibromyocyte_UCell","  _Fibroblast_UCell")) {
	df = EC@meta.data[,c("short_celltype", i)]; names(df) = c("groups", "values"); values = df$values; groups = factor(df$groups)
	q = pairwise.wilcox.test(x = values, g = groups, p.adjust.method = 'bonf')
	res[[i]] = q$p.value
}


library("gplots")
library(gridGraphics)
library(grid)
library(gridExtra)

grab_grob <- function(){
  grid.echo()
  grid.grab()
}


gl = lapply(seq(length(res)), function(i) {
	p = res[[i]]; p = ifelse(p > 0.01, 0, ifelse(p>0.001, 1, ifelse(p>0.0001, 2, ifelse(p>0.00001, 3, 4)))); p[upper.tri(p)] = NA;
	heatmap.2(p, scale = "none", col = c('white',brewer.pal(n = 7, name = "Reds")), trace = "none", density.info = "none", Rowv = FALSE, Colv = FALSE, 
		colsep = 1:11, rowsep = 1:11, main=names(res)[[i]])
	grab_grob()
	})

grid.newpage()
grid.arrange(grobs=gl, ncol=2, clip=TRUE)




# load Science brain EC Seurat Obj
table(pbmc$celltype2)

pbmc = pbmc[,pbmc$celltype2 == "FB"]
length(pbmc)
pbmc = pbmc[,colnames(pbmc) %in% sample(colnames(pbmc), 1000, replace=F)]
EC$batch2 = "original"; pbmc$batch2 = "science"
pbmc.combined <- merge(EC, y = pbmc, add.cell.ids = c("original", "science"), project = "brain_organoids")
pbmc.combined <- SCTransform(pbmc.combined, verbose = FALSE, conserve.memory = T)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

pbmc.combined = RunHarmony(pbmc.combined, group.by.vars = "batch2", assay.use = "SCT")
ElbowPlot(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:15, verbose = FALSE)
pbmc.combined <- FindClusters(pbmc.combined, verbose = FALSE, resolution = 1)

#ncols = length(unique(Idents(pbmc)))
#mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.combined, label = TRUE, group.by = "batch2", cols = c('blue','red'), label.size = 6, repel = T, pt.size = 1) + NoLegend() +
    theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))


pbmc = pbmc[,pbmc$celltype2 == "SMC"]
length(pbmc)
pbmc = pbmc[,colnames(pbmc) %in% sample(colnames(pbmc), 2000, replace=F)]
EC$batch2 = "original"; pbmc$batch2 = "science"
pbmc.combined <- merge(EC, y = pbmc, add.cell.ids = c("original", "science"), project = "brain_organoids")
pbmc.combined <- SCTransform(pbmc.combined, verbose = FALSE, conserve.memory = T)
pbmc.combined <- RunPCA(pbmc.combined, verbose = FALSE)

pbmc.combined = RunHarmony(pbmc.combined, group.by.vars = "batch2", assay.use = "SCT")
ElbowPlot(pbmc.combined)
pbmc.combined <- RunUMAP(pbmc.combined, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
pbmc.combined <- FindNeighbors(pbmc.combined, dims = 1:15, verbose = FALSE)
pbmc.combined <- FindClusters(pbmc.combined, verbose = FALSE, resolution = 1)

#ncols = length(unique(Idents(pbmc)))
#mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.combined, label = TRUE, group.by = "batch2", cols = c('blue','red'), label.size = 6, repel = T, pt.size = 1) + NoLegend() +
    theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))


patchwork::wrap_plots(plotlist = plotlist, ncol = 2)







homology genes of human mouse
homology_gene_MGI.trans.txt
http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt

homo = homo[,c(2,4)]; names(homo) = c("species", "genes")

lll = list()
for (i in seq(nrow(homo)-1)) {
	if (grepl("mouse", homo$species[i])) {
		while (homo$species[i+1] == "human") {
			lll[[homo$genes[i]]] = c(lll[[homo$genes[i]]], homo$genes[i+1])
			i = i+1
		}
	}
}

df = data.frame(mouse = rep(names(lll), sapply(lll, length)), human = unlist(lll))


Murine_Cell_EC.Robj


celltypes = table(murine_EC[,murine_EC$tissue == "Brain"]$ann_sub)
celltypes = names(celltypes[which(celltypes > 3)])

murine_EC = murine_EC[,murine_EC$tissue == "Brain"]
mtx = murine_EC@assays$RNA@data
rownames(mtx) = df$human[match(rownames(mtx), df$mouse)]
mtx = mtx[-which(is.na(rownames(mtx))), ]
mtx = mtx[rownames(mtx) %in% rownames(EC), ]
mtx2 = EC@assays$SCT@data; mtx2 = mtx2[rownames(mtx2) %in% rownames(mtx), ]
mtx = mtx[match(rownames(mtx2), rownames(mtx)), ]

data = matrix(0, nrow = 2, ncol = 10); data2 = matrix(0, nrow = 2, ncol = 10)

for (i in seq(2)) {
	g1 = apply(mtx[,murine_EC$ann_sub == celltypes[i]], 1, median)
	for (q in seq(length(unique(EC$short_celltype)))) {
		g2 = apply(mtx2[,EC$short_celltype == unique(EC$short_celltype)[q]], 1, median)
		res = cor.test(g1, g2, method="spearman", exact = F)
		data[i, q] = res$estimate; data2[i, q] = res$p.value
	}
}
rownames(data) = rownames(data2) = celltypes[1:2]; colnames(data) = colnames(data2) = unique(EC$short_celltype); 

library(fmsb)
data = as.data.frame(data);
data = rbind(rep(1,10), rep(0, 10), data)
colors_border=c('#0000ff','#ff0000')
colors_in=c('#B300ff','#B30000')
radarchart( data  , axistype=1 , 
    #custom polygon
    pcol=colors_border , plwd=4 , plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 
    )


legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=1.2, pt.cex=3)



Circulation_EC.Rdata
celltypes = unique(pbmc$celltype)

mtx = pbmc@assays$RNA@data
rownames(mtx) = df$human[match(rownames(mtx), df$mouse)]
mtx = mtx[-which(is.na(rownames(mtx))), ]
mtx = mtx[rownames(mtx) %in% rownames(EC), ]
mtx2 = EC@assays$SCT@data; mtx2 = mtx2[rownames(mtx2) %in% rownames(mtx), ]
mtx = mtx[match(rownames(mtx2), rownames(mtx)), ]

organoids = apply(mtx[,which(pbmc$tissue == "organoids")], 1, mean)
"merge_organoid_murineEC_circulation.Rdata"
celltypes = unique(pbmc$tissue); celltypes = celltypes[-11]
data = matrix(0, nrow = 1, ncol = 10); data2 = matrix(0, nrow = 1, ncol = 10)
for (i in seq(10)) {
	g1 = apply(mtx[,pbmc$tissue == celltypes[i]], 1, mean)
	g2 = organoids
	res = cor.test(g1, g2, method="spearman", exact = F)
	data[1, i] = res$estimate; data2[1, i] = res$p.value
}
rownames(data) = rownames(data2) = "organoids"; colnames(data) = colnames(data2) = celltypes; 



jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

data = matrix(0, nrow = 1, ncol = 10); data2 = matrix(0, nrow = 1, ncol = 10)
for (i in seq(10)) {
	g1 = apply(mtx[,pbmc$tissue == celltypes[i]], 1, mean); g1 = sort(g1, decreasing = T)
	g2 = sort(organoids, decreasing = T); g1 = names(g1[1:500]); g2 = names(g2[1:500])
	res = jaccard(g1, g2)
	data[1, i] = res;
}
rownames(data) = rownames(data2) = "organoids"; colnames(data) = colnames(data2) = celltypes; 

signatures = lapply(as.character(unique(markers$cluster)), function(i) {
	df = markers$gene[which(markers$cluster == i)][1:50]
	return(df)
	})
names(signatures) = as.character(unique(markers$cluster))
pbmc <- AddModuleScore_UCell(pbmc, features = signatures)




data = as.data.frame(data);
data = rbind(rep(1,10), rep(0, 10), data)
colors_border=brewer.pal(n=9, "Set1")
radarchart( data  , axistype=1 , 
    #custom polygon
    pcol=colors_border , plwd=4 , plty=1,
    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
    #custom labels
    vlcex=0.8 
    )


legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_border , text.col = "grey", cex=1.2, pt.cex=3)




mm2 = c("FABP7","PTN","CALB1","VEGFA","SLC1A3","SIX3","EN2","GAD2","DLX2","GABBR2",
	"GRIN2B","NEUROD6","HES6","STMN2","EOMES","PDGFRB","LUM","COL3A1",
	"TAC1","C7","ACTA2","TNMD","THBS4","MECOM","IGF1","SLC2A1","AUTS2","ENO1","PGK1","MECOM")
markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
top10 <- newmarkers %>% arrange(p_val_adj) %>%
    group_by(cluster) %>% 
    slice(1:5)


Idents(pbmc) = factor(pbmc$short_celltype, levels = c("Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs 1", "NPs 2", "PAst", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5", 
"Peri", "VLMCs 1", "VLMCs 2", "VLMCs 3", "SMCs", "Venous_ECs 1", "Venous_ECs 2", "Venous_ECs 3", "TCs", "Vas_ECs 1", "Brain_ECs", "PCs"   
))

top10$cluster = factor(top10$cluster, levels = haso)
top10 = top10[order(top10$cluster), ];
DoHeatmap(pbmc, features = top10$gene, raster=T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))



a = c("Peri", "Brain_ECs", "Ast 2", "VLMCs 1", "Ast 1", "PAst", "GluNs 4", "Venous_ECs 2", "GluNs 1", "TCs", "Vas_ECs 1", "GABAs", "PCs", 
"GluNs 3", "NPs 1", "GluNs 2", "VLMCs 2", "VLMCs 3", "Ast 4", "Ast 3", "GluNs 5", "SMCs", "Venous_ECs 1", "NPs 2") 
names(a) = celltypes
pbmc$short_celltype = a[pbmc$celltype2]
pbmc$short_celltype2 = gsub("\\ \\d*", "", pbmc$short_celltype)

haso = c("Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "Astrocytes 4", "Neural progenitors 1", "Neural progenitors 2", "Proliferative  astrocytes", 
	"GABAergic neurons", "Glutamatergic neurons 1", "Glutamatergic neurons 2", "Glutamatergic neurons 3", "Glutamatergic neurons 4", "Glutamatergic neurons 5", 
"Pericytes", "Vascular leptomeningeal cells 1(VLMCs)", "Vascular leptomeningeal cells 2(VLMCs)", "Vascular leptomeningeal cells 3(VLMCs)", 
"Smooth muscle cells (SMCs)", "Venous Endothelials 1", "Venous Endothelials 2", "Tendon cells", "Vascular endothelials 1", 
"Brain vascular endothelials", "Proliferative  cells")
Idents(pbmc) = factor(pbmc$celltype2, levels = haso)
a = paste(seq(24), haso, sep=". ")
names(a) = haso
Idents(pbmc) = a[pbmc$celltype2]



#gsea
library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens")

library(fgsea)
library(dplyr)
library(ggplot2)

pbmc.genes <- wilcoxauc(pbmc, 'celltype2')

m_df<- msigdbr(species = "Homo sapiens", category = "C5") # C5 for BP  H for hallmark
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)



rankGene <- function(cluster, df) {
	cluster.genes<- df %>%
  dplyr::filter(group == cluster) %>% arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
  ranks<- deframe(cluster.genes)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 200)
  fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
  fgseaResTidy <- fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj)
  return(fgseaResTidy)
}

celltypes = unique(pbmc$celltype2)

myres = lapply(celltypes, rankGene, pbmc.genes)

names(myres) = celltypes



plotEnrichment(fgsea_sets[["GOBP_CELL_PROJECTION_ORGANIZATION"]],
               ranks) + labs(title="GOBP_CELL_PROJECTION_ORGANIZATION")
























################ Sat May 21 15:24:34 2022
## Cell ECs
Cell_ECs_all_tissues.Rdata


pbmc = pbmc[,pbmc$EC_yn == "Yes"]


markers = FindAllMarkers(pbmc, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.05), ]
save(markers, file = "EC_markers.Rdata")

load('EC_markers.Rdata')

signatures = lapply(unique(markers$cluster), function(i) {
	df = markers$gene[which(markers$cluster == i)]; df = df[df %in% rownames(EC)]; df = df[1:5]; return(df)
	})

names(signatures) = unique(markers$cluster)


organoids = apply(EC@assays$SCT@data, 1, median); organoids = sort(organoids, decreasing = T)
organoids = organoids[-grep("^MT-", names(organoids))]; organoids = organoids[-grep("^RP[L|S]\\d*", names(organoids))]

score1 = c()
for (i in names(signatures)) {
	print(i)
	score1 = c(score1, 10000/sum(which(names(organoids) %in% signatures[[i]])[1:20], na.rm=T))
	signatures[[i]] = names(organoids)[which(names(organoids) %in% signatures[[i]])[1:20]]
}


library(UCell)
Idents(EC) = "orig"
EC <- AddModuleScore_UCell(EC, features = signatures)


library(reshape2)
library(ggpubr)
df = EC@meta.data[,27:37]
df = melt(df)
df$variable = gsub("_UCell", "", df$variable)
ggboxplot(df, x = "variable", y = "value", color = "variable", palette = mycolors) + theme(axis.text.y = element_text(size = 14)) + 
    theme(axis.text.x = element_text(size = 14, angle=45, hjust=1, vjust=1)) + theme(legend.position = "none") + xlab('') + ylab("")

EC1 = EC[,EC$short_celltype2 %in% c("Brain_ECs","Venous_ECs","Vas_ECs")]
EC1 <- AddModuleScore_UCell(EC1, features = signatures)
df = EC1@meta.data[,27:37]
df = melt(df)
df$variable = gsub("_UCell", "", df$variable)

ggboxplot(df, x = "variable", y = "value", color = "variable")
res = pairwise.wilcox.test(x = df$value, g = df$variable)






setwd("/Users/youzhen/.mounty/Inventory/Brain_organoid_scRNA/Brain_EC_ref/Normalized_EC_Circulation/")
load("EC.TSNE.markers.Robj")
load("EC.TSNE.Robj")



load('../../Third-filtered_feature_bc_matrix/filtered_feature_bc_matrix/filtered_combined_organoids_res.Rdata')
EC = pbmc[,pbmc$short_celltype2 %in% c('Peri',"Brain_ECs","VLMCs","Venous_ECs","Vas_ECs","SMCs")]
EC1 = pbmc[,pbmc$short_celltype2 %in% c("Brain_ECs","Venous_ECs","Vas_ECs")]

Idents(EC) = "orig"
Idents(EC1) = "orig"

load("EC.TSNE.markers.Robj")
homo = read.table("../../Third-filtered_feature_bc_matrix/filtered_feature_bc_matrix/homology_gene_MGI.trans.txt", header=T)
markers$gene = homo$human[match(markers$gene, homo$mouse)]
markers = na.omit(markers)



signatures = lapply(unique(markers$cluster), function(i) {
	df = markers$gene[which(markers$cluster == i)]; df = df[df %in% rownames(EC)]; df = df[-which(df %in% markers$gene[which(markers$cluster != i)])]
	})

names(signatures) = unique(markers$cluster)


organoids = apply(EC@assays$SCT@data, 1, median); organoids = sort(organoids, decreasing = T)
organoids = organoids[-grep("^MT-", names(organoids))]; organoids = organoids[-grep("^RP[L|S]\\d*", names(organoids))]
organoids1 = organoids[1:300]
score1 = c()
for (i in names(signatures)) {
	print(i)
	score1 = c(score1, 10000/sum(which(names(organoids) %in% signatures[[i]])[1:30], na.rm=T))
	signatures[[i]] = names(organoids)[which(names(organoids) %in% signatures[[i]])[1:30]]
}



EC <- AddModuleScore_UCell(EC, features = signatures)


library(reshape2)
library(ggpubr)
library(RColorBrewer)

df = EC@meta.data[,27:36]
df = melt(df)
df$variable = gsub("_UCell", "", df$variable)
mycolors = colorRampPalette(brewer.pal(8, "Dark2"))(10)
ggboxplot(df, x = "variable", y = "value", color = "variable", palette = mycolors) + theme(axis.text.y = element_text(size = 14)) + 
    theme(axis.text.x = element_text(size = 14, angle=45, hjust=1, vjust=1)) + theme(legend.position = "none") + xlab('') + ylab("")




#####
Idents(EC1) = "orig"
EC1 <- AddModuleScore_UCell(EC1, features = signatures)
df = EC1@meta.data[,27:36]
df = melt(df)
df$variable = gsub("_UCell", "", df$variable)

ggboxplot(df, x = "variable", y = "value", color = "variable")
res = pairwise.wilcox.test(x = df$value, g = df$variable)

## No ideal results





