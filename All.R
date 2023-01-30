
options(stringsAsFactors=F)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(DoubletFinder)

args = commandArgs(trailingOnly=TRUE)

runDoubletFinder <- function(pbmc, group, batch) {
	DefaultAssay(pbmc) = "SCT"
	print("DefaultAssay ...")
	sweep.res.pbmc <- paramSweep_v3(pbmc, PCs = 1:15, sct = TRUE, num.cores=4) # optimize params
	sweep.stats_kidney <- summarizeSweep(sweep.res.pbmc, GT = FALSE)
	print("Finding pK ...")
	bcmvn_pbmc <- find.pK(sweep.stats_kidney)
	## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
	annotations <- pbmc@meta.data$seurat_clusters
	homotypic.prop <- modelHomotypic(annotations)      
	nExp_poi <- round(0.075*nrow(pbmc@meta.data))
	nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
	## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
	print("Run DoubletFinder with varying classification stringencies ...")
	optimalpK = as.numeric(as.character(bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)]))
	seu_kidney <- doubletFinder_v3(pbmc, PCs = 1:15, pN = 0.25, pK = optimalpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
	print(table(seu_kidney@meta.data[,ncol(seu_kidney@meta.data)]))
	colnames(seu_kidney@meta.data)[ncol(seu_kidney@meta.data)] = "DF_results"
	print("saving Seurat Obj ... ")
	saveRDS(seu_kidney, file = "SeuratObj_removeDB.RDS")
}


runSeurat = function(dir, group, batch) {
	setwd(dir)
	x = Read10X(".", unique.features = TRUE, strip.suffix = FALSE) # change "." to the directory of each sample where raw data (barcode, features, matrix) were stored
	pbmc <- CreateSeuratObject(counts = x, project = "brain_organoids", min.cells = 3, min.features = 200)
	pbmc$batch = batch # assign batch to each sample
	pbmc$group = group # assign group to each sample
	pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
	pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	print("SCTransforming ...")
	pbmc <- SCTransform(pbmc, verbose = FALSE)
	print("CellCycling score ...")
	pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)
	print("Running PCA ...")
	pbmc <- RunPCA(pbmc, verbose = FALSE)
	print("Running UMAP ...")
	pbmc <- RunUMAP(pbmc, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
	pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose = FALSE)
	pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 1)
	print("Running DoubletFinder ...")
	pbmc <- runDoubletFinder(pbmc, group, batch)
	return(pbmc)
}

workdir = "/public/home/gonglianggroup/zheny/temp/Brain_Organoid/"
group = args[1]
batch = args[2]
dirs = paste0(workdir, "filtered_feature_bc_matrix_", group, "_", batch)
runSeurat(dirs, group, batch)



groups=(CON CON CON CCM CCM CCM UHC UHC UHC)
batches=(One Sec Third One Sec Third One Sec Third)

groups=("CON", "CCM", "UHC")

for (group in groups) {
	ONE = readRDS(paste0(workdir, "filtered_feature_bc_matrix_", group, "_One/SeuratObj_removeDB.RDS"))
	ONE = subset(ONE, subset = DF_results == "Singlet")
	SEC = readRDS(paste0(workdir, "filtered_feature_bc_matrix_", group, "_Sec/SeuratObj_removeDB.RDS"))
	SEC = subset(SEC, subset = DF_results == "Singlet")
	THI = readRDS(paste0(workdir, "filtered_feature_bc_matrix_", group, "_Third/SeuratObj_removeDB.RDS"))
	THI = subset(THI, subset = DF_results == "Singlet")
	print(length(ONE$orig.ident))
	print(length(SEC$orig.ident))
	print(length(THI$orig.ident))
	print("merging ...")
	pbmc.big <- merge(ONE, y = c(SEC, THI), add.cell.ids = c("One", "Sec", "Third"), project = group)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes

	pbmc.big <- SCTransform(pbmc.big, verbose = FALSE)
	pbmc.big <- CellCycleScoring(pbmc.big, s.features = s.genes, g2m.features = g2m.genes)
	print("Running PCA ...")
	pbmc.big <- RunPCA(pbmc.big, verbose = FALSE)

	pbmc.big = RunHarmony(pbmc.big, group.by.vars = "batch")
	pbmc.big <- RunUMAP(pbmc.big, dims = 1:15, verbose = FALSE, n.neighbors = 30, spread = 1, min.dist = 0.1)
	pbmc.big <- FindNeighbors(pbmc.big, dims = 1:15, verbose = FALSE)
	pbmc.big <- FindClusters(pbmc.big, verbose = FALSE, resolution = 1)
	saveRDS(pbmc.big, file = paste0("/public/home/gonglianggroup/zheny/temp/Brain_Organoid/results/", group, "_merged_SeuObj.RDS"))
	print("Done")
}






# step 3.
workdir = "/public/home/gonglianggroup/zheny/temp/Brain_Organoid/results/"
setwd(workdir)

CCM = readRDS("CCM_merged_SeuObj.RDS")
CON = readRDS("CON_merged_SeuObj.RDS")
UHC = readRDS("UHC_merged_SeuObj.RDS")
pbmc.big <- merge(CON, y = c(CCM, UHC), add.cell.ids = c("CON", "CCM", "UHC"), project = "brain_organoids")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

threds = mean(pbmc.big$nCount_RNA) + 3*sd(pbmc.big$nCount_RNA)
pbmc.big = subset(pbmc.big, subset = nCount_RNA < threds)
pbmc.big <- SCTransform(pbmc.big, verbose = FALSE)
pbmc.big <- CellCycleScoring(pbmc.big, s.features = s.genes, g2m.features = g2m.genes)
print("Running PCA ...")
pbmc.big <- RunPCA(pbmc.big, verbose = FALSE)

pbmc.big <- RunUMAP(pbmc.big, dims = 1:15, verbose = FALSE, n.neighbors = 20, spread = 1, min.dist = 0.2)
pbmc.big <- FindNeighbors(pbmc.big, dims = 1:15, verbose = FALSE)
pbmc.big <- FindClusters(pbmc.big, verbose = FALSE, resolution = 1.5)


ncols = length(unique(Idents(pbmc.big)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
pdf("All_merged_DimPlot.pdf", w=10, h=8)
DimPlot(pbmc.big, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc.big, cols = c("#1C82AD","#03C988", "#F273E6"), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc.big, cols = c("red", "navy", "darkgreen"), group.by = 'group') + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()

library(future)
options(future.globals.maxSize=18485760000) ## 10000*1024^2
plan("multicore", workers = 8)
markers = FindAllMarkers(pbmc.big, logfc.threshold = 0.2, min.pct = 0.2, only.pos = F)
markers = markers[which(markers$p_val_adj < 0.01), ]
saveRDS(pbmc.big, file = "All_merged_SeuObj.RDS")
saveRDS(markers, file = "All_merged_raw_markers.RDS")


pdf("Feature_CellCycleScore.pdf", w = 16, h=7)
FeaturePlot(pbmc.big, features = c("S.Score", "G2M.Score"), pt.size = 0, cols = c("grey", "#db8495", "#b31031"), ncol = 2, raster = TRUE)
dev.off()


ct = c(
2,6,
1,8,
4,
23,
38, 
31,
29,
37, 13,
30,
32,
21,
28,
14,19,
15,5, 40,
11,
9,
12,33,
18,
22,
24,
39,
20,35
)


ctn = c(
	rep("Astrocytes 4", 2), # CXCL14 NOG
	rep("Astrocytes 1", 2), ## FABP7 PTN PTPRZ1 SLC1A3
	"Astrocytes 2", 	# CALB1 CRYAB 
	"Astrocytes 3",  # VEGFA BNIP3 SLC16A3 SLC1A3 IGFBP2 PGK1 FAM162   cell migration/protein binding

	"Neural progenitors 3", # PAX3 CNPY1 EN2
	"Neural progenitors 1", # CENPU HELLS GINS2 HES5   C1QTNF4 CD9 HES1 SOX2
	"Neural progenitors 2", # PTTG1 SOX2 HES5

	rep("G2M phase Astrocytes 1", 2), # CCNB2 BIRC5 CENPF SLC1A3
	rep("G2M phase Astrocytes 2", 1), # CCNB2 BIRC5 CENPF SLC1A3
	rep("S phase Astrocytes 1", 1), # UBE2C TOP2A SLC1A3
	rep("S phase Astrocytes 2", 1), # UBE2C TOP2A SLC1A3

	"GABAergic neurons", # GAD1 GAD2 DLX2 DLX1
	
	rep('Glutamatergic neurons 1', 2), # GABBR2  GRIN2B BHLHE22
	rep('Glutamatergic neurons 2', 3), #  GRIN2B NEUROD6 NEUROD2 BHLHE22 Precursors of Glut
	rep('Glutamatergic neurons 3', 1), # NEUROD6 NEUROD2 BHLHE22 SLA Precursors of Glut
	rep('Glutamatergic neurons 4', 1), # DDIT4   cell migration/protein binding
	rep('Glutamatergic neurons 5', 2), # NHLH1 EOMES   UBC_like
	rep('Glutamatergic neurons 6', 1), # STMN2, DCX
	rep('Glutamatergic neurons 7', 1), # TMSB4X, ACTG1
	
	"Vascular endothelials 1", # CD34 CD93 CLDN5
	"Vascular endothelials 2", # KDR CLDN5 IL32 CD93

	rep("Proliferative cells", 2)) # TOP2A MKI67 CENPF

	"Vascular endothelials 1", ## IGFBP3 ALDOA, ENO1
	"Perivascular adipocyte", ## IGF1, MEOX2, MECOM
	"Tendon cells",  ##  TNMD, THBS4
	"Vascular endothelials 2" ## VEGFA, ALDOA, ENO1

 	rep("Vascular leptomeningeal cells 1(VLMCs)", 2), ## LUM, COL3A1, COL1A1, IGF2, PLAT
 	"Vascular leptomeningeal cells 2(VLMCs)", # GUCY1A3 C7 CELF3 
 	"Vascular endothelials 3" ## EIF1, ALDOA, ENO1 FTL

	
celltype1 = c("Tendon cells", "Brain vascular endothelials", "Astrocytes 2", "Pericytes 2", "Astrocytes 4", "Smooth muscle cells (SMCs)", 
	"G2M phase Astrocytes 2", "Glutamatergic neurons 1", "Vascular endothelials 3", "Pericytes 1", "Vascular leptomeningeal cells 3 (VLMCs)", 
	"Vascular endothelials 4", "GABAergic neurons", "Vascular leptomeningeal cells 2 (VLMCs)", "Proliferative cells", "Glutamatergic neurons 6", 
	"S phase Astrocytes 2", "Pericytes 3", "Astrocytes 1", "Glutamatergic neurons 2", "Glutamatergic neurons 4", "Neural progenitors 2", 
	"Glutamatergic neurons 3", "Glutamatergic neurons 5", "Vascular leptomeningeal cells 1 (VLMCs)", "G2M phase Astrocytes 1", 
	"Vascular endothelials 1", "Tendon cells", "Glutamatergic neurons 7", "Vascular leptomeningeal cells 1 (VLMCs)", "S phase Astrocytes 1", 
	"noise", "Neural progenitors 1", "Neural progenitors 3", "Smooth muscle cells (SMCs)", "Vascular endothelials 2", "Vascular endothelials 5", 
	"Proliferative cells", "Vascular endothelials 6", "Brain vascular endothelials")

celltype2 = c("TCs", "Brain VEs", "Ast 2", "Peri 2", "Ast 4", "SMCs", 
	"G2M Ast 2", "GluNs 1", "VEs 3", "Peri 1", "VLMCs 3", 
	"VEs 4", "GABAs", "VLMCs 2", "PCs", "GluNs 6", 
	"S Ast 2", "Peri 3", "Ast 1", "GluNs 2", "GluNs 4", "NPs 2", 
	"GluNs 3", "GluNs 5", "VLMCs 1", "G2M Ast 1", 
	"VEs 1", "TCs", "GluNs 7", "VLMCs 1", "S Ast 1", 
	"noise", "NPs 1", "NPs 3", "SMCs", "VEs 2", "VEs 5", 
	"PCs", "VEs 6", "Brain VEs")

celltype3 = c("TCs", "Brain VEs", "Ast", "Peri", "Ast", "SMCs", 
	"G2M Ast", "GluNs", "VEs", "Peri", "VLMCs", 
	"VEs", "GABAs", "VLMCs", "PCs", "GluNs", 
	"S Ast", "Peri", "Ast", "GluNs", "GluNs", "NPs", 
	"GluNs", "GluNs", "VLMCs 1", "G2M Ast", 
	"VEs", "TCs", "GluNs", "VLMCs", "S Ast", 
	"noise", "NPs", "NPs", "SMCs", "VEs", "VEs", 
	"PCs", "VEs", "Brain VEs")


feaPlo <- function(gene) {
	pdf(paste0(gene, "_FeaturePlot.pdf"), w=8, h=6)
	p = FeaturePlot(pbmc.big, gene, cols = c("grey", "#dc143c"), raster = TRUE)
	print(p)
	dev.off()
}	

# assigned full and short cell type names to each cluster
names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc.big))]
names(cellT) = names(Idents(pbmc.big))
Idents(pbmc.big) = factor(cellT)
pbmc.big$celltype = as.character(Idents(pbmc.big))

pbmc.big$full_celltype = factor(pbmc.big$full_celltype, levels = c(
"Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "Astrocytes 4", "Neural progenitors 1", "Neural progenitors 2", "Neural progenitors 3",
"GABAergic neurons","Glutamatergic neurons 1", "Glutamatergic neurons 2", "Glutamatergic neurons 3", "Glutamatergic neurons 4", 
"Glutamatergic neurons 5", "Glutamatergic neurons 6", "Glutamatergic neurons 7", "G2M phase Astrocytes 1", "G2M phase Astrocytes 2",
"S phase Astrocytes 1", "S phase Astrocytes 2", "Vascular leptomeningeal cells 1 (VLMCs)", "Vascular leptomeningeal cells 2 (VLMCs)",
"Vascular leptomeningeal cells 3 (VLMCs)", "Vascular endothelials 1", "Vascular endothelials 2", "Vascular endothelials 3", "Vascular endothelials 4",
"Vascular endothelials 5", "Vascular endothelials 6", "Brain vascular endothelials", "Pericytes 1", "Pericytes 2", "Pericytes 3",  "Smooth muscle cells (SMCs)",
"Tendon cells", "Proliferative cells"      
	))
pdf("All_full_celltype_DimPlot.pdf", w = 14, h=6)
ncols = length(unique(pbmc.big$full_celltype))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.big, group.by = "full_celltype", cols = mycolors, label.size = 6, repel = T) + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()

pbmc.big$short_celltype = factor(pbmc.big$short_celltype, levels = c(
"Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs 1", "NPs 2", "NPs 3",
"GABAs","GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", 
"GluNs 5", "GluNs 6", "GluNs 7", "G2M Ast 1", "G2M Ast 2",
"S Ast 1", "S Ast 2", "VLMCs 1", "VLMCs 2",
"VLMCs 3", "VEs 1", "VEs 2", "VEs 3", "VEs 4",
"VEs 5", "VEs 6", "Brain VEs", "Peri 1", "Peri 2", "Peri 3",  "SMCs",
"TCs", "PCs"      
	))

pdf("All_short_celltype_DimPlot.pdf", w = 9, h=6)
ncols = length(unique(pbmc.big$short_celltype))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mycolors = sample(mycolors)
DimPlot(pbmc.big, group.by = "short_celltype", cols = mycolors, label = TRUE, label.size = 6, repel = T) + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()


pbmc.big$merged_celltype = factor(pbmc.big$merged_celltype, levels = c(
"Ast", "NPs", 
"GABAs","GluNs", 
"G2M Ast", 
"S Ast", "VLMCs", "VEs", 
"Brain VEs", "Peri", "SMCs",
"TCs", "PCs"      
	))

pdf("All_merged_celltype_DimPlot.pdf", w = 8, h=6)
ncols = length(unique(pbmc.big$merged_celltype))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.big, group.by = "merged_celltype", cols = mycolors, label = TRUE, label.size = 6, repel = T) + NoLegend() +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()

Idents(pbmc.big) = pbmc.big$full_celltype
markers = FindAllMarkers(pbmc.big, logfc.threshold = 0.2, min.pct = 0.2, only.pos = F)
markers = markers[which(markers$p_val_adj < 0.01), ]

saveRDS(pbmc.big, file = "combined_organoids_celltype.Rdata")
write.table(markers, "celltype_enriched_markers.xls", sep="\t", quote=F, col.names=T, row.names=F)





#Wed Jan 25 12:00:43 2023
# making heatmap
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
markers = readRDS("CCM_UHC_merged_celltype_markers.RDS")
markers = markers[which(markers$avg_log2FC > 0), ]

markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top10

top10 = top10[order(top10$cluster), ]


, raster = T

pdf("celltype_HeatMap.pdf", w = 15, h=17)
DoHeatmap(pbmc, features = top10$gene) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))
dev.off()




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

mycolors = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

pdf("CCM_UHC_markers_stackedVln.pdf", w = 15, h = 30)
mymarker = c('SLC1A3', 'HES5', 'VEGFA', 'BNIP3', 'CXCL14', 'CCNB2','GAD2','DLX1', 'GABBR2', 'GRIN2B', 'BHLHE22', 'NEUROD2', 'EOMES', 'DDIT4', 'LUM', 'COL3A1', 'TSHZ2', 
	'PDGFRB', 'RGS5', 'REN', 'THBS4', 'SLC2A1', 'CLDN5', 'CD34', 'KDR', 'PRND', 'ACTA2', 'TAGLN', 'HIST1H1B', 'PTTG1', 'TOP2A')
StackedVlnPlot(obj = pbmc, features = mymarker[1:16], cols = mycolors)
StackedVlnPlot(obj = pbmc, features = mymarker[17:31], cols = mycolors)
dev.off()

### cell type proportions
df = pbmc@meta.data
df = df[,c('Sub_celltype', 'group')]
plotdf = as.data.frame(df %>% count(Sub_celltype, group, name="Count"))

plotdf$Count[which(plotdf$group == "CCM")] = plotdf$Count[which(plotdf$group == "CCM")]/20870
plotdf$Count[which(plotdf$group == "UHC")] = plotdf$Count[which(plotdf$group == "UHC")]/21729
plotdf$Count = round(plotdf$Count * 100, digit = 3)

plotdf$Sub_celltype = factor(plotdf$Sub_celltype, levels = levels(Idents(pbmc)))

pdf("CCM_UHC_cell_propotions.pdf", w=15, h=8)
ggplot(plotdf, aes(x = Sub_celltype, y = Count, fill = group)) + geom_bar(stat="identity", width = 0.5, position=position_dodge()) + 
scale_fill_manual(values = c("navy","red")) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))

# level2
df = pbmc@meta.data
df = df[,c('short_GlobalCT', 'group')]
plotdf = as.data.frame(df %>% count(short_GlobalCT, group, name="Count"))

plotdf$Count[which(plotdf$group == "CCM")] = plotdf$Count[which(plotdf$group == "CCM")]/20870
plotdf$Count[which(plotdf$group == "UHC")] = plotdf$Count[which(plotdf$group == "UHC")]/21729
plotdf$Count = round(plotdf$Count * 100, digit = 3)

plotdf$short_GlobalCT = factor(plotdf$short_GlobalCT, levels = c(
	"Ast", "NPs", "GABAs", "GluNs", "Pro Asts", "Pro GABAs", "VLMCs", "VEs", "Peri", "SMCs", "fibroblast", "PCs"
	))


ggplot(plotdf, aes(x = short_GlobalCT, y = Count, fill = group)) + geom_bar(stat="identity", width = 0.5, position=position_dodge()) + 
scale_fill_manual(values = c("navy","red")) + 
theme_classic() + theme(axis.ticks = element_line(size=1, color="black"), axis.line = element_line(size=1, color="black"), 
	axis.text = element_text(size=18, color="black"), axis.title = element_text(size=22, color='black')) + xlab("") + 
ylab("Percentage (%)") + scale_y_continuous(expand=c(0,0))

dev.off()


