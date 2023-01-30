
# scRNA-Seq Analysis merged three groups version
## Author Zhen Y

library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(harmony)
library(DoubletFinder)
options(stringsAsFactors = F)


### Step 1. Run the following codes on each sample of three groups (Control CCM UHC) 
x = Read10X(".", unique.features = TRUE, strip.suffix = FALSE) # change "." to the directory of each sample where raw data (barcode, features, matrix) were stored
pbmc <- CreateSeuratObject(counts = x, project = "brain_organoids", min.cells = 3, min.features = 200)
pbmc$batch = "batch2" # assign batch to each sample
pbmc$group = "Control_rep1" # assign group to each sample
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 10)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc <- SCTransform(pbmc, verbose = FALSE)
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes)

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


### DoubletFinder part cannot be run in my own PC
sweep.res.pbmc <- paramSweep_v3(pbmc, PCs = 1:15, sct = TRUE) # optimize params
sweep.stats_kidney <- summarizeSweep(sweep.res.pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_kidney)
barplot(bcmvn_pbmc$BCmetric, names.arg = bcmvn_pbmc$pK, las=2)
## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)      
nExp_poi <- round(0.075*nrow(pbmc@meta.data))  ## Assuming 7.5% doublet formation rate 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
saveRDS(seu_kidney, file = "PCA_SeuObj.RDS")
# ********************************************************



## Step 2. load Seurat object saved at step 1 from each sample
pbmc = readRDS("PCA_SeuObj.RDS")
pbmc <- FindNeighbors(pbmc, dims = 1:15, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE, resolution = 1)
ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ncols)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))

# optional
markers = FindAllMarkers(pbmc, logfc.threshold = 0.2, min.pct = 0.2, only.pos = T)
markers = markers[which(markers$p_val_adj < 0.01), ]
# ********************************************************


# step 3.
# repeat the step 2 on each group and merge biological replicates for example:
# CCM group
ONE = readRDS("./filtered_feature_bc_matrix_CCM_One/PCA_SeuObj.RDS")
SEC = readRDS("./filtered_feature_bc_matrix_CCM_Sec/PCA_SeuObj.RDS")
THI = readRDS("./filtered_feature_bc_matrix_CCM_Third/PCA_SeuObj.RDS")

# run step 2 and merge replicates
pbmc.big <- merge(ONE, y = c(SEC, THI), add.cell.ids = c("One", "Sec", "Third"), project = "CCM")

# rerun PCA and normalization and remove batch effects
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




# step 4.
# after generating all seurat objects of three groups:
# merge three groups and rerun PCA 
# for example
CCM = readRDS("CCM_merged_SeuObj.RDS")
UHC = readRDS("UHC_merged_SeuObj.RDS")
CON = readRDS("CON_merged_SeuObj.RDS")

# run step 2 and merge replicates
pbmc.big <- merge(CON, y = c(CCM, UHC), add.cell.ids = c("Control", "CCM", "UHC"), project = "Organoids")

# rerun PCA and normalization and remove batch effects
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
saveRDS(pbmc.big, file = "All_merged_SeuObj.RDS")
saveRDS(markers, file = "All_merged_raw_markers.RDS")




# step 5.
# annotate cell clusters based on marker genes provided below
# clustering will change at different working environment; So the following clusters were calculated on my side 
# assigned cell type to clusters based on your own calculation
pbmc.combined = readRDS("All_merged_SeuObj.RDS")
markers = readRDS("All_merged_raw_markers.RDS")
# this clustering were generated on my side. It would be changed in different working environment
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


# cell type with marker genes
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

# full list cell types and marker genes (Control + CCM + UHC); This cell types and markers should cover all cells from three groups
	"Proliferative Astrocytes", # SLC1A3 HIST1H4C
	"Astrocytes 1", ## FABP7 PTN PTPRZ1 SLC1A3
	"Astrocytes 2",  # VEGFA BNIP3 SLC16A3 SLC1A3 IGFBP2 PGK1 FAM162   cell migration/protein binding
	"Astrocytes 3",  # CXCL14
	"G2M phase Astrocytes", # CCNB2 BIRC5 CENPF SLC1A3
	"GABAergic neurons 1", # GAD1 GAD2 DLX2 DLX1
	"Proliferative GABAergic neurons 2",1), # GAD1 GAD2 DLX2 DLX1
	'Glutamatergic neurons 1', # GABBR2  GRIN2B BHLHE22
	'Glutamatergic neurons 2', #  GRIN2B NEUROD6 NEUROD2 BHLHE22 Precursors of Glut
	'Glutamatergic neurons 3', # NEUROD6 NEUROD2 BHLHE22 SLA Precursors of Glut
	'Glutamatergic neurons 4', # NHLH1 EOMES IGFBPL1  UBC_like
	'Glutamatergic neurons 5', # DDIT4 PGK1 STMN2 IGFBP2 MALAT1 
	"Neural progenitors 1", # C1QTNF4 CD9 HES1 SOX2
	"Neural progenitors 2", # PAX3 CNPY1 EN2
	"Brain vascular endothelials",2), ## SLC2A1 ENO1, PGK1, IGFBP2   cell migration/protein binding
	"Proliferative cells", # TOP2A MKI67 CENPF
	"Vascular leptomeningeal cells 1(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM  TAC1
	"Vascular leptomeningeal cells 2(VLMCs)", 1), # COL1A1 COL3A1 C7 LUM TSHZ2 IGF2 TGFBI
	"Tendon cells", # THBS4  SCXA
	'Pericytes 1', # PDGFRB RGS5
	'Pericytes 2', # PDGFRB BST2 MEF2C REN
	"Vascular endothelials 1", # CD34 CD93 CLDN5
	"Vascular endothelials 2", # KDR CLDN5 IL32 CD93
	"Vascular endothelials 3",  # PRND CLDN5  ESM1
	"Smooth muscle cells (SMCs)" # CCNB2 BIRC5 CENPF SLC1A3 
	

# assigned full and short cell type names to each cluster
names(ctn) = ct
cellT = ctn[as.character(Idents(pbmc.combined))]
names(cellT) = names(Idents(pbmc.combined))
Idents(pbmc.combined) = factor(cellT)
pbmc.combined$celltype = as.character(Idents(pbmc.combined))
pbmc.combined$celltype2 = pbmc.combined$celltype
pbmc.combined$celltype2 = gsub("\\ \\d", "", pbmc.combined$celltype)
shortN = unique(pbmc.combined$celltype);
# change full name to short name based on the full names generate from command `unique(pbmc.combined$celltype)`
names(shortN) = c("PAs", "Brain VEs", "Ast 2", "VLMCs 1", "Ast 1", "PAst", "GluNs 1", "TCs", "VEs 3", "GluNs 2", "GABAs", "PCs", "Ast 4", "GluNs 3", "NPs 1", 
	"Ast 3", "VLMCs 2", "NPs 2", "ECs")
pbmc.combined$short_celltype = names(shortN)[match(pbmc.combined$celltype, shortN)]

ncols = length(unique(Idents(pbmc.combined)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
DimPlot(pbmc.combined, label = TRUE, cols = mycolors, label.size = 6, repel = T) + NoLegend() + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))


# making heatmap
markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10




# short celltype changed order and cell types based on your own clusters
Idents(pbmc.combined) = factor(Idents(pbmc.combined), levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))


top10$cluster = factor(top10$cluster, levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs'))
top10 = top10[order(top10$cluster), ]

DoHeatmap(pbmc.combined, features = top10$gene, raster = T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))


saveRDS(pbmc.combined, file = "combined_organoids_celltype.Rdata")
write.table(markers, "celltype_enriched_markers.xls", sep="\t", quote=F, col.names=T, row.names=F)


# step 6.
### cell cell communication
### analysis
require(pbmcapply)
pbmc = readRDS("combined_organoids_celltype.Rdata")
markers = read.table("celltype_enriched_markers.xls", header=T, sep="\t") # 
Idents(pbmc) = pbmc$short_celltype
Idents(pbmc) = factor(Idents(pbmc), levels = c("NPs 1", "NPs 2", "Ast 1", "Ast 2", "Ast 3", "Ast 4", "PAst", 'GABAs', 
	'GluNs 1', 'GluNs 2', 'GluNs 3', 'VLMCs 1', 'VLMCs 2', 'Brain VEs', 'VEs 3', 'PAs','TCs', 'PCs', 'ECs')) # factorize cell types change cell types accordingly
shortN = unique(Idents(pbmc))
names(shortN) = unique(pbmc$short_celltype)
markers$cluster = names(shortN)[match(markers$cluster, shortN)]

# normalized
countTable = pbmc@assays$RNA@counts
LR = read.csv("mouse_human_lrdb_version2.csv") # ligand receptor database
countTable = as.matrix(countTable[rownames(countTable) %in% unique(c(LR$Ligand.Human, LR$Receptor.Human)), ])
cTMM = TMMnormalization(countTable) # TMM normalization
celltypes = levels(unique(Idents(pbmc)))


celltype1 = ""; celltype2 = ""

condition = as.character(Idents(pbmc));



#  make all combination df
workdf = expand.grid(celltypes, celltypes)
names(workdf) = c("Sender", "Receiver")

# run main func of CCC
totaldf = data.frame()
CCTObj = list()
for (i in seq(nrow(workdf))) {
    print(paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> "))
    tempdf = cellcrosstalker(cTMM, workdf$Sender[i], workdf$Receiver[i], 0.95, condition, 1)
    tdf1 = tempdf@result_df; tdf1$Sender = workdf$Sender[i]; tdf1$Receiver = workdf$Receiver[i];
    totaldf = rbind(totaldf, tdf1)
    CCTObj[[paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> ")]] = c(tempdf@paramL, tempdf@paramR)
}



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

save(totaldf, CCTObj, file = "CCC_CellCrossTalker2.Rdata") # save CCC results




