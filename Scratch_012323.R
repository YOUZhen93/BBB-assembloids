clusterProfiler



# CCM UHC:
options(stringsAsFactors=F)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(multienrichjam)
library(R.utils)
library(openxlsx)

dofunctionenrich <- function(pbmc, cluster) {
	print(cluster)
	y = markers[which(markers$cluster == cluster), ]
	geneList = y$avg_log2FC[which(y$avg_log2FC>0)]
	names(geneList) = y$gene[which(y$avg_log2FC>0)]
	print(length(geneList))
	backgrounds = bitr(rownames(pbmc), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

	eg = bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
	eg$FC = geneList[match(eg$SYMBOL, names(geneList))]
	geneList = eg$FC; names(geneList) = eg$ENTREZID
	keggup1 = enrichKEGG(names(geneList), universe = backgrounds$ENTREZID, organism = "hsa", 
	                    keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 7, 
	                    qvalueCutoff = 0.05, use_internal_data = F)
	keggup1 = keggup1@result
	keggup1$regulation = "UP"; keggup1$class = "KEGG"
	ego1 <- enrichGO(gene = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05, 
                readable      = TRUE)
    ego1 = ego1@result
    ego1$regulation = "UP"; ego1$class = "GO"
    res1 = rbind(ego1, keggup1)


	y = markers[which(markers$cluster == cluster), ]
	geneList = y$avg_log2FC[which(y$avg_log2FC < 0)]
	names(geneList) = y$gene[which(y$avg_log2FC < 0)]
	backgrounds = bitr(rownames(pbmc), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

	eg = bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
	eg$FC = geneList[match(eg$SYMBOL, names(geneList))]
	geneList = eg$FC; names(geneList) = eg$ENTREZID
	keggup2 = enrichKEGG(names(geneList), universe = backgrounds$ENTREZID, organism = "hsa", 
	                    keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = 'BH', minGSSize = 7, 
	                    qvalueCutoff = 0.05, use_internal_data = F)
	keggup2 = keggup2@result
	keggup2$regulation = "DOWN"; keggup2$class = "KEGG"
	ego2 <- enrichGO(gene = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05, 
                readable      = TRUE)
    ego2 = ego2@result
    ego2$regulation = "DOWN"; ego2$class = "GO"
    res2 = rbind(ego2, keggup2)

    res = rbind(res1, res2)
    return(res)    
}

markers$cluster = factor(markers$cluster, levels = c(
"Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "Astrocytes 4", "Neural progenitors 1", "Neural progenitors 2", "GABAergic neurons", 
"Glutamatergic neurons 1", "Glutamatergic neurons 2", "Glutamatergic neurons 3", "Glutamatergic neurons 4", "Glutamatergic neurons 5", 
"Proliferative  astrocytes", "Vascular leptomeningeal cells 1(VLMCs)", "Vascular leptomeningeal cells 2(VLMCs)", 
"Vascular leptomeningeal cells 3(VLMCs)", "vascular endothelial cells 1", "venous endothelial cells 1", "venous endothelial cells 2",
"Brain vascular endothelials", "Pericytes", "Smooth muscle cells (SMCs)", "fibroblast", "Proliferative  cells"
))

# Control
x = lapply(levels(markers$cluster), function(i) dofunctionenrich(pbmc, i))
names(x) = c("Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs 1", "NPs 2", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "VLMCs 1", "VLMCs 2", "VLMCs 3", "VEs 1", "Venous ECs 1", "Venous ECs 2", "Brain VEs", "Peri", "SMCs", "Fibro", "PCs")
openxlsx::write.xlsx(x, "GO_KEGG_analysis_result.xlsx")

# UHC CCM

x = lapply(levels(markers$cluster), function(i) dofunctionenrich(pbmc, i))
names(x) = c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "Pro GABAs", "VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fibro", "PCs")
openxlsx::write.xlsx(x, "CCM_UHC_GO_KEGG_analysis_result.xlsx")











# UHC CCM
pbmc = readRDS("./Figures/CCM_UHC_merged_SeuratObj.RDS")
pbmc$Sub_celltype = factor(pbmc$Sub_celltype, levels = c(
"Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "neural progenitors", "GABAergic neurons 1", "Glutamatergic neurons 1", "Glutamatergic neurons 2", 
"Glutamatergic neurons 3", "Glutamatergic neurons 4", "Glutamatergic neurons 5", "Proliferative Astrocytes", "Proliferative GABAergic neurons 2",
"Vascular leptomeningeal cells 1(VLMCs)", "Vascular leptomeningeal cells 2(VLMCs)", "vascular endothelial cells 1", "vascular endothelial cells 2",
"vascular endothelial cells 3", "Pericytes 1", "Pericytes 2", "Smooth muscle cells (SMCs)", "fibroblast", "Proliferative cells"   
))
Idents(pbmc)= pbmc$Sub_celltype

ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mycolors = sample(mycolors)
pdf("./Figures/CCM_UHC_merged_DimPlot.pdf", w=15, h=8)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("#276893","#c65306", "#dbce54", "#5a5c5b", "#c35655", "#2ec3e7"), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("red", "navy"), group.by = 'group') + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()


# Control
pbmc$celltype2 = factor(pbmc$celltype2, levels = c(
"Astrocytes 1", "Astrocytes 2", "Astrocytes 3", "Astrocytes 4", "Neural progenitors 1", "Neural progenitors 2", "GABAergic neurons", 
"Glutamatergic neurons 1", "Glutamatergic neurons 2", "Glutamatergic neurons 3", "Glutamatergic neurons 4", "Glutamatergic neurons 5", 
"Proliferative  astrocytes", "Vascular leptomeningeal cells 1(VLMCs)", "Vascular leptomeningeal cells 2(VLMCs)", 
"Vascular leptomeningeal cells 3(VLMCs)", "vascular endothelial cells 1", "venous endothelial cells 1", "venous endothelial cells 2",
"Brain vascular endothelials", "Pericytes", "Smooth muscle cells (SMCs)", "fibroblast", "Proliferative  cells"
))
Idents(pbmc)= pbmc$celltype2


ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mycolors = sample(mycolors)
pdf("./Figures/Control_merged_DimPlot.pdf", w=15, h=8)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("#276893","#c65306", "#dbce54"), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
dev.off()



# Tue Jan 24 15:08:45 2023
# heatmap
markers = markers[which(markers$avg_log2FC > 0), ]
markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10

pdf("Control_merged_HeatMap.pdf", w=14, h=14)
DoHeatmap(pbmc, features = top10$gene, raster = T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))
dev.off()










# cellcrosstalker
#Tue Jan 24 15:51:29 2023
require(pbmcapply)
pbmc = readRDS("Control_merged_SeuratObj.RDS")
markers = readRDS("Control_merged_celltype_markers.RDS")

source("DGBD_ready2use.R")
# normalized
countTable = pbmc@assays$RNA@counts
LR = read.csv("/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/mouse_human_lrdb_version2.csv")
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








pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
markers = readRDS("CCM_UHC_merged_celltype_markers.RDS")
markers = markers[which(markers$avg_log2FC > 0), ]
markers1 = markers
source("DGBD_ready2use.R")
# normalized
LR = read.csv("/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/mouse_human_lrdb_version2.csv")

for (group in c("CCM", "UHC")) {
	print(group)
	pbmc1 = pbmc[,pbmc$group == group]
	celltypes = data.frame(table(Idents(pbmc1)))
	celltypes = celltypes$Var1[which(celltypes$Freq > 10)]
	celltypes = factor(celltypes, levels = celltypes)
	pbmc1 = pbmc1[,Idents(pbmc1) %in% celltypes]
	Idents(pbmc1) = factor(Idents(pbmc1), levels = celltypes)
	markers = markers1[which(markers1$cluster %in% celltypes), ]
	markers$cluster = factor(markers$cluster, levels = celltypes)
	countTable = pbmc1@assays$RNA@counts
	countTable = as.matrix(countTable[rownames(countTable) %in% unique(c(LR$Ligand.Human, LR$Receptor.Human)), ])
	cTMM = TMMnormalization(countTable)
	celltype1 = ""; celltype2 = ""
	condition = as.character(Idents(pbmc1));
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
	saveRDS(totaldf, file = paste0("CellCrossTalker_", group, "_res.RDS"))
	saveRDS(CCTObj, file = paste0("CellCrossTalker_", group, "_CCTObj.RDS"))
}




#  make all combination df
# Wed Jan 25 20:43:33 2023
options(stringsAsFactors=F)
library(tidyverse)
library(reshape2)
library(circlize)
library(ComplexHeatmap)

CCTObj = readRDS("CellCrossTalker_Control_CCTObj.RDS")
totaldf = readRDS("CellCrossTalker_Control_res.RDS")

# Control
x = c("Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs 1", "NPs 2", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "VLMCs 1", "VLMCs 2", "VLMCs 3", "VEs 1", "Venous ECs 1", "Venous ECs 2", "Brain VEs", "Peri", "SMCs", "Fib", "PCs")
names(x) = levels(totaldf$Sender)

totaldf$Sender = as.character(totaldf$Sender)
totaldf$Receiver = as.character(totaldf$Receiver)

totaldf$Sender = x[match(totaldf$Sender, names(x))]
totaldf$Receiver = x[match(totaldf$Receiver, names(x))]



totaldf$cell_cell = paste(totaldf$Sender, totaldf$Receiver, sep = ">>>")
totaldf$LR = paste(totaldf$Ligand, totaldf$Receptor, sep = "---")
totaldf$LRscore = -log10(totaldf$KLscore) ## The smaller the KL divergence the greater connections
totaldf$LRscore = totaldf$LRscore + abs(min(totaldf$LRscore)) + 0.1
totaldf = arrange(totaldf, cell_cell, -LRscore)

saveRDS(totaldf, file = "CellCrossTalker_Control_res.RDS")





CCC_table = aggregate(totaldf$LRscore, by = list(totaldf$cell_cell), FUN=sum)
names(CCC_table) = c("cell_cell", "Aggr_score")
CCC_table = CCC_table[order(CCC_table$Aggr_score, decreasing=T), ]
# CCC_table$Aggr_score = CCC_table$Aggr_score + abs(min(CCC_table$Aggr_score)) + 1

x = CCC_table
x$Sender = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
x$Receiver = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
mtx1 = acast(x, Sender~Receiver, value.var="Aggr_score") # convert to nonsymmetric matrix
fords = c("Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs 1", "NPs 2", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "VLMCs 1", "VLMCs 2", "VLMCs 3", "VEs 1", "Venous ECs 1", "Venous ECs 2", "Brain VEs", "Peri", "SMCs", "Fib", "PCs")
reindex = match(fords, colnames(mtx1))
mtx1 = mtx1[reindex, reindex]
write.csv(mtx1, "Control_aggregated_CCC_matrix.csv", quote = F, col.names=T, row.names=T)



# UHC CCM
CCTObj = readRDS("CellCrossTalker_CCM_CCTObj.RDS")
totaldf = readRDS("CellCrossTalker_CCM_res.RDS")

x = c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "Pro GABAs", "VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "Peri 1", "Peri 2", "Fib", "PCs")
names(x) = levels(totaldf$Sender)

totaldf$Sender = as.character(totaldf$Sender)
totaldf$Receiver = as.character(totaldf$Receiver)

totaldf$Sender = x[match(totaldf$Sender, names(x))]
totaldf$Receiver = x[match(totaldf$Receiver, names(x))]

totaldf$cell_cell = paste(totaldf$Sender, totaldf$Receiver, sep = ">>>")
totaldf$LR = paste(totaldf$Ligand, totaldf$Receptor, sep = "---")
totaldf$LRscore = -log10(totaldf$KLscore + 0.0001) ## The smaller the KL divergence the greater connections
totaldf$LRscore = totaldf$LRscore + abs(min(totaldf$LRscore)) + 0.1
totaldf = arrange(totaldf, cell_cell, -LRscore)

saveRDS(totaldf, file = "CellCrossTalker_CCM_res.RDS")

CCC_table = aggregate(totaldf$LRscore, by = list(totaldf$cell_cell), FUN=sum)
names(CCC_table) = c("cell_cell", "Aggr_score")
CCC_table = CCC_table[order(CCC_table$Aggr_score, decreasing=T), ]
# CCC_table$Aggr_score = CCC_table$Aggr_score + abs(min(CCC_table$Aggr_score)) + 1

x = CCC_table
x$Sender = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
x$Receiver = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
mtx1 = acast(x, Sender~Receiver, value.var="Aggr_score") # convert to nonsymmetric matrix
fords = c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "Pro GABAs", "VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "Peri 1", "Peri 2", "Fib", "PCs")
reindex = match(fords, colnames(mtx1))
mtx1 = mtx1[reindex, reindex]
write.csv(mtx1, "CCM_aggregated_CCC_matrix.csv", quote = F, col.names=T, row.names=T)






CCTObj = readRDS("CellCrossTalker_UHC_CCTObj.RDS")
totaldf = readRDS("CellCrossTalker_UHC_res.RDS")

x = c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "Pro GABAs", "VLMCs 1", "VLMCs 2", "VEs 1", "VEs 3", "Peri 1", "Peri 2", "SMCs","Fib", "PCs")
names(x) = levels(totaldf$Sender)

totaldf$Sender = as.character(totaldf$Sender)
totaldf$Receiver = as.character(totaldf$Receiver)

totaldf$Sender = x[match(totaldf$Sender, names(x))]
totaldf$Receiver = x[match(totaldf$Receiver, names(x))]

totaldf$cell_cell = paste(totaldf$Sender, totaldf$Receiver, sep = ">>>")
totaldf$LR = paste(totaldf$Ligand, totaldf$Receptor, sep = "---")
totaldf$LRscore = -log10(totaldf$KLscore + 0.0001) ## The smaller the KL divergence the greater connections
totaldf$LRscore = totaldf$LRscore + abs(min(totaldf$LRscore)) + 0.1
totaldf = arrange(totaldf, cell_cell, -LRscore)

saveRDS(totaldf, file = "CellCrossTalker_UHC_res.RDS")

CCC_table = aggregate(totaldf$LRscore, by = list(totaldf$cell_cell), FUN=sum)
names(CCC_table) = c("cell_cell", "Aggr_score")
CCC_table = CCC_table[order(CCC_table$Aggr_score, decreasing=T), ]
# CCC_table$Aggr_score = CCC_table$Aggr_score + abs(min(CCC_table$Aggr_score)) + 1

x = CCC_table
x$Sender = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
x$Receiver = sapply(x$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
mtx1 = acast(x, Sender~Receiver, value.var="Aggr_score") # convert to nonsymmetric matrix
fords = c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"Pro Ast", "Pro GABAs", "VLMCs 1", "VLMCs 2", "VEs 1", "VEs 3", "Peri 1", "Peri 2", "SMCs","Fib", "PCs")
reindex = match(fords, colnames(mtx1))
mtx1 = mtx1[reindex, reindex]
write.csv(mtx1, "UHC_aggregated_CCC_matrix.csv", quote = F, col.names=T, row.names=T)




#Thu Jan 26 12:05:46 2023

control = read.csv("Control_aggregated_CCC_matrix.csv", check.names=F, row.names=1)
CCM = read.csv("CCM_aggregated_CCC_matrix.csv", check.names=F, row.names=1)
UHC = read.csv("UHC_aggregated_CCC_matrix.csv", check.names=F, row.names=1)

a = c(as.numeric(as.matrix(control)), as.numeric(as.matrix(UHC)), as.numeric(as.matrix(CCM)))
control = as.matrix(control)
CCM = as.matrix(CCM)
UHC = as.matrix(UHC)

col_fun = colorRamp2(c(min(a)+median(a)/sd(a), median(a), max(a)-median(a)/sd(a)), c("#0000ff", "#fffafa", "#ff0000"))
pdf("Aggregated_CCC_matrix_three_groups.pdf", w = 10, h = 8)
Heatmap(control, name = "Control",col=col_fun,column_names_side = "top", rect_gp = gpar(col = "black", lwd = 2),
            
            column_names_gp = gpar(fontsize = 12,fontface="bold"),
            column_names_rot = -45,row_names_side = "left", row_names_gp = gpar(fontsize = 12,fontface="bold"),
            
            cluster_rows = FALSE,cluster_columns = FALSE,   heatmap_legend_param = list(
              title = "Score",grid_width=unit(0.6,"cm")
            ))
Heatmap(CCM, name = "CCM",col=col_fun,column_names_side = "top", rect_gp = gpar(col = "black", lwd = 2),
            
            column_names_gp = gpar(fontsize = 12,fontface="bold"),
            column_names_rot = -45,row_names_side = "left", row_names_gp = gpar(fontsize = 12,fontface="bold"),
            
            cluster_rows = FALSE,cluster_columns = FALSE,   heatmap_legend_param = list(
              title = "Score",grid_width=unit(0.6,"cm")
            ))
Heatmap(UHC, name = "UHC",col=col_fun,column_names_side = "top", rect_gp = gpar(col = "black", lwd = 2),
            
            column_names_gp = gpar(fontsize = 12,fontface="bold"),
            column_names_rot = -45,row_names_side = "left", row_names_gp = gpar(fontsize = 12,fontface="bold"),
            
            cluster_rows = FALSE,cluster_columns = FALSE,   heatmap_legend_param = list(
              title = "Score",grid_width=unit(0.6,"cm")
            ))
dev.off()




# step 2

control = readRDS("CellCrossTalker_Control_res.RDS")
CCM = readRDS("CellCrossTalker_CCM_res.RDS")
UHC = readRDS("CellCrossTalker_UHC_res.RDS")

neurons = c("Ast 1", "Ast 2", "Ast 3", "Ast 4", "NPs", "NPs 1", "NPs 2", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5")
vascular = c("VLMCs 1", "VLMCs 2", "VLMCs 3", "VEs 1", "VEs 2", "VEs 3", "Venous ECs 1", "Venous ECs 2", "Brain VEs", "Peri", "Peri 1", "Peri 2",
	"SMCs", "Fib")

# neurons2neurons
con1 = control[which(control$Sender %in% neurons & control$Receiver %in% neurons), 6:8]
CCM1 = CCM[which(CCM$Sender %in% neurons & CCM$Receiver %in% neurons), 6:8]
UHC1 = UHC[which(UHC$Sender %in% neurons & UHC$Receiver %in% neurons), 6:8]

m1 = list(UHC1, CCM1, con1) %>% Reduce(function(x, y) merge(x, y, by = c("cell_cell","LR"), all = T), .)
colnames(m1)[3:5] = c("UHC_LRscore", "CCM_LRscore", "control_LRscore")
m2 = m1[,3:5]; m2 = as.matrix(m2); m2[which(is.na(m2))] = 0; m2 = data.frame(m2)
m1 = cbind(m1[,1:2], m2)
m1$Sender = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
m1$Receiver = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
m1$Ligand = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[1])
m1$Receptor = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[2])

m2 = as.matrix(m1[,3:5])
m2 = scale(m2)
con_specific = which(m2[,3]>0 & m2[,1]<0 & m2[,2]<0)
uhc_specific = which(m2[,1]>0 & m2[,3]<0 & m2[,2]<0)
ccm_specific = which(m2[,2]>0 & m2[,1]<0 & m2[,3]<0)

m1$group_specific = ""
m1$group_specific[con_specific] = "control"
m1$group_specific[uhc_specific] = "UHC"
m1$group_specific[ccm_specific] = "CCM"

write.table(m1, "neurons2neurons_LR_specific.xls", quote=F, col.names=T, row.names=F, sep='\t')

library(ggsankey)
plotdf = m1[which(m1$group_specific %in% c("UHC", "CCM")), ]
plotdf = plotdf[,c("Ligand","Sender","Receiver","Receptor", "group_specific")]
plotdf$Ligand = "Sender"
plotdf$Receptor = plotdf$group_specific
plotdf = plotdf[,1:4]
plotdf1 = plotdf
plotdf1$Ligand = plotdf$Receptor
plotdf1$Receptor = plotdf$Ligand
plotdf1$Sender = plotdf$Receiver
plotdf1$Receiver = plotdf$Sender
plotdf = rbind(plotdf, plotdf1)
plotdf$Ligand = factor(plotdf$Ligand, levels = c("UHC", "Sender", "CCM"))
plotdf$Receptor = factor(plotdf$Receptor, levels = c("UHC", "Sender", "CCM"))
plotdf$Sender = factor(plotdf$Sender, levels = 
	c("Ast 1","Ast 2","Ast 3","NPs","GABAs","GluNs 1","GluNs 2","GluNs 3","GluNs 4", "GluNs 5")
	)
plotdf$Receiver = factor(plotdf$Receiver, levels = 
	c("Ast 1","Ast 2","Ast 3","NPs","GABAs","GluNs 1","GluNs 2","GluNs 3","GluNs 4", "GluNs 5")
	)
mycolor = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

pdf("neurons2neurons_sankey.pdf", w=12, h=8)
ggplot(plotdf, aes(x = Ligand, next_x = Receptor, node = Sender, next_node = Receiver, fill = Sender)) + 
	geom_sankey(flow.alpha = 0.5) + theme_sankey(base_size = 16) + scale_fill_manual(values = mycolor[1:12])
dev.off()                  




# neurons2vascular
con1 = control[which(control$Sender %in% neurons & control$Receiver %in% vascular), 6:8]
CCM1 = CCM[which(CCM$Sender %in% neurons & CCM$Receiver %in% vascular), 6:8]
UHC1 = UHC[which(UHC$Sender %in% neurons & UHC$Receiver %in% vascular), 6:8]

m1 = list(UHC1, CCM1, con1) %>% Reduce(function(x, y) merge(x, y, by = c("cell_cell","LR"), all = T), .)
colnames(m1)[3:5] = c("UHC_LRscore", "CCM_LRscore", "control_LRscore")
m2 = m1[,3:5]; m2 = as.matrix(m2); m2[which(is.na(m2))] = 0; m2 = data.frame(m2)
m1 = cbind(m1[,1:2], m2)
m1$Sender = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
m1$Receiver = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
m1$Ligand = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[1])
m1$Receptor = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[2])

m2 = as.matrix(m1[,3:5])
m2 = scale(m2)
con_specific = which(m2[,3]>0 & m2[,1]<0 & m2[,2]<0)
uhc_specific = which(m2[,1]>0 & m2[,3]<0 & m2[,2]<0)
ccm_specific = which(m2[,2]>0 & m2[,1]<0 & m2[,3]<0)

m1$group_specific = ""
m1$group_specific[con_specific] = "control"
m1$group_specific[uhc_specific] = "UHC"
m1$group_specific[ccm_specific] = "CCM"

write.table(m1, "neurons2vascular_LR_specific.xls", quote=F, col.names=T, row.names=F, sep='\t')

library(ggsankey)
plotdf = m1[which(m1$group_specific %in% c("UHC", "CCM")), ]
plotdf = plotdf[,c("Ligand","Sender","Receiver","Receptor", "group_specific")]
plotdf$Ligand = "Sender"
plotdf$Receptor = plotdf$group_specific
plotdf = plotdf[,1:4]
plotdf1 = plotdf
plotdf1$Ligand = plotdf$Receptor
plotdf1$Receptor = plotdf$Ligand
plotdf1$Sender = plotdf$Receiver
plotdf1$Receiver = plotdf$Sender
plotdf = rbind(plotdf, plotdf1)
plotdf$Ligand = factor(plotdf$Ligand, levels = c("UHC", "Sender", "CCM"))
plotdf$Receptor = factor(plotdf$Receptor, levels = c("UHC", "Sender", "CCM"))
plotdf$Sender = factor(plotdf$Sender, levels = 
	c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
plotdf$Receiver = factor(plotdf$Receiver, levels = 
	c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
mycolor = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

pdf("neurons2vascular_sankey.pdf", w=12, h=8)
ggplot(plotdf, aes(x = Ligand, next_x = Receptor, node = Sender, next_node = Receiver, fill = Sender)) + 
	geom_sankey(flow.alpha = 0.5) + theme_sankey(base_size = 16) + scale_fill_manual(values = mycolor[1:19])
dev.off()                  




# vascular2vascular
con1 = control[which(control$Sender %in% vascular & control$Receiver %in% vascular), 6:8]
CCM1 = CCM[which(CCM$Sender %in% vascular & CCM$Receiver %in% vascular), 6:8]
UHC1 = UHC[which(UHC$Sender %in% vascular & UHC$Receiver %in% vascular), 6:8]

m1 = list(UHC1, CCM1, con1) %>% Reduce(function(x, y) merge(x, y, by = c("cell_cell","LR"), all = T), .)
colnames(m1)[3:5] = c("UHC_LRscore", "CCM_LRscore", "control_LRscore")
m2 = m1[,3:5]; m2 = as.matrix(m2); m2[which(is.na(m2))] = 0; m2 = data.frame(m2)
m1 = cbind(m1[,1:2], m2)
m1$Sender = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
m1$Receiver = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
m1$Ligand = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[1])
m1$Receptor = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[2])

m2 = as.matrix(m1[,3:5])
m2 = scale(m2)
con_specific = which(m2[,3]>0 & m2[,1]<0 & m2[,2]<0)
uhc_specific = which(m2[,1]>0 & m2[,3]<0 & m2[,2]<0)
ccm_specific = which(m2[,2]>0 & m2[,1]<0 & m2[,3]<0)

m1$group_specific = ""
m1$group_specific[con_specific] = "control"
m1$group_specific[uhc_specific] = "UHC"
m1$group_specific[ccm_specific] = "CCM"

write.table(m1, "vascular2vascular_LR_specific.xls", quote=F, col.names=T, row.names=F, sep='\t')

library(ggsankey)
plotdf = m1[which(m1$group_specific %in% c("UHC", "CCM")), ]
plotdf = plotdf[,c("Ligand","Sender","Receiver","Receptor", "group_specific")]
plotdf$Ligand = "Sender"
plotdf$Receptor = plotdf$group_specific
plotdf = plotdf[,1:4]
plotdf1 = plotdf
plotdf1$Ligand = plotdf$Receptor
plotdf1$Receptor = plotdf$Ligand
plotdf1$Sender = plotdf$Receiver
plotdf1$Receiver = plotdf$Sender
plotdf = rbind(plotdf, plotdf1)
plotdf$Ligand = factor(plotdf$Ligand, levels = c("UHC", "Sender", "CCM"))
plotdf$Receptor = factor(plotdf$Receptor, levels = c("UHC", "Sender", "CCM"))
plotdf$Sender = factor(plotdf$Sender, levels = 
	c(
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
plotdf$Receiver = factor(plotdf$Receiver, levels = 
	c(
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
mycolor = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

pdf("vascular2vascular_sankey.pdf", w=12, h=8)
ggplot(plotdf, aes(x = Ligand, next_x = Receptor, node = Sender, next_node = Receiver, fill = Sender)) + 
	geom_sankey(flow.alpha = 0.5) + theme_sankey(base_size = 16) + scale_fill_manual(values = mycolor[11:22])
dev.off()                  





# vascular2neurons
con1 = control[which(control$Sender %in% vascular & control$Receiver %in% neurons), 6:8]
CCM1 = CCM[which(CCM$Sender %in% vascular & CCM$Receiver %in% neurons), 6:8]
UHC1 = UHC[which(UHC$Sender %in% vascular & UHC$Receiver %in% neurons), 6:8]

m1 = list(UHC1, CCM1, con1) %>% Reduce(function(x, y) merge(x, y, by = c("cell_cell","LR"), all = T), .)
colnames(m1)[3:5] = c("UHC_LRscore", "CCM_LRscore", "control_LRscore")
m2 = m1[,3:5]; m2 = as.matrix(m2); m2[which(is.na(m2))] = 0; m2 = data.frame(m2)
m1 = cbind(m1[,1:2], m2)
m1$Sender = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[1])
m1$Receiver = sapply(m1$cell_cell, function(i) unlist(strsplit(i, ">>>"))[2])
m1$Ligand = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[1])
m1$Receptor = sapply(m1$LR, function(i) unlist(strsplit(i, "---"))[2])

m2 = as.matrix(m1[,3:5])
m2 = scale(m2)
con_specific = which(m2[,3]>0 & m2[,1]<0 & m2[,2]<0)
uhc_specific = which(m2[,1]>0 & m2[,3]<0 & m2[,2]<0)
ccm_specific = which(m2[,2]>0 & m2[,1]<0 & m2[,3]<0)

m1$group_specific = ""
m1$group_specific[con_specific] = "control"
m1$group_specific[uhc_specific] = "UHC"
m1$group_specific[ccm_specific] = "CCM"

write.table(m1, "vascular2neurons_LR_specific.xls", quote=F, col.names=T, row.names=F, sep='\t')

library(ggsankey)
plotdf = m1[which(m1$group_specific %in% c("UHC", "CCM")), ]
plotdf = plotdf[,c("Ligand","Sender","Receiver","Receptor", "group_specific")]
plotdf$Ligand = "Sender"
plotdf$Receptor = plotdf$group_specific
plotdf = plotdf[,1:4]
plotdf1 = plotdf
plotdf1$Ligand = plotdf$Receptor
plotdf1$Receptor = plotdf$Ligand
plotdf1$Sender = plotdf$Receiver
plotdf1$Receiver = plotdf$Sender
plotdf = rbind(plotdf, plotdf1)
plotdf$Ligand = factor(plotdf$Ligand, levels = c("UHC", "Sender", "CCM"))
plotdf$Receptor = factor(plotdf$Receptor, levels = c("UHC", "Sender", "CCM"))
plotdf$Sender = factor(plotdf$Sender, levels = 
	c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
plotdf$Receiver = factor(plotdf$Receiver, levels = 
	c("Ast 1", "Ast 2", "Ast 3", "NPs", "GABAs", "GluNs 1", "GluNs 2", "GluNs 3", "GluNs 4", "GluNs 5",
	"VLMCs 1", "VLMCs 2", "VEs 1", "VEs 2", "VEs 3", "Peri 1", "Peri 2", "SMCs", "Fib")
	)
mycolor = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

pdf("vascular2neurons_sankey.pdf", w=12, h=8)
ggplot(plotdf, aes(x = Ligand, next_x = Receptor, node = Sender, next_node = Receiver, fill = Sender)) + 
	geom_sankey(flow.alpha = 0.5) + theme_sankey(base_size = 16) + scale_fill_manual(values = mycolor[1:19])
dev.off()                  






####
#Thu Jan 26 17:13:26 2023
library(reshape2)
m1 = read.table("neurons2neurons_LR_specific.xls", header=T, sep="\t")
m1$DIFF = m1$UHC_LRscore - m1$CCM_LRscore
m1 = m1[order(m1$DIFF), ]
m1$ID = paste0(m1$Sender,":",m1$Ligand,">>",m1$Receiver,":",m1$Receptor)

# CCM enriched
CCMe = m1[1:10,-c(6:11)]
UHCe = m1[(nrow(m1)-9):nrow(m1),-c(6:11)]
CCMe = melt(CCMe, id = c("cell_cell", "LR", "ID"))
UHCe = melt(UHCe, id = c("cell_cell", "LR", "ID"))
CCMe$variable = gsub("_LRscore", "", CCMe$variable)
UHCe$variable = gsub("_LRscore", "", UHCe$variable)
CCMe$variable = factor(CCMe$variable, levels = c("UHC", "CCM", "control"))
UHCe$variable = factor(UHCe$variable, levels = c("UHC", "CCM", "control"))


pdf("neurons2neurons_specific_bubble.pdf", w = 6, h=5)
ggplot(data = CCMe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
ggplot(data = UHCe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
dev.off()




m1 = read.table("neurons2vascular_LR_specific.xls", header=T, sep="\t")
m1$DIFF = m1$UHC_LRscore - m1$CCM_LRscore
m1 = m1[order(m1$DIFF), ]
m1$ID = paste0(m1$Sender,":",m1$Ligand,">>",m1$Receiver,":",m1$Receptor)

# CCM enriched
CCMe = m1[1:10,-c(6:11)]
UHCe = m1[(nrow(m1)-9):nrow(m1),-c(6:11)]
CCMe = melt(CCMe, id = c("cell_cell", "LR", "ID"))
UHCe = melt(UHCe, id = c("cell_cell", "LR", "ID"))
CCMe$variable = gsub("_LRscore", "", CCMe$variable)
UHCe$variable = gsub("_LRscore", "", UHCe$variable)
CCMe$variable = factor(CCMe$variable, levels = c("UHC", "CCM", "control"))
UHCe$variable = factor(UHCe$variable, levels = c("UHC", "CCM", "control"))


pdf("neurons2vascular_specific_bubble.pdf", w = 6, h=5)
ggplot(data = CCMe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
ggplot(data = UHCe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
dev.off()




m1 = read.table("vascular2vascular_LR_specific.xls", header=T, sep="\t")
m1$DIFF = m1$UHC_LRscore - m1$CCM_LRscore
m1 = m1[order(m1$DIFF), ]
m1$ID = paste0(m1$Sender,":",m1$Ligand,">>",m1$Receiver,":",m1$Receptor)

# CCM enriched
CCMe = m1[1:10,-c(6:11)]
UHCe = m1[(nrow(m1)-9):nrow(m1),-c(6:11)]
CCMe = melt(CCMe, id = c("cell_cell", "LR", "ID"))
UHCe = melt(UHCe, id = c("cell_cell", "LR", "ID"))
CCMe$variable = gsub("_LRscore", "", CCMe$variable)
UHCe$variable = gsub("_LRscore", "", UHCe$variable)
CCMe$variable = factor(CCMe$variable, levels = c("UHC", "CCM", "control"))
UHCe$variable = factor(UHCe$variable, levels = c("UHC", "CCM", "control"))


pdf("vascular2vascular_specific_bubble.pdf", w = 6, h=5)
ggplot(data = CCMe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
ggplot(data = UHCe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
dev.off()




m1 = read.table("vascular2neurons_LR_specific.xls", header=T, sep="\t")
m1$DIFF = m1$UHC_LRscore - m1$CCM_LRscore
m1 = m1[order(m1$DIFF), ]
m1$ID = paste0(m1$Sender,":",m1$Ligand,">>",m1$Receiver,":",m1$Receptor)

# CCM enriched
CCMe = m1[1:10,-c(6:11)]
UHCe = m1[(nrow(m1)-9):nrow(m1),-c(6:11)]
CCMe = melt(CCMe, id = c("cell_cell", "LR", "ID"))
UHCe = melt(UHCe, id = c("cell_cell", "LR", "ID"))
CCMe$variable = gsub("_LRscore", "", CCMe$variable)
UHCe$variable = gsub("_LRscore", "", UHCe$variable)
CCMe$variable = factor(CCMe$variable, levels = c("UHC", "CCM", "control"))
UHCe$variable = factor(UHCe$variable, levels = c("UHC", "CCM", "control"))


pdf("vascular2neurons_specific_bubble.pdf", w = 6, h=5)
ggplot(data = CCMe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
ggplot(data = UHCe, mapping = aes_string(x = "variable", 
                                             y = "ID")) + 
  geom_point(mapping = aes_string(size="value",color="value"))+
  scale_color_gradient2(low="grey",mid="blue",high="red",midpoint=2,limits=c(0,6.5))+
  labs(title="")+
  xlab("")+ylab("")+theme_bw()+scale_size_continuous(limits=c(0,8))+
  theme(legend.title = element_text(face="bold",size=12),
        plot.title=element_text(face="bold",size=28,hjust=0.5),
        axis.text.x = element_text(angle=90,vjust=1,hjust=1,size=14,face ="bold", color="black"),
        axis.text.y = element_text(size=14,face ="bold", color = "black")
        # axis.title.x = element_text(size=20),
        # axis.title.y = element_text(size=20)
  )
dev.off()
