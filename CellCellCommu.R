# cellcrosstalker
# Author Zhen Y

# control group
# step 1. load Seurat object and cell type DEGs list
require(pbmcapply)
pbmc = readRDS("Control_merged_SeuratObj.RDS")
markers = readRDS("Control_merged_celltype_markers.RDS")

# step 2 source DGBD function
source("DGBD_ready2use.R")
countTable = pbmc@assays$RNA@counts

# step 3 load LR interactions table and normalization
LR = read.csv("mouse_human_lrdb_version2.csv")
countTable = as.matrix(countTable[rownames(countTable) %in% unique(c(LR$Ligand.Human, LR$Receptor.Human)), ])
# normalized
cTMM = TMMnormalization(countTable)

# step 4 create cell to cell dataframe
celltypes = levels(unique(Idents(pbmc)))
celltype1 = ""; celltype2 = ""
condition = as.character(Idents(pbmc));



#  make all combination df
workdf = expand.grid(celltypes, celltypes)
names(workdf) = c("Sender", "Receiver")

# step 5 calculate cell cell insteractions score
totaldf = data.frame()
CCTObj = list()
for (i in seq(nrow(workdf))) {
    print(paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> "))
    tempdf = CellCT(cTMM, workdf$Sender[i], workdf$Receiver[i], 0.95, condition, 1)
    tdf1 = tempdf@result_df; tdf1$Sender = workdf$Sender[i]; tdf1$Receiver = workdf$Receiver[i];
    totaldf = rbind(totaldf, tdf1)
    CCTObj[[paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> ")]] = c(tempdf@paramL, tempdf@paramR)
}







# CCM UHC groups
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
markers = readRDS("CCM_UHC_merged_celltype_markers.RDS")
markers = markers[which(markers$avg_log2FC > 0), ]
markers1 = markers
source("DGBD_ready2use.R")
# normalized
LR = read.csv("mouse_human_lrdb_version2.csv")

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
	    tempdf = CellCT(cTMM, workdf$Sender[i], workdf$Receiver[i], 0.95, condition, 1)
	    tdf1 = tempdf@result_df; tdf1$Sender = workdf$Sender[i]; tdf1$Receiver = workdf$Receiver[i];
	    totaldf = rbind(totaldf, tdf1)
	    CCTObj[[paste(workdf$Sender[i], workdf$Receiver[i], sep=" >>> ")]] = c(tempdf@paramL, tempdf@paramR)
	}
	saveRDS(totaldf, file = paste0("CellCrossTalker_", group, "_res.RDS"))
	saveRDS(CCTObj, file = paste0("CellCrossTalker_", group, "_CCTObj.RDS"))
}




