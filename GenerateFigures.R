# making related Figures
# Author Zheny
options(stringsAsFactors = F)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


1. UMAP
# CCM UHC UMAP:
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")

ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mycolors = sample(mycolors)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("#276893","#c65306", "#dbce54", "#5a5c5b", "#c35655", "#2ec3e7"), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("red", "navy"), group.by = 'group') + 
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))


# Control UMAP:
pbmc = readRDS("Control_merged_SeuratObj.RDS")
ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)
mycolors = sample(mycolors)
DimPlot(pbmc, label = TRUE, cols = mycolors, label.size = 6, repel = T) +
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))
DimPlot(pbmc, cols = c("#276893","#c65306", "#dbce54"), group.by = 'batch') +  
theme(axis.text = element_text(size=18), axis.title = element_text(size=20), axis.line = element_line(size=1), axis.ticks = element_line(size=1))



2. Violin Plot

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
# CCM UHC
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
mycolors = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

mymarker = c('SLC1A3', 'HES5', 'VEGFA', 'BNIP3', 'CXCL14', 'CCNB2','GAD2','DLX1', 'GABBR2', 'GRIN2B', 'BHLHE22', 'NEUROD2', 'EOMES', 'DDIT4', 'LUM', 'COL3A1', 'TSHZ2', 
	'PDGFRB', 'RGS5', 'REN', 'THBS4', 'SLC2A1', 'CLDN5', 'CD34', 'KDR', 'PRND', 'ACTA2', 'TAGLN', 'HIST1H1B', 'PTTG1', 'TOP2A')
StackedVlnPlot(obj = pbmc, features = mymarker[1:16], cols = mycolors)
StackedVlnPlot(obj = pbmc, features = mymarker[17:31], cols = mycolors)

# Control
pbmc = readRDS("Control_merged_SeuratObj.RDS")
ncols = length(unique(Idents(pbmc)))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(ncols)

mymarker = c('SLC1A3', 'NOTCH1', 'SOX2', 'NES', 'TOP2A', 'CENPF','MKI67', 'GAD1', 'GAD2', 'DLX2', 'SLC17A6', 'SLC17A7', 'TBR1', 'GABBR2', 'GRIN2B', 'GABBR1', 
	'LUM', 'COL3A1', 'COL1A1', 'ALDOA', 'SLC2A1', 'FN1', 'RSPO3', 'IGF1', 'MEOX2', 'MECOM', 'TNMD', 'THBS4', 'CLDN5', 'ESAM')
StackedVlnPlot(obj = pbmc.combined, features = mymarker[1:15], cols = mycolors)
StackedVlnPlot(obj = pbmc.combined, features = mymarker[16:30], cols = mycolors)


3. cell proporation

### cell type proportions
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
df = pbmc@meta.data
df = df[,c('Sub_celltype', 'group')]
plotdf = as.data.frame(df %>% count(Sub_celltype, group, name="Count"))

plotdf$Count[which(plotdf$group == "CCM")] = plotdf$Count[which(plotdf$group == "CCM")]/20870
plotdf$Count[which(plotdf$group == "UHC")] = plotdf$Count[which(plotdf$group == "UHC")]/21729
plotdf$Count = round(plotdf$Count * 100, digit = 3)

plotdf$Sub_celltype = factor(plotdf$Sub_celltype, levels = levels(Idents(pbmc)))

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



4. heatmap

# same for control and CCM_UHC 
pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")
markers = readRDS("CCM_UHC_merged_celltype_markers.RDS")
markers = markers[which(markers$avg_log2FC > 0), ]
markers %>% group_by(cluster)  %>% arrange(p_val_adj, desc(avg_log2FC)) -> newmarkers
newmarkers = as.data.frame(newmarkers)
newmarkers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10

DoHeatmap(pbmc, features = top10$gene, raster = T) + NoLegend() + scale_fill_gradientn(colors = c("#00ffff", "black", "#ffff00"))





5. cell-cell aggregated communications
options(stringsAsFactors=F)
library(tidyverse)
library(reshape2)
library(circlize)
library(ComplexHeatmap)

control = read.csv("Control_aggregated_CCC_matrix.csv", check.names=F, row.names=1)
CCM = read.csv("CCM_aggregated_CCC_matrix.csv", check.names=F, row.names=1)
UHC = read.csv("UHC_aggregated_CCC_matrix.csv", check.names=F, row.names=1)

a = c(as.numeric(as.matrix(control)), as.numeric(as.matrix(UHC)), as.numeric(as.matrix(CCM)))
control = as.matrix(control)
CCM = as.matrix(CCM)
UHC = as.matrix(UHC)

col_fun = colorRamp2(c(min(a)+median(a)/sd(a), median(a), max(a)-median(a)/sd(a)), c("#0000ff", "#fffafa", "#ff0000"))

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



5. cell-cell specifc communications

control = read.table("Control_CCC_results.xls", header=T, sep="\t")
CCM = readRDS("CCM_CCC_results.xls", header=T, sep="\t")
UHC = readRDS("UHC_aggregated_CCC_matrix.csv", header=T, sep="\t")

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

write.table(m1, "neurons2neurons_LR_specific.xls", quote=F, col.names=T, row.names=F, sep='\t')

m2 = as.matrix(m1[,3:5])
m2 = scale(m2)
con_specific = which(m2[,3]>0 & m2[,1]<0 & m2[,2]<0)
uhc_specific = which(m2[,1]>0 & m2[,3]<0 & m2[,2]<0)
ccm_specific = which(m2[,2]>0 & m2[,1]<0 & m2[,3]<0)

m1$group_specific = ""
m1$group_specific[con_specific] = "control"
m1$group_specific[uhc_specific] = "UHC"
m1$group_specific[ccm_specific] = "CCM"
#neuron to neuron
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

ggplot(plotdf, aes(x = Ligand, next_x = Receptor, node = Sender, next_node = Receiver, fill = Sender)) + 
	geom_sankey(flow.alpha = 0.5) + theme_sankey(base_size = 16) + scale_fill_manual(values = mycolor[1:12])




6. cell-cell top 10 LR pairs

library(reshape2)
m1 = read.table("neurons2neurons_LR_specific.xls", header=T, sep="\t")
m1$DIFF = m1$UHC_LRscore - m1$CCM_LRscore
m1 = m1[order(m1$DIFF), ]
m1$ID = paste0(m1$Sender,":",m1$Ligand,">>",m1$Receiver,":",m1$Receptor)

CCMe = m1[1:10,-c(6:11)]
UHCe = m1[(nrow(m1)-9):nrow(m1),-c(6:11)]
CCMe = melt(CCMe, id = c("cell_cell", "LR", "ID"))
UHCe = melt(UHCe, id = c("cell_cell", "LR", "ID"))
CCMe$variable = gsub("_LRscore", "", CCMe$variable)
UHCe$variable = gsub("_LRscore", "", UHCe$variable)
CCMe$variable = factor(CCMe$variable, levels = c("UHC", "CCM", "control"))
UHCe$variable = factor(UHCe$variable, levels = c("UHC", "CCM", "control"))

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
  )







