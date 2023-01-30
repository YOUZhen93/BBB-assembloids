# Stream
Tue Jan 24 09:23:58 2023
options(stringsAsFactors=F)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

pbmc = readRDS("CCM_UHC_merged_SeuratObj.RDS")

# prepare input

mycolor = c('#689C6F','#547DB4','#E8B0A9','#874F9F','#D9C541','#C3718A','#F0D33D','#728D68','#5E8D92','#DD7E1E','#C36E3C','#6A5C81','#DC80BC','#99BBCF','#FEFF4E','#A55E6B','#E5A72D','#C32021','#93572F','#B58D37','#AB645B','#7A6D83')

names(mycolor) = levels(Idents(pbmc))

neuron = pbmc[, Idents(pbmc) %in% levels(Idents(pbmc))[1:12]]
vascula = pbmc[, Idents(pbmc) %in% levels(Idents(pbmc))[13:22]]

meta = neuron@meta.data[,c("group", "Sub_celltype")]
meta$label_color = mycolor[match(meta$Sub_celltype, names(mycolor))]
meta1 = meta[which(meta$group == "CCM"), ]
meta2 = meta[which(meta$group == "UHC"), ]
meta1 = meta1[,-1]; meta2 = meta2[,-1]
colnames(meta1)[1] = colnames(meta2)[1] = "label"
rownames(meta1) = gsub("-", "\\.", rownames(meta1))
rownames(meta2) = gsub("-", "\\.", rownames(meta2))
write.table(meta1, "CCM_neurons_metadata.tsv", sep="\t", quote = F, col.names=T, row.names=T)
write.table(meta2, "UHC_neurons_metadata.tsv", sep="\t", quote = F, col.names=T, row.names=T)

meta = vascula@meta.data[,c("group", "Sub_celltype")]
meta$label_color = mycolor[match(meta$Sub_celltype, names(mycolor))]
meta1 = meta[which(meta$group == "CCM"), ]
meta2 = meta[which(meta$group == "UHC"), ]
meta1 = meta1[,-1]; meta2 = meta2[,-1]
colnames(meta1)[1] = colnames(meta2)[1] = "label"
rownames(meta1) = gsub("-", "\\.", rownames(meta1))
rownames(meta2) = gsub("-", "\\.", rownames(meta2))
write.table(meta1, "CCM_vascular_metadata.tsv", sep="\t", quote = F, col.names=T, row.names=T)
write.table(meta2, "UHC_vascular_metadata.tsv", sep="\t", quote = F, col.names=T, row.names=T)


nc = neuron@assays$SCT@counts
nc1 = nc[,which(neuron$group == "CCM")]
nc2 = nc[,which(neuron$group == "UHC")]
write.table(nc1, "CCM_neurons_RawCounts.tsv", sep="\t", quote=F, col.names=T, row.names=T)
write.table(nc2, "UHC_neurons_RawCounts.tsv", sep="\t", quote=F, col.names=T, row.names=T)


nc = vascula@assays$SCT@counts
nc1 = nc[,which(vascula$group == "CCM")]
nc2 = nc[,which(vascula$group == "UHC")]
write.table(nc1, "CCM_vascular_RawCounts.tsv", sep="\t", quote=F, col.names=T, row.names=T)
write.table(nc2, "UHC_vascular_RawCounts.tsv", sep="\t", quote=F, col.names=T, row.names=T)





# running stream
import stream as st
st.set_figure_params(dpi=300,style='white',figsize=[5.4,4.8],
	rc={'image.cmap': 'viridis'})

adata=st.read(file_name='./CCM_vascular_RawCounts.tsv.gz',workdir='/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/stream')
adata.var_names_make_unique()

adata

st.add_metadata(adata,file_name='./CCM_vascular_metadata.tsv')

adata.obs.head()

st.cal_qc(adata,assay='rna')
st.plot_qc(adata,jitter=0.3,)
st.plot_qc(adata,jitter=0.3,log_scale=[0,1,4,5],hist_plot=[0,1,4,5]) 
st.filter_cells(adata,min_n_features= 0)
st.filter_features(adata,min_n_cells = 0)

st.normalize(adata,method='lib_size')
st.log_transform(adata)
st.remove_mt_genes(adata)
st.select_variable_genes(adata,loess_frac=0.01,n_genes=500)
st.dimension_reduction(adata,method='mlle',feature='var_genes',n_components=10,nb_pct=0.025,n_jobs=8)
st.plot_visualization_2D(adata,nb_pct=0.1,color=['label'],use_precomputed=False)

adata_low = st.switch_to_low_dimension(adata,n_components=3)
st.plot_dimension_reduction(adata_low,n_components=3)

st.seed_elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=2,show_graph=True,show_text=False)
st.elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True)

st.plot_stream_sc(adata_low,root='S0',color=['label'], dist_scale=0.2,show_graph=True,show_text=True)

init_nodes_pos,init_edges = st.infer_initial_structure(adata_low)
st.seed_elastic_principal_graph(adata,init_nodes_pos=init_nodes_pos,init_edges=init_edges)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'],fig_ncol=4)
st.elastic_principal_graph(adata,incr_n_nodes=10)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'])

st.extend_elastic_principal_graph(adata)
st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=True)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=False,show_text=False)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=True,show_text=True)

st.plot_stream_sc(adata,root='S0',dist_scale=0.5)
st.plot_stream(adata,root='S0')
st.plot_stream(adata,root='S0',dist_scale=1.2,factor_min_win=1.2)
st.plot_stream(adata,root='S0',log_scale=True,factor_min_win=1.2)

st.write(adata,file_name='CCM_vascular_STREAM.pkl')
st.write(adata_low,file_name='CCM_vascular_STREAM_low.pkl')




# running stream
import stream as st
st.set_figure_params(dpi=300,style='white',figsize=[5.4,4.8],
	rc={'image.cmap': 'viridis'})

adata=st.read(file_name='./CCM_neurons_RawCounts.tsv.gz',workdir='/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/stream')
adata.var_names_make_unique()

adata

st.add_metadata(adata,file_name='./CCM_neurons_metadata.tsv')

adata.obs.head()

st.cal_qc(adata,assay='rna')
st.plot_qc(adata,jitter=0.3,)
st.plot_qc(adata,jitter=0.3,log_scale=[0,1,4,5],hist_plot=[0,1,4,5]) 
st.filter_cells(adata,min_n_features= 0)
st.filter_features(adata,min_n_cells = 0)

st.normalize(adata,method='lib_size')
st.log_transform(adata)
st.remove_mt_genes(adata)
st.select_variable_genes(adata,loess_frac=0.01,n_genes=500)
#st.dimension_reduction(adata,method='mlsle',feature='var_genes',n_components=10,nb_pct=0.025,n_jobs=10)
st.dimension_reduction(adata,method='se',feature='var_genes', n_neighbors=50,n_components=10,n_jobs=10)
st.plot_visualization_2D(adata,nb_pct=0.1,color=['label'],use_precomputed=False)

adata_low = st.switch_to_low_dimension(adata,n_components=3)
st.plot_dimension_reduction(adata_low,n_components=3)

st.seed_elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=2,show_graph=True,show_text=False)
st.elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True)

st.plot_stream_sc(adata_low,root='S0',color=['label'], dist_scale=0.2,show_graph=True,show_text=True)

init_nodes_pos,init_edges = st.infer_initial_structure(adata_low)
st.seed_elastic_principal_graph(adata,init_nodes_pos=init_nodes_pos,init_edges=init_edges)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'],fig_ncol=4)
st.elastic_principal_graph(adata,incr_n_nodes=10)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'])

st.extend_elastic_principal_graph(adata)
st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=True)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=False,show_text=False)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=True,show_text=True)

st.plot_stream_sc(adata,root='S0',dist_scale=0.5)
st.plot_stream(adata,root='S0')
st.plot_stream(adata,root='S0',dist_scale=1.2,factor_min_win=1.2)
st.plot_stream(adata,root='S0',log_scale=True,factor_min_win=1.2)

st.write(adata,file_name='CCM_vascular_STREAM.pkl')
st.write(adata_low,file_name='CCM_vascular_STREAM_low.pkl')






# running stream
import stream as st
st.set_figure_params(dpi=300,style='white',figsize=[5.4,4.8],
	rc={'image.cmap': 'viridis'})

adata=st.read(file_name='./UHC_vascular_RawCounts.tsv.gz',workdir='/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/stream')
adata.var_names_make_unique()

adata

st.add_metadata(adata,file_name='./UHC_vascular_metadata.tsv')

adata.obs.head()

st.cal_qc(adata,assay='rna')
st.plot_qc(adata,jitter=0.3,)
st.plot_qc(adata,jitter=0.3,log_scale=[0,1,4,5],hist_plot=[0,1,4,5]) 
st.filter_cells(adata,min_n_features= 0)
st.filter_features(adata,min_n_cells = 0)

st.normalize(adata,method='lib_size')
st.log_transform(adata)
st.remove_mt_genes(adata)
st.select_variable_genes(adata,loess_frac=0.01,n_genes=500)
st.dimension_reduction(adata,method='mlle',feature='var_genes',n_components=10,nb_pct=0.025,n_jobs=8)
st.plot_visualization_2D(adata,nb_pct=0.1,color=['label'],use_precomputed=False)

adata_low = st.switch_to_low_dimension(adata,n_components=3)
st.plot_dimension_reduction(adata_low,n_components=3)

st.seed_elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=2,show_graph=True,show_text=False)
st.elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True)

st.plot_stream_sc(adata_low,root='S2',color=['label'], dist_scale=0.2,show_graph=True,show_text=True)

init_nodes_pos,init_edges = st.infer_initial_structure(adata_low)
st.seed_elastic_principal_graph(adata,init_nodes_pos=init_nodes_pos,init_edges=init_edges)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'],fig_ncol=4)
st.elastic_principal_graph(adata,incr_n_nodes=10)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'])

st.extend_elastic_principal_graph(adata)
st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=True)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=False,show_text=False)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=True,show_text=True)

st.plot_stream_sc(adata,root='S6',dist_scale=0.5)
st.plot_stream(adata,root='S6')
st.plot_stream(adata,root='S6',dist_scale=1.2,factor_min_win=1.2)
st.plot_stream(adata,root='S6',log_scale=True,factor_min_win=1.2)

st.write(adata,file_name='UHC_vascular_STREAM.pkl')
st.write(adata_low,file_name='UHC_vascular_STREAM_low.pkl')






# running stream
import stream as st
st.set_figure_params(dpi=300,style='white',figsize=[5.4,4.8],
	rc={'image.cmap': 'viridis'})

adata=st.read(file_name='./UHC_neurons_RawCounts.tsv.gz',workdir='/public/home/gonglianggroup/zheny/temp/Brain_Organoid/Figures/stream')
adata.var_names_make_unique()

adata

st.add_metadata(adata,file_name='./UHC_neurons_metadata.tsv')

adata.obs.head()

st.cal_qc(adata,assay='rna')
st.plot_qc(adata,jitter=0.3,)
st.plot_qc(adata,jitter=0.3,log_scale=[0,1,4,5],hist_plot=[0,1,4,5]) 
st.filter_cells(adata,min_n_features= 0)
st.filter_features(adata,min_n_cells = 0)

st.normalize(adata,method='lib_size')
st.log_transform(adata)
st.remove_mt_genes(adata)
st.select_variable_genes(adata,loess_frac=0.01,n_genes=500)
st.dimension_reduction(adata,method='mlle',feature='var_genes',n_components=10,nb_pct=0.025,n_jobs=8)
st.plot_visualization_2D(adata,nb_pct=0.1,color=['label'],use_precomputed=False)

adata_low = st.switch_to_low_dimension(adata,n_components=3)
st.plot_dimension_reduction(adata_low,n_components=3)

st.seed_elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=2,show_graph=True,show_text=False)
st.elastic_principal_graph(adata_low)
st.plot_dimension_reduction(adata_low,color=['label'],n_components=3,show_graph=True,show_text=True)

st.plot_stream_sc(adata_low,root='S0',color=['label'], dist_scale=0.2,show_graph=True,show_text=True)

init_nodes_pos,init_edges = st.infer_initial_structure(adata_low)
st.seed_elastic_principal_graph(adata,init_nodes_pos=init_nodes_pos,init_edges=init_edges)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'],fig_ncol=4)
st.elastic_principal_graph(adata,incr_n_nodes=10)
st.plot_visualization_2D(adata,color=['label','branch_id_alias'])

st.extend_elastic_principal_graph(adata)
st.plot_dimension_reduction(adata,color=['label'],n_components=3,show_graph=True,show_text=True)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=False,show_text=False)
st.plot_flat_tree(adata,color=['label','branch_id_alias','S0_pseudotime'], dist_scale=0.3,show_graph=True,show_text=True)

st.plot_stream_sc(adata,root='S0',dist_scale=0.5)
st.plot_stream(adata,root='S0')
st.plot_stream(adata,root='S0',dist_scale=1.2,factor_min_win=1.2)
st.plot_stream(adata,root='S0',log_scale=True,factor_min_win=1.2)

st.write(adata,file_name='UHC_neurons_STREAM.pkl')
st.write(adata_low,file_name='UHC_neurons_STREAM_low.pkl')




