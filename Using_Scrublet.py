#####11/14/2022 10:30:59 AM
# Scrublet find doublet

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import pandas as pd

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42


inpd='.'

def doScrublet(inpd, output):
	counts_matrix = scipy.io.mmread(inpd + '/matrix.mtx').T.tocsc()
	genes = np.array(scr.load_genes(inpd + '/features.tsv', delimiter='\t', column=2))
	print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
	print('Number of genes in gene list: {}'.format(len(genes)))


	# iteration of different expected doublet rate if this is unknown 0.05-0.1
	scrub5 = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)
	scrub6 = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
	scrub8 = scr.Scrublet(counts_matrix, expected_doublet_rate=0.08)
	scrub1 = scr.Scrublet(counts_matrix, expected_doublet_rate=0.1)

	doublet_scores5, predicted_doublets5 = scrub5.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
	doublet_scores6, predicted_doublets6 = scrub6.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
	doublet_scores8, predicted_doublets8 = scrub8.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
	doublet_scores1, predicted_doublets1 = scrub1.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
	
	scrub5.plot_histogram();
	plt.savefig(inpd + '/histo5.pdf')
	scrub6.plot_histogram();
	plt.savefig(inpd + '/histo6.pdf')
	scrub8.plot_histogram();
	plt.savefig(inpd + '/histo8.pdf')
	scrub1.plot_histogram();
	plt.savefig(inpd + '/histo1.pdf')

	print('Running UMAP with par 0.05...')
	scrub5.set_embedding('UMAP', scr.get_umap(scrub5.manifold_obs_, 10, min_dist=0.3))
	print('Done.')
	print('Print detected doublet rate:')
	print(scrub5.detected_doublet_rate_)
	scrub5.plot_embedding('UMAP', order_points=True);
	plt.savefig(inpd + '/UMAP5.pdf')

	print('Running UMAP with par 0.06...')
	scrub6.set_embedding('UMAP', scr.get_umap(scrub6.manifold_obs_, 10, min_dist=0.3))
	print('Done.')
	print('Print detected doublet rate:')
	print(scrub6.detected_doublet_rate_)
	scrub6.plot_embedding('UMAP', order_points=True);
	plt.savefig(inpd + '/UMAP6.pdf')

	print('Running UMAP with par 0.08...')
	scrub8.set_embedding('UMAP', scr.get_umap(scrub8.manifold_obs_, 10, min_dist=0.3))
	print('Done.')
	print('Print detected doublet rate:')
	print(scrub8.detected_doublet_rate_)
	scrub8.plot_embedding('UMAP', order_points=True);
	plt.savefig(inpd + '/UMAP8.pdf')

	print('Running UMAP with par 0.1...')
	scrub1.set_embedding('UMAP', scr.get_umap(scrub1.manifold_obs_, 10, min_dist=0.3))
	print('Done.')
	print('Print detected doublet rate:')
	print(scrub1.detected_doublet_rate_)
	scrub1.plot_embedding('UMAP', order_points=True);
	plt.savefig(inpd + '/UMAP1.pdf')

	out_df=pd.read_csv(inpd + "/barcodes.tsv", header = None, index_col = None)
	out_df['doublet_scores0.05'] = doublet_scores5
	out_df['doublet_scores0.06'] = doublet_scores6
	out_df['doublet_scores0.08'] = doublet_scores8
	out_df['doublet_scores0.1'] = doublet_scores1
	out_df['predicted_doublets0.05'] = predicted_doublets5
	out_df['predicted_doublets0.07'] = predicted_doublets8
	out_df['predicted_doublets0.08'] = predicted_doublets8
	out_df['predicted_doublets0.1'] = predicted_doublets1
	out_df.to_csv(inpd + output, index=False, header=True, sep="\t")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='find doublet with scublet')
    parser.add_argument('--input', metavar='./', required=True,  
    	help='<input dir>')
    parser.add_argument('--output', metavar='barcodes.tsv', required=True,  
    	help='<output barcodes.tsv>')
    args = parser.parse_args()

    doScrublet(args.input, args.output)


