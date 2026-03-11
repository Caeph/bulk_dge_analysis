import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from dge_visualization import *

import re
import os
inference = DefaultInference(n_cpus=8)

def quality_control(count_table_df,
                    samples_df,
                    reporter,
                    figure_directory,
                    pseudocount = 0.01,
                    ntop=500):
    os.makedirs(figure_directory, exist_ok=True)

    print("Performing quality control...")
    reporter.start_section("Quality control")

    log_counts = np.log2(count_table_df + pseudocount)
    log_counts_annot = count_table_df.index

    fig = plot_clustering(log_counts, log_counts_annot)
    reporter.send_figure(
        fig,
        "Hierarchical clustering on correlation between replicates",
        os.path.join(figure_directory, "hierarchical_clustering.png")
    )
    plt.close(fig)


    # PCA
    fig = plot_pca(count_table_df,
                   log_counts,
                   log_counts_annot,
                   samples_df,
                   ntop)
    reporter.send_figure(
        fig,
        f"PCA on {ntop} most variable genes on all replicates\n(both plots whos the same information)",
        os.path.join(figure_directory, "pca.png")
    )
    plt.close(fig)

    print("--")

def load_gene_annotation(gtf_path):
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#', header=None, low_memory=False)
    attribute_from, attribute_to = 'gene_id', 'gene_name'
    gtf['from'] = gtf[8].str.extract(rf'{re.escape(attribute_from)}\s+"([^"]+)"')
    gtf['to'] = gtf[8].str.extract(rf'{re.escape(attribute_to)}\s+"([^"]+)"')
    gtf['to2'] = gtf[8].str.extract(rf'{re.escape("gene")}\s+"([^"]+)"')
    gene_names = pd.concat(
        [gtf[['from', 'to']].drop_duplicates().dropna().rename(columns={"from": "gene_id", "to": "gene_name"}),
         gtf[['from', 'to2']].drop_duplicates().dropna().rename(columns={"from": "gene_id", "to2": "gene_name"})]
    )
    gene_names['organism'] = np.select([~gene_names['gene_id'].str.startswith("ENSG")], ['virus'], default="human")
    gene_names['gene_name'] = gene_names['organism'] + '-' + gene_names['gene_name']
    gene_names.index = gene_names['gene_id']
    gene_names = gene_names.drop(columns=['organism'])
    return gene_names

def combine_to_input(count_table_df,
                     samples_df,
                     minimal_gene_count):
    samples_df.index = samples_df['sample']
    samples_df.sort_index(inplace=True)

    print(f"Loaded samples, seeing {len(samples_df)} with {
        len(samples_df["condition"].unique())} unique conditions.")

    count_table_df.index = count_table_df['gene_id']
    count_table_df = count_table_df.drop(
        columns=['gene_id', 'gene_name']).astype(int).T
    count_table_df.sort_index(inplace=True)

    max_read_count_seen = pd.DataFrame(count_table_df.max(axis=0),
                                       columns=['max_read_count'])
    genes_to_keep = max_read_count_seen[max_read_count_seen['max_read_count'] >= minimal_gene_count].index
    filtered_table = count_table_df[genes_to_keep]
    print(f"Filtered table, keeping {len(genes_to_keep)} genes as at least one replicate has read count higher than {minimal_gene_count}.")
    print("--")

    return filtered_table, samples_df

def get_deseq_object(count_table_df, samples_df):
    dds = DeseqDataSet(
        counts=count_table_df,
        metadata=samples_df,
        design="condition",
        refit_cooks=True,
        inference=inference,
    )
    dds.fit_size_factors()
    dds.deseq2()
    return dds

def normalize_table(dds, count_table_df): # dds is a DeseqDataSet
    normed_counts = dds.layers['normed_counts']

    # scale to z-score on gene id
    normed_zscores = pd.DataFrame(normed_counts, 
                              columns=count_table_df.columns, 
                              index=count_table_df.index)
    for col in normed_zscores.columns:
        normed_zscores[col] = (normed_zscores[col] - normed_zscores[col].mean()) / normed_zscores[col].std()

    normed_counts = pd.DataFrame(normed_counts, 
                                 columns=count_table_df.columns, 
                                 index=count_table_df.index)

    return normed_counts, normed_zscores

def pairwise_dge_analysis(dds,
                      contrast_sample1,
                      contrast_sample2):
    ds = DeseqStats(dds, contrast=["condition",
                                   contrast_sample1,
                                   contrast_sample2],
                    inference=inference)

    # Run the Wald tests, with default Cook's distance filtering and independent filtering
    ds.summary()
    return ds

