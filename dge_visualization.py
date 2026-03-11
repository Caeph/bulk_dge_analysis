import os
os.environ["MPLBACKEND"] = "Agg"

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
plt.ioff()  # don’t open interactive windows

import seaborn as sns
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np

def plot_pca(filtered_table, log_counts, log_counts_annot,
             samples, ntop):
    # gene-wise variance, get ntop most variable genes
    genewise_std = np.std(filtered_table, axis=0
                          ).sort_values(ascending=False)
    most_variable_genes = genewise_std[:ntop].index

    pca = PCA(n_components=2)
    X = pca.fit_transform(log_counts[most_variable_genes])
    sample_map = {row["sample"]: row["condition"] for _, row in samples.iterrows()}
    fig, axes = plt.subplots(nrows=1,ncols=2, figsize=(10,7))
    sns.scatterplot(
        x=X[:, 0],
        y=X[:, 1],
        hue=[sample_map[x] for x in log_counts_annot],
        palette="Dark2",
        #  s=10
        ax=axes[0]
    )
    p1 = sns.scatterplot(
        x=X[:, 0],
        y=X[:, 1],
        hue=[sample_map[x] for x in log_counts_annot],
        palette="Dark2",
        #  s=10
        ax=axes[1]
    )

    left, right = axes[1].get_xlim()
    move = (right - left) * 0.01
    for pc1, pc2, label in zip(X[:,0],
                               X[:,1],
                               log_counts_annot):
        ha = "left"
        rot = 45
        p1.text(pc1 + move,
                pc2,
                label,
                horizontalalignment=ha,
                size='small',
                color='black',
                weight='semibold',
                rotation=rot)

    pc1_explained_ratio, pc2_explained_ratio = pca.explained_variance_ratio_
    for ax in axes:
        ax.set_xlabel(
            f"Principal component 1\nexplains {np.round(100*pc1_explained_ratio, decimals=2)}% variance"
        )
        ax.set_ylabel(
            f"Principal component 2\nexplains {np.round(100 * pc2_explained_ratio, decimals=2)}% variance"
        )

    axes[0].set_title("PCA without replicate annotation")
    axes[1].set_title("PCA with replicate annotation")

    plt.tight_layout()
    return fig

def plot_clustering(log_counts, log_counts_annot):
    # clustering table
    cor = np.corrcoef(log_counts)
    cor = pd.DataFrame(cor,
                       columns=log_counts_annot,
                       index=log_counts_annot)
    fig = sns.clustermap(cor.sort_index(),
                         cmap='Spectral')
    return fig.fig

def plot_MA(res, alpha, annotate_extremes_no):
    res = res.sort_values("padj", ascending=True)

    logmeancount = np.log2(res['baseMean'])
    res["log2Counts"] = logmeancount
    significant = (res["padj"] < alpha).map({True: "significant", False: "not significant"})

    fig, ax = plt.subplots(1,1,figsize=(10,7))

    p1 = sns.scatterplot(
        x=logmeancount,
        y=res['log2FoldChange'],
        hue=significant,
        s=5,
        palette={"significant": 'red', 'not significant': 'gray'},
        ax=ax
    )
    ax.set_xlabel("log2 mean count")
    ax.set_ylabel("log2 fold change")
    ax.legend().set_title(f'Decision for alpha={alpha}')

    plt.suptitle("MA plot")

    to_annotate = res.head(annotate_extremes_no)
    x_col, y_col = 'log2Counts', 'log2FoldChange'
    for line in range(0, to_annotate.shape[0]):
        p1.text(to_annotate[x_col].iloc[line], to_annotate[y_col].iloc[line],
                to_annotate["gene_name"].iloc[line], horizontalalignment='left',rotation=45,
                size='small', color='black')
    plt.tight_layout()
    return fig


def plot_volcano(res, alpha, fc_thr, annotate_extremes_no):
    #res = res.sort_values("padj", ascending=True)
    res = res.copy().dropna()

    res["gene_color"] = "under alpha, under FC thr."
    res.loc[(res['padj'] < alpha), "gene_color"] = "significant, under FC thr."
    # results_df['log2FoldChange']
    negative_mask = (res['log2FoldChange'] < -fc_thr) & (
            res['padj'] < alpha)
    res.loc[negative_mask, "gene_color"] = "significant, negative FC"
    positive_mask = (res['log2FoldChange'] > fc_thr) & (
            res['padj'] < alpha)
    res.loc[positive_mask, "gene_color"] = "significant, positive FC"
    palette = {
        "under alpha, under FC thr.": 'grey',
        "significant, under FC thr.": 'black',
        "significant, positive FC": "green",
        "significant, negative FC": "red"
    }

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    res["-logpadj"] = np.fmin(-np.log2(res['padj']), 1000)

    p1 = sns.scatterplot(
        x=res['log2FoldChange'],
        y=res["-logpadj"],
        hue=res["gene_color"],
        palette=palette,
        ax=ax,
        s=5
    )
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("- log2 adjusted p-value")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[1:], labels=labels[1:])


    to_annotate = res.head(annotate_extremes_no)
    x_col, y_col = 'log2FoldChange', '-logpadj'
    for line in range(0, to_annotate.shape[0]):
        ax.text(to_annotate[x_col].iloc[line], to_annotate[y_col].iloc[line],
                to_annotate["gene_name"].iloc[line], horizontalalignment='left',rotation=45,
                size='small', color='black')
    plt.tight_layout()
    return fig


_go_analysis_labeling = {
    "GO:BP": ["GO biological process", "blue"],
    "GO:MF": ["GO molecular function", 'green'],
    "GO:CC": ["GO cellular component", 'orange'],
    "KEGG": ["KEGG", "violet"]
}

def plot_go_enrichment(enriched_processes_df, label, ntop):
    # counts = enriched_processes_df[['source', 'name']].groupby('source').count().to_dict()['name']
    if len(enriched_processes_df) == 0:
        return None
    nrows = len(enriched_processes_df['source'].unique())
    if nrows > 1:
        fig, axes = plt.subplots(nrows=nrows,
                                 ncols=1,
                                 sharex=True,
                                 figsize=(15, 3.5*nrows),)
    else:
        fig, ax = plt.subplots(nrows=nrows,
                                 ncols=1,
                                 sharex=True,
                                 figsize=(15, 3.5*nrows), )
        axes = [ax]
    i = 0
    for name, group in enriched_processes_df.sort_values('name').groupby("source"):
        ax = axes[i]
        i += 1

        title, color = _go_analysis_labeling[name]

        to_plot = group.sort_values("precision", ascending=False)
        if len(to_plot) > ntop:
            to_plot = to_plot.head(ntop)
        sns.barplot(
            data=to_plot,
            y='name',
            x='precision',
           #  hue='source',
            ax=ax,
            color=color
        )
        ax.set_ylabel("")
        ax.set_title(title)
        ax.set_xlabel("gene ratio")

    fig.suptitle(
        f"Enrichment analysis: {label}",
        fontsize=20,  # larger font
        fontweight="bold",  # bold
        x=0.01,  # move to left (0 = far left, 0.5 = center)
        ha="left"  # left align text anchor
    )

    plt.tight_layout()
    return fig

def plot_informed_pairwise_heatmap(res,
                                   normed_counts, # only containing the needed samples
                                   normed_zscores, # dtto
                                   alpha,
                                   fc_thr,
                                   these_samples,
                                   max_to_show=100,
                                   average_read_counts=True
                                   ):
    # TODO check and make sure the things align row-wise

    res = res.copy()
    res['abs_log2foldchange'] = res["log2FoldChange"].abs()
    comparison_data = res[(res["padj"] < alpha) & (
                    res['abs_log2foldchange'] >= fc_thr)
        ].sort_values("abs_log2foldchange", ascending=False)
    comparison_data.index = comparison_data["gene_id"]

    if len(comparison_data) > max_to_show:
        print(f"{len(comparison_data)} genes changed significantly: clipping for visualization to top {max_to_show} by absolute log fold change")
        comparison_data = comparison_data.head(max_to_show)

    gene_ids = comparison_data['gene_id']

    figsize_y = np.fmax(len(comparison_data) * 0.2 + 5, 5)
    fig, axes = plt.subplots(nrows=1, ncols=3,
                             sharey=True,
                             figsize=(20, figsize_y),
                             width_ratios=[10, 1, 1])
    if len(comparison_data) == 0:
        return fig

    axes[0].set_title("Read counts (annotation) and z-score (color)\nz-score from entire experiment")
    current_tab = normed_counts.loc[:, gene_ids].T.copy()
    renaming_dict = pd.DataFrame(these_samples['condition'] + ':' + these_samples["sample"]).to_dict()[0]

    current_tab.columns = [renaming_dict[x] for x in current_tab.columns]
    current_tab = current_tab[sorted(current_tab.columns)]

    current_tab = pd.merge(current_tab,
                    comparison_data[['gene_name']], left_index=True, right_index=True).drop_duplicates()
    duplicated_gene_names = current_tab[current_tab["gene_name"].duplicated()]["gene_name"].unique()
    current_tab.loc[current_tab["gene_name"].isin(duplicated_gene_names), "gene_name"] = (current_tab.loc[
        current_tab["gene_name"].isin(duplicated_gene_names), "gene_name"] + ' (' +
        current_tab[current_tab["gene_name"].isin(duplicated_gene_names)].index + ")")
    # this must be done for comparison_data as well
    comparison_data.loc[comparison_data["gene_name"].isin(duplicated_gene_names), "gene_name"] = (comparison_data.loc[
                                                                                              comparison_data[
                                                                                                  "gene_name"].isin(
                                                                                                  duplicated_gene_names), "gene_name"] + ' (' +
                                                                                          comparison_data[comparison_data[
                                                                                              "gene_name"].isin(
                                                                                              duplicated_gene_names)].index + ")")

    current_tab.index = current_tab["gene_name"]
    current_tab = current_tab.drop(columns=["gene_name"])

    current_zscore_tab = current_tab.copy().T
    for col in current_zscore_tab.columns:
        current_zscore_tab[col] = (current_zscore_tab[col] - current_zscore_tab[col].mean()
                                   ) / current_zscore_tab[col].std()


    if average_read_counts:
        current_zscore_tab = current_zscore_tab.copy()
        gene_cols = current_zscore_tab.columns
        current_zscore_tab["samples"] = current_zscore_tab.index.str.split(":").str[0]
        current_zscore_tab = current_zscore_tab.groupby("samples").mean()

        genes_to_show = sorted(list(set(current_zscore_tab.columns)))
        current_zscore_tab = current_zscore_tab[genes_to_show].T

        current_tab = current_tab.T[genes_to_show]
        current_tab["samples"] = current_tab.index.str.split(":").str[0]

        def create_annot(x):
            mean = np.mean(x)
            std = np.std(x)
            return f"{mean:.2f} (s.d. {std:.2f})"

        current_tab = current_tab.groupby("samples").agg(
            {col: create_annot for col in gene_cols}
        )
        current_tab = current_tab[sorted(current_tab.columns)].T

        sns.heatmap(
            current_zscore_tab.T[comparison_data['gene_name']].T,
            ax=axes[0],
            annot=current_tab.T[comparison_data['gene_name']].T, fmt="",
            vmin=-2, vmax=2, cmap='Spectral_r',
            cbar_kws=dict(use_gridspec=True, location="left")
        )

    else:
        current_zscore_tab = current_zscore_tab[comparison_data['gene_name']].T
        sns.heatmap(
            current_zscore_tab,
            ax=axes[0],
            annot=np.round(current_tab.T[comparison_data['gene_name']].T), fmt=".5g",
            vmin=-2, vmax=2, cmap='Spectral_r',
            cbar_kws=dict(use_gridspec=True, location="left")
        )

    axes[0].set_ylabel("")

    comparison_data.index = comparison_data['gene_name']

    axes[1].set_title("log2 fold change")
    sns.heatmap(
        comparison_data[['log2FoldChange']],
        cmap='coolwarm',
        annot=True,
        ax=axes[1],
        cbar=False,
        vmin=-res['log2FoldChange'].abs().max(),
        vmax=res['log2FoldChange'].abs().max() + 0.01
    )

    axes[2].set_title("adjusted pvalue")
    sns.heatmap(
        np.fmin(-np.log2(comparison_data[['padj']]), 1000
                ),
        cmap='terrain',
        annot=comparison_data[['padj']],
        ax=axes[2],
        cbar=False,
    )
    plt.tight_layout()

    return fig

