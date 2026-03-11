import argparse
from itertools import combinations

from parameter_storing import DGE_parameters
from dge_toolbox import *

from enrichment_analysis import enrichment_analysis

from reporting_html import HTMLReporter
import os

parser = argparse.ArgumentParser()

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("--count_matrix", type=str,
                           help="""Path to tab delimited file with the following structure:  first columns contains GeneIDs, 
                    and each other column represents read counts (NOT NORMALIZED) for each sample.
                    """,
                           required=True
                           )
requiredNamed.add_argument("--sample_file", type=str,
                           help="""
                    A path to a tab delimited file with groupID in first column and sampleID in the second.
                    """,
                           required=True
                           )
requiredNamed.add_argument("--contrasts", type=str,
                           help="""
                    A path to a tab delimited file indicating which groups should be compared in a pairwise manner. 
                    Keep always the same structure: 1st column - treatment, 2nd column - control
                    """,
                           required=True
                           )
requiredNamed.add_argument(
    "--gtf_file", default=None, type=str, help="""
    File for GTF resource, will be used only if gene_annotation_resource is set to GTF
    """
)

parser.add_argument("--padj_alpha", default=0.05, type=float,
                    help="""
                    FDR maximum threshold for results selection.
                    """)
parser.add_argument("--fold_change_threshold", default=1, type=float,
                    help="""
                    Minimal log fold change threshold
                    """)
parser.add_argument("--output_directory_path", default=None, type=str,
                    help="""
                    Path to store the directory. If the directory exists, the run will fail.
                    """
                    )
parser.add_argument("--minimal_gene_count", default=15, type=int,
                    help="""
                    Count table is filtered to contain the genes that
    have at least {minimal gene count} in at least one of the replicates.
                    """)

parser.add_argument("--interesting_genes_file", default=None, type=str,
                    help="""
                    Path to ENSEMBL IDs that we want to show in more detail. The file should contain one gene ID per line, no header.
                    Detailed information will be shown only for those that remain after minimal gene number filtering.
                    """)
parser.add_argument("--organism_for_gprofiler", default=None, type=str, help="""
                    Code of organism for gprofiler (hsapiens for human). If not set, enrichment analysis is not done.
                    """)

def pairwise_enrichment_analysis(res, parameters,
                                 reporter, all_tested, contrast_sample1,
                                 contrast_sample2, directory):
    if parameters.organism_for_gprofiler is not None:
        figure_dir = os.path.join(directory, "enrichment_figures")
        os.makedirs(directory, exist_ok=True)
        os.makedirs(figure_dir, exist_ok=True)

        alpha = parameters.padj_alpha
        fc_thr = parameters.fold_change_threshold
        print("Performing gene enrichment analysis...")

        upregulated_genes = list(res[(res['padj'] < alpha) & (res['log2FoldChange'] > fc_thr)]['gene_id'].unique())
        downregulated_genes = list(res[(res['padj'] < alpha) & (res['log2FoldChange'] < -fc_thr)]['gene_id'].unique())

        upregulated_enrichment = enrichment_analysis(upregulated_genes,
                                                     all_tested,
                                                     parameters.organism_for_gprofiler)
        downregulated_enrichment = enrichment_analysis(downregulated_genes,
                                                       all_tested,
                                                       parameters.organism_for_gprofiler)

        # save enrichment to file
        upregulated_enrichment.to_csv(os.path.join(directory,
                                                   "upregulated_genes_enrichment_analysis.tsv"), sep="\t")
        downregulated_enrichment.to_csv(os.path.join(directory,
                                                   "downregulated_genes_enrichment_analysis.tsv"), sep="\t")

        # plot
        reporter.start_section("Enrichment analysis")
        ntop = parameters.annotate_extremes_no
        reporter.send_text(f"""
        Gene ratio is calculated as the ratio between genes corresponding with the GO term AND
        all upregulated/downregulated genes. The closer to one, the more important the result.
        Only top {ntop} results with p-value > 0.05 are shown. All significant are stored in tsv files.
        """)
        fig = plot_go_enrichment(upregulated_enrichment, "upregulated genes", ntop)
        reporter.send_text(f"As upregulated genes, we understand all those with adjusted p-val < {alpha},"
                           f"log2 fold change > {fc_thr}.")
        reporter.send_figure(fig, "Upregulated genes: gene ratio bar plot",
                             save_figure_path=os.path.join(figure_dir,
                                    f"upregulated_enrichment__{contrast_sample1}_vs_{contrast_sample2}.png"))
        plt.close(fig)
        fig = plot_go_enrichment(downregulated_enrichment, "downregulated genes", ntop)
        reporter.send_text(f"As downregulated genes, we understand all those with adjusted p-val < {alpha},"
                           f"log2 fold change < -{fc_thr}.")
        reporter.send_figure(fig, "Downregulated genes: gene ratio bar plot",
                             save_figure_path=os.path.join(figure_dir,
                                                           f"downregulated_enrichment__{contrast_sample1}_vs_{contrast_sample2}.png"))
        plt.close(fig)


def pairwise_analysis(reporter, dds,
                      contrast_sample1,
                      contrast_sample2,
                      gene_names,
                      count_table_df,
                      normalized_count_table,
                      normalized_zscore_table,
                      interesting_genes,
                      directory,
                      parameters
                      ):
    os.makedirs(directory, exist_ok=True)

    print("Performing pairwise comparison...")

    ds = pairwise_dge_analysis(dds, contrast_sample1, contrast_sample2)
    res = ds.results_df  # columns include: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    res["contrast_sample1"] = contrast_sample1
    res["contrast_sample2"] = contrast_sample2

    res = res.merge(gene_names, how="inner",
                    left_index=True, right_index=True).sort_values(by='padj', ascending=True)
    # save to csvs
    res.to_csv(os.path.join(directory, "pairwise_analysis_FULL.tsv"), sep="\t")
    alpha = parameters.padj_alpha
    res[res["padj"] < alpha].to_csv(os.path.join(directory,
                                                 "pairwise_analysis_SIGNIFICANT.tsv"), sep="\t")
    fc_thr = parameters.fold_change_threshold
    res[(res["padj"] < parameters.padj_alpha) & (
        res["log2FoldChange"].abs() > fc_thr)].to_csv(os.path.join(directory,
                f"pairwise_analysis_SIGNIFICANT_absFC>{fc_thr}.tsv"), sep="\t")

    figure_dir = os.path.join(directory, "figures")
    os.makedirs(figure_dir, exist_ok=True)

    # MA plot
    samples = parameters.sample_file
    this_samples = samples[samples["condition"].isin([contrast_sample1, contrast_sample2])
        ].sort_values("condition")
    subset_table = count_table_df.loc[this_samples.index]
    quality_control(subset_table,
                    this_samples,
                    reporter,
                    figure_dir
                    )

    reporter.start_section("MA plot")
    fig = plot_MA(res, alpha, parameters.annotate_extremes_no)
    reporter.send_figure(fig, caption="MA plot",
                         save_figure_path=os.path.join(figure_dir,
                                        f"MAplot__{contrast_sample1}_vs_{contrast_sample2}.png"))
    plt.close(fig)

    # volcano plot
    reporter.start_section("Volcano plot")
    fig = plot_volcano(res, alpha, fc_thr, parameters.annotate_extremes_no)
    reporter.send_figure(fig, caption="Volcano plot (-log2 p-value is capped at 1000 to show all examples)",
                         save_figure_path=os.path.join(figure_dir,
                                                       f"Volcanoplot__{contrast_sample1}_vs_{contrast_sample2}.png"))
    plt.close(fig)

    reporter.start_section("General results: most changed genes")
    # # heatmap with normed counts
    subset_normalized_count_table = normalized_count_table.loc[this_samples.index]
    subset_normalized_zscore_table = normalized_zscore_table.loc[this_samples.index]
    fig = plot_informed_pairwise_heatmap(res,
                                         subset_normalized_count_table,
                                         subset_normalized_zscore_table,
                                         alpha, fc_thr, this_samples, average_read_counts=False)
    reporter.send_figure(fig, "Heatmap with normalized read counts per sample, "
                              "colored by per-comparison z-score, sorted by absolute log2 fold change",
                         save_figure_path=os.path.join(figure_dir,
                                                       f"heatmap__{contrast_sample1}_vs_{contrast_sample2}.png"))
    plt.close(fig)

    fig = plot_informed_pairwise_heatmap(res,
                                         subset_normalized_count_table,
                                         subset_normalized_zscore_table,
                                         alpha, fc_thr, this_samples, average_read_counts=True)
    reporter.send_figure(fig, "Heatmap with normalized read counts per sample, "
                              "colored by per-comparison z-score, sorted by absolute log2 fold change",
                         save_figure_path=os.path.join(figure_dir,
                                                       f"heatmapAvg__{contrast_sample1}_vs_{contrast_sample2}.png"))
    plt.close(fig)

    # if interesting genes are defined, do that as well
    if interesting_genes is not None:

        interesting_res = pd.merge(res, interesting_genes,
                                   left_index=True,
                                   right_on='interesting_gene_id').drop(
            columns=["interesting_gene_id"])
        if len(interesting_res) > 0:
            reporter.start_section(f"General results: interesting genes ({parameters.path_to_interesting_genes})")
            fig = plot_informed_pairwise_heatmap(interesting_res,
                                             subset_normalized_count_table,
                                             subset_normalized_zscore_table,
                                             alpha,
                                        0, # we wanna show as many interesting genes as possible if the change is statistically significant
                                                 this_samples, average_read_counts=False)
            reporter.send_figure(fig, "Heatmap with normalized read counts per sample,\n"
                                      "colored by per-comparison z-score, sorted by absolute log2 fold change.\n"
                                      "Only significantly changed genes are shown, no fold change filter is applied apart from visualization clipping.",
                                 save_figure_path=os.path.join(figure_dir,
                                                               f"heatmapOnInteresting__{contrast_sample1}_vs_{contrast_sample2}.png"))
            plt.close(fig)

            fig = plot_informed_pairwise_heatmap(interesting_res,
                                                 subset_normalized_count_table,
                                                 subset_normalized_zscore_table,
                                                 alpha,
                                                 0, this_samples, average_read_counts=True)
            reporter.send_figure(fig, "Heatmap with normalized read counts per sample,\n"
                                      "colored by per-comparison z-score, sorted by absolute log2 fold change.\n"
                                      "Only significantly changed genes are shown, no fold change filter is applied apart from visualization clipping.",
                                 save_figure_path=os.path.join(figure_dir,
                                                               f"heatmapOnInterestingAvg__{contrast_sample1}_vs_{contrast_sample2}.png"))
            plt.close(fig)
    return res[(res["padj"] < parameters.padj_alpha) & (
        res["log2FoldChange"].abs() > fc_thr)]


def plot_bulk_heatmap(normed_counts, normed_zscores, selected_gene_ids, samples, gene_names_mapper, average_values=True):
    subset_counts = normed_counts[selected_gene_ids]

    renaming_dict = pd.DataFrame(samples['condition'] + ':' + samples["sample"]).to_dict()[0]

    current_tab = subset_counts.T.copy()
    current_tab.columns = [renaming_dict[x] for x in current_tab.columns]
    current_tab = current_tab[sorted(current_tab.columns)]

    current_tab.index = [gene_names_mapper[x] if (x in gene_names_mapper) else x for x in current_tab.index ]

    current_zscore_tab = current_tab.copy().T
    for col in current_zscore_tab.columns:
        current_zscore_tab[col] = (current_zscore_tab[col] - current_zscore_tab[col].mean()
                                   ) / current_zscore_tab[col].std()

    fig, ax = plt.subplots(1, 1, figsize=(1 * len(samples), 0.3 * len(selected_gene_ids)))
    if not average_values:
        sns.heatmap(
            current_zscore_tab.T,
            annot=np.round(current_tab),
            vmin=-3, vmax=3, fmt=".8g",
            ax=ax, cmap="Spectral_r"
        )
    else:
        current_zscore_tab['samples'] = current_zscore_tab.index.str.split(':').str[0]
        current_zscore_tab = current_zscore_tab.groupby("samples").mean()

        current_tab = current_tab.T
        cols = current_tab.columns
        current_tab['samples'] = current_tab.index.str.split(':').str[0]

        def create_annot(x):
            mean = np.mean(x)
            std = np.std(x)
            return f"{mean:.2f} (s.d. {std:.2f})"

        current_tab = current_tab.groupby("samples").agg({x: create_annot for x in cols})

        sns.heatmap(
            current_zscore_tab.T,
            annot=current_tab.T,
            vmin=-3, vmax=3, fmt="",
            ax=ax, cmap="Spectral_r"
        )
    return fig

def plot_individual_gene_boxplot(gene_id, gene_name, samples,
                                 normed_counts, collected_res, contrasts): # these should not be saved
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(20, 8))
    on_gene_counts = normed_counts[[gene_id]].merge(samples, left_index=True, right_index=True)
    on_gene_res = collected_res[collected_res.index == gene_id]

    ax = axes[0]

    sns.boxplot(
        data=on_gene_counts,
        x='condition',
        y=gene_id,
        ax=ax
    )
    sns.stripplot(
        data=on_gene_counts,
        x='condition',
        y=gene_id,
        ax=ax
    )
    plt.sca(ax)
    plt.xticks(rotation=90)
    ax.set_ylabel("normalized read counts")
    ax.set_xlabel("")
    ax.set_title(f"Normalized counts for gene {gene_name}")

    ax = axes[1]
    all_conditions = samples['condition'].unique()

    lfc_tab = np.zeros((len(all_conditions), len(all_conditions)))
    lfc_tab_annot = np.zeros((len(all_conditions), len(all_conditions))).astype(str)

    for i, cond1 in enumerate(all_conditions):
        for j, cond2 in enumerate(all_conditions):
            # exclude those that are not in contrasts
            seen_in_contrasts = contrasts[(contrasts['treatment'] == cond1) & (contrasts['control'] == cond2)]
            if len(seen_in_contrasts) == 0:
                lfc_tab[i, j] = np.nan
                lfc_tab_annot[i, j] = ""

            lfc_val = on_gene_res[(on_gene_res['contrast_sample1'] == cond1
                                   ) & (on_gene_res['contrast_sample2'] == cond2)]['log2FoldChange']
            if len(lfc_val) > 0:  # a value is available
                lfc_tab[i, j] = lfc_val.iloc[0]
                lfc_tab_annot[i, j] = np.round(lfc_val.iloc[0], decimals=1)
            else:  # unsignificant -- keep color as 0
                lfc_tab_annot[i, j] = 'ns./small'

    lfc_tab = pd.DataFrame(lfc_tab, columns=all_conditions, index=all_conditions)

    max_observed = lfc_tab.abs().max().max()
    sns.heatmap(
        lfc_tab,
        annot=lfc_tab_annot, fmt='',
        cmap='PuOr',
        vmin=-max_observed-0.01, vmax=max_observed+0.01,
        ax=ax
    )
    ax.set_title(f"Log2 fold change in pairwise contrasts")
    ax.set_xlabel("treatment")
    ax.set_ylabel("control")

    plt.suptitle(gene_name)
    plt.tight_layout()
    return fig

def plot_lfc_heatmap(collected_res, ntop):
    collected_res["comparison"] = collected_res["contrast_sample1"] + '_vs_' + collected_res["contrast_sample2"]
    counts_per_comparison = collected_res.reset_index(drop=True)[['gene_id', 'comparison', 'log2FoldChange']
    ].groupby('gene_id').agg(
        {"comparison": "count", "log2FoldChange": lambda x: x.abs().max()}
    )

    if len(counts_per_comparison) > ntop:
        top_frequent = counts_per_comparison.sort_values(
        ["comparison", 'log2FoldChange'], ascending=False).head(ntop)
    else:
        top_frequent = counts_per_comparison.sort_values(["comparison", 'log2FoldChange'], ascending=False)

    genes_to_show = top_frequent.index.values

    subset = collected_res[collected_res.index.isin(genes_to_show)]

    tab = pd.pivot_table(
        subset,
        values='log2FoldChange',
        columns='comparison',
        index='gene_name', aggfunc='max')
    max_observed = tab.abs().max().max()

    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20, 0.5 * len(tab.T)))
    sns.heatmap(tab.T,
                cmap='PuOr',
                vmin=-max_observed-0.01, vmax=max_observed+0.01,
                ax=axes, annot=True
                )

    return fig, genes_to_show


def bulk_overview_results(collected_res, count_table_df, normed_counts, normed_zscores,
                          interesting_genes, gene_names, parameters, reporter, ntop=20):
    samples = parameters.sample_file
    genewise_std = np.std(normed_counts, axis=0
                          ).sort_values(ascending=False)
    most_variable_genes = genewise_std[:ntop].index
    gene_names_mapper = gene_names[['gene_name']].to_dict()['gene_name']

    reporter.start_section("Illustratory read count heatmaps")
    fig = plot_bulk_heatmap(normed_counts, normed_zscores, most_variable_genes, samples, gene_names_mapper,
                            average_values=False)
    reporter.send_figure(fig,
                         f"Top {ntop} most variable genes: read count heatmap",
                         save_figure_path=os.path.join(parameters.output_dir,
                                                       'central_analysis_figures',
                                                       f"top{ntop}_read_count_heatmap.png"))
    plt.close(fig)

    fig = plot_bulk_heatmap(normed_counts, normed_zscores, most_variable_genes, samples, gene_names_mapper,
                            average_values=True)
    reporter.send_figure(fig,
                         f"Top {ntop} most variable genes: read count heatmap",
                         save_figure_path=os.path.join(parameters.output_dir,
                                                       'central_analysis_figures',
                                                       f"top{ntop}_read_count_heatmap_AVG.png")
                         )
    plt.close(fig)

    if interesting_genes is not None:
        genewise_std = np.std(normed_counts[interesting_genes['interesting_gene_id'].values], axis=0
                              ).sort_values(ascending=False)
        most_variable_interesting_genes = genewise_std[:ntop].index

        reporter.start_section(
            f"Illustratory read count heatmaps: interesting genes ({parameters.path_to_interesting_genes})")
        fig = plot_bulk_heatmap(normed_counts, normed_zscores, most_variable_interesting_genes, samples,
                                gene_names_mapper, average_values=False)
        reporter.send_figure(fig,
                             f"Top {ntop} most variable interesting genes: read count heatmap",
                             save_figure_path=os.path.join(parameters.output_dir,
                                                           'central_analysis_figures',
                                                           f"top{ntop}_interesting_read_count_heatmap.png")
                             )
        plt.close(fig)

        fig = plot_bulk_heatmap(normed_counts, normed_zscores, most_variable_interesting_genes, samples,
                                gene_names_mapper, average_values=True)
        reporter.send_figure(fig,
                             f"Top {ntop} most variable interesting genes: read count heatmap",
                             save_figure_path=os.path.join(parameters.output_dir,
                                                           'central_analysis_figures',
                                                           f"top{ntop}_interesting_read_count_heatmap_AVG.png")
                             )
        plt.close(fig)


    # plot overview of most frequently changed genes
    reporter.start_section("Illustration of most frequently changed genes")
    fig, genes_to_show = plot_lfc_heatmap(collected_res, ntop)
    reporter.send_figure(fig,"Heatmap of most frequently changed genes vs. comparison, "
                             "each position is (statistically significant) log2 fold change",
                         save_figure_path=os.path.join(parameters.output_dir,
                                                       'central_analysis_figures',
                                                       f"top{ntop}_freq_changed_heatmap.png")
                         )
    plt.close(fig)

    if interesting_genes is not None:
        reporter.start_section(f"Illustration of most frequently changed genes: interesting genes ({parameters.path_to_interesting_genes})")
        fig, interesting_genes_to_show = plot_lfc_heatmap(collected_res[collected_res.index.isin(interesting_genes['interesting_gene_id'].values)],
                               ntop)
        reporter.send_figure(fig, "Heatmap of most frequently changed genes vs. comparison, "
                                  "each position is (statistically significant) log2 fold change",
                             save_figure_path=os.path.join(parameters.output_dir,
                                                           'central_analysis_figures',
                                                           f"top{ntop}_interesting_freq_changed_heatmap.png")
                             )
        plt.close(fig)
    else:
        interesting_genes_to_show = []


    reporter.start_section("Visualization for individual genes: most frequently changed genes (table above)")
    reporter.send_text("For individual genes, a pair of plots is shown:"
                       "on the left, normed read count per condition"
                       "on the right, heatmap of log2 fold changes (ns. if no significant change was seen)"
                       "")
    for gene_id in genes_to_show:
        gene_name = gene_id
        if gene_id in gene_names_mapper:
            gene_name = gene_names_mapper[gene_id]
        fig = plot_individual_gene_boxplot(gene_id, gene_name, samples,
                                 normed_counts, collected_res,
                                           parameters.contrasts)

        reporter.send_figure(fig,
                             f"Gene {gene_name} ({gene_id})")
        plt.close(fig)

    if interesting_genes is not None:
        reporter.start_section(
            f"Visualization for individual genes: most frequently changed interesting genes ({parameters.path_to_interesting_genes})")
        reporter.send_text("For individual genes, a pair of plots is shown:"
                           "on the left, normed read count per condition"
                           "on the right, heatmap of log2 fold changes (ns. if no significant change was seen)"
                           "")
        for gene_id in interesting_genes_to_show:
            gene_name = gene_id
            if gene_id in gene_names_mapper:
                gene_name = gene_names_mapper[gene_id]
            fig = plot_individual_gene_boxplot(gene_id, gene_name, samples,
                                               normed_counts, collected_res,
                                               parameters.contrasts)

            reporter.send_figure(fig,
                                 f"Gene {gene_name} ({gene_id}), from the \"interesting\" set")
            plt.close(fig)







def main(args):
    parameters = DGE_parameters(args)
    parameters.report()

    central_reporter = HTMLReporter(f"Bulk RNA-seq report on {args.count_matrix}")

    # this probably needs more variability but is fine for human
    gene_names = load_gene_annotation(
        parameters.gene_annotation_resource)

    count_table_df, samples_df = combine_to_input(
        parameters.count_matrix,
        parameters.sample_file,
        parameters.minimal_gene_count
    )

    interesting_genes = None
    if parameters.path_to_interesting_genes is not None:
        interesting_genes = pd.read_csv(parameters.path_to_interesting_genes, header=None)
        interesting_genes.columns = ['interesting_gene_id']

        interesting_genes = interesting_genes[interesting_genes['interesting_gene_id'].isin(
            count_table_df.columns
        )]

    quality_control(count_table_df,
                    samples_df,
                    central_reporter,
                    os.path.join(parameters.output_dir,
                                 "central_analysis_figures")
                    )

    dds = get_deseq_object(count_table_df, samples_df)
    normalized_count_table, normalized_zscore_table = normalize_table(dds, count_table_df)

    significant_results = []

    for i, row in parameters.contrasts.iterrows():
        contrast_sample1, contrast_sample2 = row['treatment'], row['control']
        print(f"Now processing {contrast_sample1} and {contrast_sample2} pairwise analysis")
        reporter = HTMLReporter(f"Pairwise analysis on {contrast_sample1} and {contrast_sample2}")
        directory = os.path.join(parameters.output_dir,
                                       f"pairwise_analysis__{contrast_sample1}__vs__{contrast_sample2}")
        res = pairwise_analysis(reporter,
                             dds,
                          contrast_sample1,
                          contrast_sample2,
                          gene_names,
                          count_table_df,
                          normalized_count_table,
                          normalized_zscore_table,
                          interesting_genes,
                          directory,
                          parameters
                          )
        pairwise_enrichment_analysis(res,
                                     parameters,
                                     reporter,
                                     list(count_table_df.columns),
                                     contrast_sample1, contrast_sample2,
                                     directory)



        reporter.save(os.path.join(parameters.output_dir,
                                   f"report__{contrast_sample1}_vs_{contrast_sample2}.html"))
        significant_results.append(res)


    # save reporters
    significant_results = pd.concat(significant_results)
    print("Creating bulk visualizations...")
    bulk_overview_results(significant_results, count_table_df,
                          normalized_count_table,
                          normalized_zscore_table,
                          interesting_genes, gene_names, parameters,
                          central_reporter
                          )
    significant_results.to_csv(os.path.join(parameters.output_dir,
                                            "all_pairwise_comparison_significant_highabsFC.tsv.gz"), sep='\t')
    central_reporter.save(os.path.join(parameters.output_dir,
                                       "central_report.html"))

if __name__ == '__main__':
    args = parser.parse_args()
    main(args)