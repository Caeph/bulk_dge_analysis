from gprofiler import GProfiler


def enrichment_analysis(query_ensembl_ids, all_tested_ensembl_ids, organism):
    gp = GProfiler(return_dataframe=True)

    enr = gp.profile(
        organism=organism,
        query=query_ensembl_ids,
        background=all_tested_ensembl_ids,
        sources=["GO:BP", "GO:MF", "GO:CC", "KEGG"],
    )
    return enr