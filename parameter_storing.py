import pandas as pd
import datetime
import os
from pathlib import Path


class DGE_parameters:
    def __init__(self, args):
        self.count_matrix_filename = args.count_matrix
        self.count_matrix = pd.read_csv(args.count_matrix, sep='\t', comment='#')
        cm_columns = self.count_matrix.columns
        self.count_matrix.columns = ["gene_id", *cm_columns[1:]]

        self.sample_filename = args.sample_file
        self.sample_file = pd.read_csv(args.sample_file, sep='\t', header=None)
        # self.sample_file.columns = ["groupID", "sampleID"]
        self.sample_file.columns =  ['condition', 'sample']

        self.contrasts_filename = args.contrasts
        self.contrasts = pd.read_csv(args.contrasts, sep='\t', header=None)
        self.contrasts.columns = ["treatment", "control"]

        self.gene_annotation_resource = args.gtf_file

        # no of extreme genes to show annotation for
        self.annotate_extremes_no = 20
        self.padj_alpha=args.padj_alpha
        self.fold_change_threshold = args.fold_change_threshold
        self.minimal_gene_count = args.minimal_gene_count
        self.path_to_interesting_genes = args.interesting_genes_file
        self.organism_for_gprofiler = args.organism_for_gprofiler

        # if new parameter is added, it should be added to these lists too
        # they serve for reporting, no practical measure though
        self.__optional_parameters = [# self.organism_info,
                                      self.padj_alpha,
                                      self.fold_change_threshold,
                                      self.gene_annotation_resource,
                                      self.annotate_extremes_no,
            self.minimal_gene_count,
            self.path_to_interesting_genes,
        ]
        self.__optional_parameter_labels = [# "organism info",
                                            "Adjusted p-value maximum threshold",
                                            "fold change minimum threshold",
                                            "gene annotation resource",
                                            "no of extreme genes to annotate",
            "minimal gene count", "path to interesting genes"
        ]

        # output prep
        if args.output_directory_path is None:
            now = datetime.datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
            build = ["output-dge-analysis", now]
            for param, val in zip(["padj", "fc"],
                                  [self.padj_alpha, self.fold_change_threshold]):
                build.append(f"{param}={val}")
            args.output_directory_path = "_".join(build)

        self.output_dir = args.output_directory_path
        if Path(self.output_dir).exists():
            raise FileExistsError(f"The path {self.output_dir} already exists.")
        os.makedirs(self.output_dir, exist_ok=False)

        # recap params to the output
        with open(os.path.join(self.output_dir, "parameters.tsv"), mode='w') as writer:
            print("Parameter\tvalue\tinfo", file=writer)
            for filename, label in zip(
                    [self.count_matrix_filename, self.sample_filename, self.contrasts_filename],
                    ["count matrix", "sample file", "contrasts"]):
                print(f"{label}\t{filename}\trequired", file=writer)
            for value, label in zip(self.__optional_parameters, self.__optional_parameter_labels):
                stringified_values = str(value).replace("\n", " ")
                print(f"{label}\t{stringified_values}\t", file=writer)

    def report(self):
        # required -- filenames and other stuff
        print("Input files overview:")
        for filename, dataframe, label in zip(
                [self.count_matrix_filename, self.sample_filename, self.contrasts_filename],
                [self.count_matrix, self.sample_file, self.contrasts],
                ["count matrix", "sample file", "contrasts"]
        ):
            dfsize = len(dataframe)
            columns = ", ".join(dataframe.columns)
            print(f"{label}: taken from {filename}, it has {dfsize} data rows")
            print(f"detected columns: {columns}")
            print("--")

        print()
        print("Other parameters:")
        for value, label in zip(self.__optional_parameters, self.__optional_parameter_labels):
            stringified_values = str(value).replace("\n", " ")
            print(f"{label}: {stringified_values}")
        print("--")
        print()