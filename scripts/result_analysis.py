# result_analysis: Script that creates the log2 fold changes between target_samples and reference_samples and
# facilitates the creation of gene- and guide-level significance plots
# Last modified 02.10.2024
# ------------------------------------
import logging as log
from pathlib import Path
import pandas as pd
import numpy as np
from analysis_tools import run_script, assign_type

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)


class ResultAnalysis:
    def __init__(self):
        pass

    def create_drugz_log2fc(self, drugz_input, target_samples, reference_samples, essential_genes, non_essential_genes,
                            x_axis,
                            threshold_fdr, top):
        """
        Calculate the drugz log2 fold-changes between the target samples and reference samples

        :param drugz_input: Input file used for running the DrugZ analysis
        :param target_samples: Names of target samples
        :param reference_samples: Names of reference/control samples
        :param essential_genes: Csv file that contains list of essential genes to mark them in dataset (type='p')
        :param non_essential_genes: Csv file that contains list of non-essential genes to mark them in dataset 
        (type='n')
        :param x_axis: X-axis value for the final significance plots (either 'normZ' or 'log2 fold-change')
        :param threshold_fdr: Maximum false discovery rate that hits can have to be considered significant
        :param top: Maximum amount of significant hits that are displayed in the plot

        :return:
        """

        log_.debug("Computing log2 fold changes from drugz results")

        # Read the provided input files
        drugz_input = pd.read_csv(drugz_input, sep="\t")
        essential_genes = pd.read_csv(essential_genes)
        non_essential_genes = pd.read_csv(non_essential_genes)
        target_samples = target_samples.split(',')
        reference_samples = reference_samples.split(',')

        # Select the time points (conditions without the replicate number, e.g. t1_2d,...)
        timepoints = list(drugz_input.select_dtypes(include='number').columns.str.rsplit("_", n=1).str[0].unique())

        # Calculate the geometric means
        geometric_means = self.calculate_geometric_means(drugz_input, timepoints)

        log2df_all = None
        # Loop over all elements in target_samples
        for i in range(len(target_samples)):

            # Select the corresponding DrugZ output file
            target_sample = target_samples[i]
            reference_sample = reference_samples[i]
            comparison = f"{target_sample}-{reference_sample}"
            input_file_format = fr"drugz_{comparison}.txt"

            # Read the DrugZ output file
            input_file = pd.read_csv(fr".\drugz\{input_file_format}", sep="\t")

            # Rename the first column to ensure consistent naming
            input_file.rename(columns={input_file.columns[0]: "Gene"}, inplace=True)

            # Calculate the log2 fold-change between the current target and reference sample
            log2df = self.calculate_log2_fold_change(drugz_input, geometric_means, target_sample, reference_sample)

            # Place all log2 fold-changes per sgRNA in one row, separated by '|'
            log2_per_sgrna = self.get_log2_per_sgrna(log2df)

            # Take the mean of the log2 fold-changes for all replicates per gene
            log2_per_gene = self.get_log2_per_gene(log2df)

            # Merge the dataset with the log2 fold-changes per gRNA & gene
            log2df = pd.merge(input_file, log2_per_sgrna, on="Gene")
            log2df = pd.merge(log2df, log2_per_gene, on="Gene")

            # Add a column with the type of each gene to the dataframe
            log2df.insert(1, "type", assign_type(log2df["Gene"], essential_genes, non_essential_genes))

            # Sort the dataset by "rank_supp"
            log2df = log2df.sort_values(by="rank_supp")

            # Add the corresponding filename to the individual read_data columns
            new_column_names = [f"{comparison}.{col}" for col in log2df.columns[2:]]

            # Update the column names in the DataFrame
            log2df.columns = log2df.columns[:2].tolist() + new_column_names

            if log2df_all is None:
                log2df_all = log2df
            else:
                log2df_all = pd.merge(log2df_all, log2df, on=["Gene", "type"])

        # Store the dataframe in a csv file
        log2df_all.to_csv(r".\drugz_log2fcs_all.csv", index=False, sep=";")

        # Start the execution of R_analysis_2
        log_.info("Creating significance plots\n")
        log2fc_all = r".\drugz_log2fcs_all.csv"
        run_script(Path(__file__).parents[0] / "R_analysis_2.R",
                   additional_args=[log2fc_all, ",".join(target_samples), ",".join(reference_samples), x_axis,
                                    str(threshold_fdr),
                                    str(top)])

    @staticmethod
    def calculate_geometric_means(read_data, timepoints):
        """
        Calculate the geometric mean for each timepoint

        :param read_data: Cleaned data file with read counts per condition
        :param timepoints: List of timepoints
        :return: Geometric mean for each timepoint
        """

        geometric_means = {}
        for tp in timepoints:
            # Select all columns that belong to this time point
            tp_columns = [column for column in read_data.columns if tp == column.rsplit('_', 1)[0]]
            subset = read_data[tp_columns]

            # Calculate the geometric means
            log_values = np.log(subset.values)
            geometric_mean = np.exp(np.mean(log_values, axis=1))
            geometric_means[tp] = geometric_mean

        geometric_means = pd.DataFrame(geometric_means, index=read_data["sgRNA"])

        return geometric_means

    @staticmethod
    def calculate_log2_fold_change(read_data, geometric_means, target_sample, reference_sample):
        """
        Calculate the log2 fold change for the given target and reference sample

        :param read_data: Cleaned data file with read counts per condition
        :param geometric_means: Geometric means for each timepoint
        :param target_sample: Names of target sample
        :param reference_sample: Names of reference/control sample
        :return: Log2 fold change for the given target and reference sample
        """

        log2fc = np.log2(geometric_means[target_sample] / geometric_means[reference_sample])
        log2fc.reset_index(drop=True, inplace=True)
        log2df = pd.DataFrame({"Gene": read_data["Gene"], "l2fc": round(log2fc, 3)})

        return log2df

    @staticmethod
    def get_log2_per_sgrna(log2df):
        """
        Calculate the log2 fold change per gRNA

        :param log2df: Log2 fold change
        :return: Log2 fold change per gRNA
        """

        log2_per_sgrna = log2df.groupby("Gene")["l2fc"].apply(lambda x: "|".join(map(str, x)))

        return log2_per_sgrna

    @staticmethod
    def get_log2_per_gene(log2df):
        """
        Calculate the log2 fold change per gene

        :param log2df: Log2 fold change
        :return: Log2 fold change per gene
        """

        log2_per_gene = log2df.groupby("Gene")["l2fc"].mean().round(3).reset_index()
        log2_per_gene.columns = ["Gene", "Gene.log2fc"]

        return log2_per_gene
