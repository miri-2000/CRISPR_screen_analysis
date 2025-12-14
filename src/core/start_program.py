# Python VERSION = "3.12"
# ---------------------------------------------------------------------------------------------------
# Script for the identification of drug-gene interactions in CRISPR screens using DrugZ
# Last modified 02.10.2024
# Free to modify and redistribute
# ---------------------------------------------------------------------------------------------------
# This file provides the full control over the entire core. Any inputs that need to be defined are stated here.
# For details on how to adapt this file for new CRISPR screen analyses, please refer to the documentation.
# ---------------------------------------------------------------------------------------------------
# Parameters

# target_samples: Treated conditions (e.g. with exposure to drugs,...) to be compared against baseline conditions
# reference_samples: Baseline conditions (e.g. with no added drugs)
# input_file: Text file containing the raw read counts for each gRNA from the CRISPR screen
# essential_genes: CSV file that contains list of essential genes to mark them in dataset (type='p')
# non_essential_genes: CSV file that contains list of non-essential genes to mark them in dataset (type='n')
# library_file: Text file that is used as a reference to compare current dataset with
# unwanted_columns: Columns that should be removed from the input file
# unwanted_rows: Rows that should be removed from the input file
# unwanted_row_substrings: Row substrings that should be removed from gRNA names
# threshold_reads: Minimum number of total reads/guide so that the guide will not be discarded
# x_axis: Metric that should be taken for the x axis of the distribution plot
# threshold_fdr: Maximum FDR for a gene to be considered significant
# top: Maximum number of significant genes that will be displayed
# distribution_condition1: Sample that should be taken as the positive control of the screen
# distribution_condition2: Sample that should be taken as the negative control of the screen
# working_dir: Directory where the results should be stored
# replicate_type:Defines whether replicate samples are biological or technical replicates
# -----------------------------------------------------------------------------------------------------
# Expected structure of the input data file
# all annotation columns contain non-numerical data while data columns only have numerical values
# the first two columns contain the sgRNA names and sequences
# the sgRNA name and the sequence column from the table need to be unique
# the file requires a nohit row
# rows with non-targeting controls (if any) are expected to be called "Non Targeting Control" (independent on casing
# and separation symbol (possible symbols: "_", "-", ".")
# the input file is expected to have guide_mm1_ columns: columns with 1bp mismatch on the
#       sgRNA (any possible 1bp mismatch, and if only possible from this sgRNA)
#       -> required for creating counts_summary file
# the input file has annotation columns that are identified by being the only non-numerical columns
# ---------------------------------------------------------------------------------------------------
import os
from pathlib import Path
import logging as log
from src.core.data_preparation import DataPreparation
from src.core.result_analysis import ResultAnalysis
from src.core.analysis_tools import run_script

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)


class CRISPRScreenAnalysis:
    def __init__(self):
        pass

    def run_analysis(self, working_dir, input_file, target_samples, reference_samples, essential_genes,
                     non_essential_genes, library_file, unwanted_columns, unwanted_rows, unwanted_row_substrings,
                     threshold_reads, x_axis, threshold_fdr, top, distribution_condition1, distribution_condition2,
                     replicate_type):
        """
        Executes the complete CRISPR screen core.
        """

        dataset = os.path.basename(input_file).split('_')[0]
        log_.info(f"Starting core for {dataset}\n")

        self.setup_output_directory(working_dir)

        self.prepare_data(input_file, essential_genes, non_essential_genes, library_file, unwanted_columns,
                          unwanted_rows, unwanted_row_substrings, threshold_reads, dataset, replicate_type,
                          target_samples, reference_samples, distribution_condition1, distribution_condition2)

        self.perform_drugz_analysis(working_dir, target_samples, reference_samples)

        self.calculate_log2fc(dataset, target_samples, reference_samples, essential_genes, non_essential_genes, x_axis,
                              threshold_fdr, top)

        log_.info(f"Analysis for {dataset} complete")

    @staticmethod
    def results_folder(working_dir):
        """Return the path to the results folder for a given working directory."""
        return Path(working_dir) / "results"

    @staticmethod
    def setup_output_directory(working_dir):
        """
        Sets up the output directory for storing results.
        """

        log_.info(f"Setting up the output directory\n")

        output_folder = CRISPRScreenAnalysis.results_folder(working_dir)
        os.makedirs(output_folder, exist_ok=True)
        os.chdir(output_folder)
        log_.debug(f"Output directory set to: {output_folder}")

    @staticmethod
    def prepare_data(input_file, essential_genes, non_essential_genes, library_file, unwanted_columns, unwanted_rows,
                     unwanted_row_substrings, threshold_reads, dataset, replicate_type, target_samples,
                     reference_samples, distribution_condition1, distribution_condition2):
        """
        Prepares data for DrugZ core.
        """

        log_.info(f"Preparing dataset for hit identification")

        output_file = f"{dataset}_data_prep.csv"
        summary_file = f"{dataset}_counts_summary.csv"

        data_preparation = DataPreparation()
        data_preparation.prepare_data(input_file, essential_genes, non_essential_genes, library_file, unwanted_columns,
                                      unwanted_rows, unwanted_row_substrings, summary_file, threshold_reads,
                                      output_file, dataset, replicate_type, target_samples, reference_samples,
                                      distribution_condition1, distribution_condition2)

    @staticmethod
    def perform_drugz_analysis(working_dir, target_samples, reference_samples):
        """
        Runs DrugZ core using the provided script.
        """

        log_.info(f"Performing DrugZ core")

        input_folder = CRISPRScreenAnalysis.results_folder(working_dir)

        run_script(Path(__file__).parent / "run_drugz.py",
                   additional_args=[input_folder,target_samples, reference_samples])

    @staticmethod
    def calculate_log2fc(dataset, target_samples, reference_samples, essential_genes, non_essential_genes, x_axis,
                         threshold_fdr, top):
        """
        Calculates log2 fold-changes between the target and reference samples.
        """

        log_.info(f"Calculating log2 fold-changes between the target and reference sample")

        # Call result core function here
        result_analysis = ResultAnalysis()
        result_analysis.create_drugz_log2fc(drugz_input=f"{dataset}_drugz-input.txt",
                                            target_samples=target_samples,
                                            reference_samples=reference_samples,
                                            essential_genes=essential_genes,
                                            non_essential_genes=non_essential_genes,
                                            x_axis=x_axis,
                                            threshold_fdr=threshold_fdr,
                                            top=top)
