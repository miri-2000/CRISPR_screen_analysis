# Python VERSION = "3.12"
# ---------------------------------------------------------------------------------------------------
# Script for the identification of drug-gene interactions in CRISPR screens using DrugZ
# Last modified 02.10.2024
# Free to modify and redistribute
# ---------------------------------------------------------------------------------------------------
# This file provides the full control over the entire analysis. Any inputs that need to be defined are stated here.
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
# the file requires a nohit row
# rows with non-targeting controls are expected to be called "Non-Targeting_Control"
# the sgRNA name and the sequence column from the table need to be unique
# the input file is expected to have guide_mm1_ columns: columns with 1bp mismatch on the
#       sgRNA (any possible 1bp mismatch, and if only possible from this sgRNA)
#       -> required for creating counts_summary file
# the input file has annotation columns that are identified by being the only non-numerical columns
# ---------------------------------------------------------------------------------------------------
import os
from pathlib import Path
import logging as log
from data_preparation import DataPreparation
import result_analysis
from analysis_tools import run_script

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)


class CRISPRScreenAnalysis:
    def __init__(self, args):
        """
        Initialize the CRISPRScreenAnalysis with user-specified arguments.
        Args:
            args: Argument object containing various parameters like input file, working directory, etc.
        """

        self.working_dir = args.working_dir
        self.input_file = args.input_file
        self.dataset = os.path.basename(self.input_file).split('_')[0]
        self.output_file = f"{self.dataset}_data_prep.csv"
        self.summary_file = f"{self.dataset}_counts_summary.csv"
        self.normalized_file = f"{self.dataset}_norm.csv"
        self.target_samples = args.target_samples
        self.reference_samples = args.reference_samples
        self.output_folder = Path(self.working_dir) / "results"
        self.essential_genes = args.essential_genes
        self.non_essential_genes = args.non_essential_genes
        self.library_file = args.library_file
        self.unwanted_columns = args.unwanted_columns
        self.unwanted_rows = args.unwanted_rows
        self.unwanted_row_substrings = args.unwanted_row_substrings
        self.threshold_reads = args.threshold_reads
        self.x_axis = args.x_axis
        self.threshold_fdr = args.threshold_fdr
        self.top = args.top
        self.distribution_condition1 = args.distribution_condition1
        self.distribution_condition2 = args.distribution_condition2
        self.replicate_type = args.replicate_type

        self.run_analysis()

    def run_analysis(self):
        """
        Executes the complete CRISPR screen analysis.
        """

        log_.info(f"Starting analysis for {self.dataset}\n")

        self.setup_output_directory()
        self.prepare_data()
        self.perform_drugz_analysis()
        self.calculate_log2fc()

        log_.info(f"Analysis for {self.dataset} complete")

    def setup_output_directory(self):
        """
        Sets up the output directory for storing results.
        """

        log_.info(f"Setting up the output directory\n")

        os.makedirs(self.output_folder, exist_ok=True)
        os.chdir(self.output_folder)
        log_.debug(f"Output directory set to: {self.output_folder}")

    def prepare_data(self):
        """
        Prepares data for DrugZ analysis.
        """

        log_.info(f"Preparing dataset for hit identification")

        data_preparation = DataPreparation()
        data_preparation.prepare_data(self.input_file, self.essential_genes, self.non_essential_genes,
                                           self.library_file, self.unwanted_columns, self.unwanted_rows,
                                           self.unwanted_row_substrings, self.summary_file, self.threshold_reads,
                                           self.output_file, self.dataset, self.replicate_type, self.target_samples,
                                           self.reference_samples, self.distribution_condition1,
                                           self.distribution_condition2)

    def perform_drugz_analysis(self):
        """
        Runs DrugZ analysis using the provided script.
        """

        log_.info(f"Performing DrugZ analysis")

        run_script(Path(__file__).parent / "run_drugz.py",
                   additional_args=[self.target_samples, self.reference_samples])

    def calculate_log2fc(self):
        """
        Calculates log2 fold-changes between the target and reference samples.
        """

        log_.info(f"Calculating log2 fold-changes between the target and reference sample")
        # Call result analysis function here
        result_analysis.create_drugz_log2fc(drugz_input=f"{self.dataset}_drugz-input.txt",
                                            target_samples=self.target_samples,
                                            reference_samples=self.reference_samples,
                                            essential_genes=self.essential_genes,
                                            non_essential_genes=self.non_essential_genes,
                                            x_axis=self.x_axis,
                                            threshold_fdr=self.threshold_fdr,
                                            top=self.top)
