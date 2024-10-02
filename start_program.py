# Python VERSION = "3.12"
# ---------------------------------------------------------------------------------------------------
# Script for the identification of drug-gene interactions in CRISPR screens using DrugZ
# Last modified 19 November 2023
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
# threshold: Minimum number of total reads/guide so that the guide will not be discarded
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
from pathlib import Path
import data_preparation
import result_analysis
from analysis_tools import run_script
import os
import logging as log

log.basicConfig(level=log.DEBUG)
log_ = log.getLogger(__name__)


def CRISPR_screen_analysis(args):
    # Set a working directory (the results will be inserted here)
    working_dir = args.working_dir
    os.chdir(working_dir)

    # Extract the dataset identifier
    dataset = os.path.basename(args.input_file).split('_')[0]
    log_.info(f"Starting analysis for {dataset}\n")

    # Define the names of the output files using the dataset identifier
    args.output_file = f"{dataset}_data_prep.csv"
    args.summary_file = f"{dataset}_counts_summary.csv"
    args.normalized_file = f"{dataset}_norm.csv"
    args.dataset = dataset

    # Create the output directory for the results files and plots
    args.output_folder = rf".\results"
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)
    os.chdir(args.output_folder)

    # Prepare the data for DrugZ
    log_.info(f"Preparing dataset for hit identification\n")
    data_preparation.data_preparation(args)

    # Perform the drugz analysis (additional optional parameters can be changed in 'run_drugz_try.py' if wanted)
    log_.info(f"Performing DrugZ analysis")
    run_script(rf"{Path(__file__).parents[0]}\run_drugz.py",
               additional_args=[args.target_samples, args.reference_samples])

    # Create drugZ log2 fold changes
    log_.info(f"Calculating log2 fold-changes between the target and reference sample\n")
    result_analysis.create_drugz_log2fc(drugz_input=f"{dataset}_drugz-input.txt", target_samples=args.target_samples,
                                        reference_samples=args.reference_samples,
                                        essential_genes=args.essential_genes,
                                        non_essential_genes=args.non_essential_genes, x_axis=args.x_axis,
                                        threshold_fdr=args.threshold_fdr,
                                        top=args.top)

    log_.info(f"Analysis for {dataset} complete\n")
