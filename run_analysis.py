# Python VERSION = "3.10"
# ---------------------------------------------------------------------------------------------------
# run_analysis: Script for the identification of drug-gene interactions in CRISPR screens using DrugZ
# Last modified 19 November 2023
# Free to modify and redistribute
# ---------------------------------------------------------------------------------------------------
# This file provides the full control over the entire analysis. Any inputs that need to be defined are stated here.
# For details on how to adapt this file for new CRISPR screen analyses, please refer to the documentation.
# ------------------------------------
## Parameters
# input_file: txt file containing the raw read counts for each gRNA from the CRISPR screen
# essential_genes: csv file that contains list of essential genes to mark them in dataset (type='p')
# non_essential_genes: csv file that contains list of non-essential genes to mark them in dataset (type='n')
# library_file: txt file that is used as a reference to compare current dataset with
# output_file: output filename to write the output data file into
# summary_file: output filename to write the counts summary file into
# normalized_file: output filename to write the normalized read counts into
# type:# used as basis to create correlation plots (if type=='biological' analysis is done assuming replicates
#   are biological replicates; vice versa for 'technical')
# unwanted_columns: columns that should be removed from the input file
# unwanted_rows: rows that should be removed from the input file
# unwanted_row_substrings: substrings that should be removed from gRNA names
# threshold: minimal overall read count over all conditions for one gRNA
# target_samples: the names of the target conditions
# reference_samples: the names of the reference/control samples
# -----------------------------------------------------------------------------------------------------
## Expected structure of the input data file
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

import Data_preparation
import Result_analysis
from analysis_tools import run_script
import os
import logging as log

log.basicConfig(level=log.DEBUG)
log_ = log.getLogger(__name__)


# Define the input parameters
class Args:
    working_dir = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\results"
    essential_genes = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data\list_essentials_roderick.csv"
    non_essential_genes = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data\non_essentials_roderick.csv"
    library_file = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program data\broadgpp-brunello-library-contents.txt"
    input_file = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Datasets to be analysed\dataset\input_data\7234_all_Brunello_library_target_genes-req.txt"
    type = "biological"
    unwanted_columns= "guide_mm1_mismatch1,mismatch1_,nohit_cols,guide_mm1_nohit"
    unwanted_rows = ""
    unwanted_row_substrings = ":mismatch"
    threshold_reads = 0
    target_samples = "t2_2d_tr,t2_3d_tr,t1_3d,t1_3d2d,t2_3d2d_ut,t2_3d2d_tr,t2_2d_tr,t2_2d_ut,t2_3d_tr,t2_3d_ut"
    reference_samples = "t2_2d_ut,t2_3d_ut,t1_2d,t1_3d,t2_3d_ut,t2_3d_tr,t1_2d,t1_2d,t1_3d,t1_3d"
    x_axis = "normZ"
    threshold_fdr = 0.25
    top = 15
    distribution_condition1 = "t1_2d"
    distribution_condition2 = "t0"

def CRISPR_screen_analysis(args):
    # Set a working directory (the results will be inserted here)
    working_dir = args.working_dir
    # working_dir = r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Datasets to be analysed\dataset"
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
    Data_preparation.data_preparation(args)

    # Perform the drugz analysis (additional optional parameters can be changed in 'run_drugz_try.py' if wanted)
    log_.info(f"Performing DrugZ analysis")
    run_script(r"C:\Users\Miriam\PycharmProjects\CRISPR_screen_analysis\run_drugz.py",
               additional_args=[args.target_samples, args.reference_samples])

    # Create drugZ log2 fold changes
    log_.info(f"Calculating log2 fold-changes between the target and reference sample\n")
    Result_analysis.create_drugz_log2fc(drugz_input=f"{dataset}_drugz-input.txt", target_samples=args.target_samples,
                                        reference_samples=args.reference_samples,
                                        essential_genes=args.essential_genes,
                                        non_essential_genes=args.non_essential_genes, x_axis=args.x_axis,
                                        threshold_fdr=args.threshold_fdr,
                                        top=args.top)

    log_.info(f"Analysis for {dataset} complete\n")

# CRISPR_screen_analysis(Args)