# run_drugz: This script is used to execute DrugZ and allows the
# modification of additional parameters for running DrugZ if wanted
# by the user.
# Last modified 02.10.2024
# ------------------------------------
import drugZ
import os
import fnmatch
import sys


class Args:
    infile = os.path.abspath(fnmatch.filter(os.listdir('..'), '*drugz-input.txt')[0])
    drugz_output_file = "./drugz/drugz_"
    gRNA_outfile = "./drugz/gRNA_"
    target_samples = sys.argv[1]
    reference_samples = sys.argv[2]
    remove_genes = None
    unpaired = True
    replicates = True
    pseudocount = 5
    half_window_size = 500


target_samples_list = Args.target_samples.split(',')
reference_samples_list = Args.reference_samples.split(',')

# Create a DrugZ output folder (if it does not exist yet)
os.makedirs("./drugz", exist_ok=True)

# Call the drugZ_analysis function for each pair of target_sample and reference_sample
for target_sample, reference_sample in zip(target_samples_list, reference_samples_list):
    args = Args()
    args.target_sample = target_sample
    args.reference_sample = reference_sample
    args.drugz_output_file = args.drugz_output_file + target_sample + "-" + reference_sample + ".txt"
    args.gRNA_outfile = args.gRNA_outfile + target_sample + "-" + reference_sample + ".txt"
    drugZ.drugz_analysis(args)
