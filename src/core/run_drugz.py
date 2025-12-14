# run_drugz: This script is used to execute DrugZ and allows the
# modification of additional parameters for running DrugZ if wanted
# by the user.
# Last modified 14.12.2025
# ------------------------------------
from src.core import drugZ
import os
import sys
from pathlib import Path


class Args:
    def __init__(self):
        self.infile = list(Path(sys.argv[1]).glob("*_drugz-input.txt"))
        self.drugz_output = "./drugz/drugz_"
        self.gRNA_output = "./drugz/gRNA_"
        self.target_samples = sys.argv[2]
        self.reference_samples = sys.argv[3]
        self.remove_genes = None
        self.unpaired = True
        self.replicates = True
        self.pseudocount = 5
        self.half_window_size = 500


if __name__ == '__main__':
    args = Args()
    target_samples_list = args.target_samples.split(',')
    reference_samples_list = args.reference_samples.split(',')

    # Create a DrugZ output folder (if it does not exist yet)
    os.makedirs("./drugz", exist_ok=True)

    # Call the drugZ_analysis function for each pair of target_sample and reference_sample
    for target_sample, reference_sample in zip(target_samples_list, reference_samples_list):
        args.target_sample = target_sample
        args.reference_sample = reference_sample
        args.drugz_output_file = args.drugz_output + target_sample + "-" + reference_sample + ".txt"
        args.gRNA_output_file = args.gRNA_output + target_sample + "-" + reference_sample + ".txt"
        drugZ.drugz_analysis(args)
