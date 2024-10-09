from start_program import CRISPRScreenAnalysis
from pathlib import Path


class Args:
    input_file = Path(__file__).parents[1] / "example" / "5994_example_screen_read_counts.txt"
    essential_genes = Path(__file__).parents[1] / "example" / "essential_genes.csv"
    non_essential_genes = Path(__file__).parents[1] / "example" / "non_essential_genes.csv"
    library_file = Path(__file__).parents[1] / "example" / "library_file.txt"
    target_samples = "t7,t29"
    reference_samples = "t0,t7"
    unwanted_columns = "guide_mm1_mismatch1,mismatch1_,nohit_cols,guide_mm1_nohit"
    unwanted_rows = ""
    unwanted_row_substrings = ":mismatch"
    threshold_reads = 0
    x_axis = "normZ"  # or "log2 fold-change"
    threshold_fdr = 0.25
    top = 15
    distribution_condition1 = "t7"
    distribution_condition2 = "t0"
    replicate_type = "biological"  # or "technical"
    working_dir = Path(__file__).parents[1] / "example"


if __name__ == "__main__":
    CRISPRScreenAnalysis(Args)
