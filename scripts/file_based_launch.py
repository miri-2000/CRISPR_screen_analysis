from start_program import CRISPRScreenAnalysis


class Args:
    input_file = r"file\path\to\the\crispr_screen_read_counts.txt"
    essential_genes = r"file\path\to\the\essential_genes.csv"
    non_essential_genes = r"file\path\to\the\non_essential_genes.csv"
    library_file = r"file\path\to\the\library_file.txt"
    target_samples = "treated_conditions_within_input_file_comma_separated"  # e.g. t1_2d_tr,t2_2d_tr
    reference_samples = "baseline_conditions_within_input_file_comma_separated"  # e.g. t1_2d_ut,t2_2d_ut
    unwanted_columns = "unwanted_columns_within_input_file_comma_separated"  # e.g. guide_mm1_mismatch1,mismatch1_
    unwanted_rows = "unwanted_rows_within_input_file_comma_separated"
    unwanted_row_substrings = "unwanted_row_substrings_within_input_file_comma_separated"
    threshold_reads = 0
    x_axis = "normZ"  # or "log2 fold-change"
    threshold_fdr = 0.25
    top = 15
    distribution_condition1 = "positive_control_sample"
    distribution_condition2 = "negative_control_sample"
    replicate_type = "biological"  # or "technical"
    working_dir = r"folder\path\to\the\crispr_screen_analysis"


if __name__ == "__main__":
    CRISPRScreenAnalysis(Args)
