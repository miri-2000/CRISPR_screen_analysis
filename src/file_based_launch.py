from src.core.start_program import CRISPRScreenAnalysis
from pathlib import Path
from src.core.input_validation_cl import InputValidatorCL


class Args:
    def __init__(self):
        self.input_file = Path(__file__).parents[1] / "example" / "5994_example_screen_read_counts.txt"
        self.essential_genes = Path(__file__).parents[1] / "example" / "essential_genes.csv"
        self.non_essential_genes = Path(__file__).parents[1] / "example" / "non_essential_genes.csv"
        self.library_file = Path(__file__).parents[1] / "example" / "library_file.txt"
        self.target_samples = "t7,t29"
        self.reference_samples = "t0,t7"
        self.unwanted_columns = "guide_mm1_mismatch1,mismatch1_,nohit_cols,guide_mm1_nohit"
        self.unwanted_rows = ""
        self.unwanted_row_substrings = ":mismatch"
        self.threshold_reads = 0
        self.x_axis = "normZ"  # or "log2 fold-change"
        self.threshold_fdr = 0.25
        self.top = 15
        self.distribution_condition1 = "t7"
        self.distribution_condition2 = "t0"
        self.replicate_type = "biological"  # or "technical"
        self.working_dir = Path(__file__).parents[1] / "temp_out"


if __name__ == "__main__":
    args = Args()

    validator = InputValidatorCL()
    validator.validate(args.input_file, args.essential_genes, args.non_essential_genes, args.library_file,
                       args.threshold_reads, args.top, args.threshold_fdr, args.x_axis, args.replicate_type,
                       args.working_dir)

    analysis = CRISPRScreenAnalysis()
    analysis.run_analysis(args.working_dir, args.input_file, args.target_samples, args.reference_samples,
                          args.essential_genes, args.non_essential_genes, args.library_file, args.unwanted_columns,
                          args.unwanted_rows, args.unwanted_row_substrings, args.threshold_reads, args.x_axis,
                          args.threshold_fdr, args.top, args.distribution_condition1, args.distribution_condition2,
                          args.replicate_type)
