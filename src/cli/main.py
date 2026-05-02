import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
from ..core.start_program import CRISPRScreenAnalysis


@dataclass
class RunConfig:
    input_file: Path
    output_dir: Path
    target_samples: List[str]
    reference_samples: List[str]
    essential_genes: Path
    non_essential_genes: Path
    library_file: Path
    unwanted_columns: List[str] = field(default_factory=list)
    unwanted_rows: List[str] = field(default_factory=list)
    unwanted_row_substrings: List[str] = field(default_factory=list)
    threshold_reads: int = 30
    x_axis: str = 'normZ'
    threshold_fdr: float = 0.05
    top: int = 10
    distribution_condition1: str = ""
    distribution_condition2: str = ""
    replicate_type: str = 'biological'


def run_pipeline(cfg: RunConfig):
    analysis = CRISPRScreenAnalysis()
    analysis.run_analysis(
        working_dir=str(cfg.output_dir),
        input_file=str(cfg.input_file),
        target_samples=cfg.target_samples,
        reference_samples=cfg.reference_samples,
        essential_genes=str(cfg.essential_genes),
        non_essential_genes=str(cfg.non_essential_genes),
        library_file=str(cfg.library_file),
        unwanted_columns=cfg.unwanted_columns,
        unwanted_rows=cfg.unwanted_rows,
        unwanted_row_substrings=cfg.unwanted_row_substrings,
        threshold_reads=cfg.threshold_reads,
        x_axis=cfg.x_axis,
        threshold_fdr=cfg.threshold_fdr,
        top=cfg.top,
        distribution_condition1=cfg.distribution_condition1,
        distribution_condition2=cfg.distribution_condition2,
        replicate_type=cfg.replicate_type
    )


def main():
    parser = argparse.ArgumentParser(description='CRISPR Screen Analysis Pipeline')
    parser.add_argument("--input", required=True, help="Input read counts file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--target-samples", required=True, nargs='+', 
                        help="Target sample names (space-separated)")
    parser.add_argument("--reference-samples", required=True, nargs='+',
                        help="Reference sample names (space-separated)")
    parser.add_argument("--essential-genes", required=True,
                        help="Path to essential genes CSV file")
    parser.add_argument("--non-essential-genes", required=True,
                        help="Path to non-essential genes CSV file")
    parser.add_argument("--library-file", required=True,
                        help="Path to library file")
    parser.add_argument("--threshold-reads", type=int, default=30,
                        help="Minimum reads per guide (default: 30)")
    parser.add_argument("--threshold-fdr", type=float, default=0.05,
                        help="Maximum FDR threshold (default: 0.05)")
    parser.add_argument("--top", type=int, default=10,
                        help="Number of top genes to display (default: 10)")
    parser.add_argument("--x-axis", default='normZ',
                        help="Metric for x-axis (default: normZ)")
    parser.add_argument("--distribution-condition1", default="",
                        help="Treated condition for distribution plot")
    parser.add_argument("--distribution-condition2", default="",
                        help="Baseline condition for distribution plot")
    parser.add_argument("--replicate-type", default='biological',
                        choices=['biological', 'technical'],
                        help="Replicate type (default: biological)")

    args = parser.parse_args()

    cfg = RunConfig(
        input_file=Path(args.input),
        output_dir=Path(args.output),
        target_samples=args.target_samples,
        reference_samples=args.reference_samples,
        essential_genes=Path(args.essential_genes),
        non_essential_genes=Path(args.non_essential_genes),
        library_file=Path(args.library_file),
        threshold_reads=args.threshold_reads,
        threshold_fdr=args.threshold_fdr,
        top=args.top,
        x_axis=args.x_axis,
        distribution_condition1=args.distribution_condition1,
        distribution_condition2=args.distribution_condition2,
        replicate_type=args.replicate_type,
    )

    run_pipeline(cfg)


if __name__ == "__main__":
    main()
