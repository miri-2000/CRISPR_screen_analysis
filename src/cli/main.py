import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--threads", type=int, default=1)

    args = parser.parse_args()

    cfg = RunConfig(
        input_file=Path(args.input),
        output_dir=Path(args.output),
        alpha=args.alpha,
        threads=args.threads,
    )

    run_pipeline(cfg)

if __name__ == "__main__":
    main()
