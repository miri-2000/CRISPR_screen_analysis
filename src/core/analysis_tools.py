# analysis_tools: Helper script required for running R src and other python files.
# Last modified 14.12.2025
# ------------------------------------
import re
from shutil import which
from pathlib import Path
import sys
import pandas as pd
import subprocess


def assign_type(genes, essential_list, non_essential_list):
    """
    This function assigns the type (target, essential, non-essential, non-targeting) of
    the targeted genes to each guide using the given lists of essential and non-essential genes.

    :param genes: Pandas dataframe column with the gene names
    :param essential_list: Pandas dataframe column with essential genes
    :param non_essential_list: Pandas dataframe column with non-essential genes
    :return: types: The types of each gene in the genes variable
    """
    # Specify how non-targeting control rows look like
    pattern = re.compile(r"Non[-_. ]Targeting[-_. ]Control", re.IGNORECASE)

    # Create a mask for essential and non-essential genes
    is_essential = genes.isin(essential_list.iloc[:, 0])
    is_non_essential = genes.isin(non_essential_list.iloc[:, 0])

    # Create a mask for non-targeting control genes
    is_non_targeting_control = genes.str.match(pattern)

    # Create a series for the new column values
    types = pd.Series('x', index=genes.index)

    # Use masks to assign types
    types[is_essential] = 'p'
    types[is_non_essential] = 'n'
    types[is_non_targeting_control] = 'o'

    return types

INTERPRETERS = {
    ".R": lambda: resolve_executable("Rscript"),
    ".py": lambda: sys.executable,
}

def resolve_executable(name: str) -> str:
    """
    Resolve an executable using the system PATH.
    """
    exe = which(name)
    if exe is None:
        raise FileNotFoundError(f"Required executable '{name}' not found in PATH")
    return exe

def run_script(script_file: Path, additional_args=None):
    """
    Execute a script using the appropriate interpreter.
    """

    script_file = Path(script_file)
    additional_args = additional_args or []

    try:
        interpreter = INTERPRETERS[script_file.suffix]()
    except KeyError:
        raise ValueError(f"Unsupported script type: {script_file.suffix}")

    command = [
        interpreter,
        str(script_file),
        *map(str, additional_args)
    ]

    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True
        )

        if result.stderr:
            print(result.stderr)

        return result.stdout

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Command failed with exit code {e.returncode}\n"
            f"STDOUT:\n{e.stdout}\n"
            f"STDERR:\n{e.stderr}"
        ) from e

