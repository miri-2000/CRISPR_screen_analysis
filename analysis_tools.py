# analysis_tools: Helper script required for running R scripts and other python files.
# Last modified 02.10.2024
# ------------------------------------
import re
import os
import sys

import pandas as pd
import subprocess

# Initialize the exe_file variable
exe_files = {}


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
    pattern = re.compile(r"Non-Targeting(_| )Control")

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


def run_script(script_file, additional_args=None):
    """
    This function is used to execute additional Python and R files using the module "subprocess".
    The exe-file is automatically identified so that the given script file can be executed and the
    error/logging messages are returned.

    :param script_file: Filepath of the script to be run
    :param additional_args: Required input-arguments to run the script

    :return: Error messages, if the script returns a non-zero exit code
    """

    # Identify where the Rscript.exe file is located on the pc
    if script_file.suffix == ".R":
        command = ["Rscript", script_file]
    else:
        command = [sys.executable, script_file]
    if additional_args is not None:
        command += additional_args
    print(command)

    # Set up a try-except statement to condcatch occurring errors while executing the R file
    try:
        result = subprocess.run(args=command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if not script_file.suffix == ".R":
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("Command returned non-zero exit status:", e.returncode)
        print("Standard Output:", e.stdout)
        print("Standard Error:", e.stderr)
        raise SystemExit()
