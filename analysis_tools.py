# analysis_tools: Helper script required for running R scripts and other python files.
# Last modified 19.11.2023
# ------------------------------------
import re
import os
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
    exe_file = get_exe(script_file)
    command = [exe_file, script_file]
    if additional_args is not None:
        command += additional_args

    # Set up a try-except statement to catch occurring errors while executing the R file
    try:
        result = subprocess.run(args=command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if not exe_file.endswith("Rscript.exe"):
            print(result.stderr)
    except subprocess.CalledProcessError as e:
        print("Command returned non-zero exit status:", e.returncode)
        print("Standard Output:", e.stdout)
        print("Standard Error:", e.stderr)
        raise SystemExit()


def get_exe(script_file):
    """
    This function is called to return the path of the R and Python exe files.

    :param script_file: File path of the script to be run

    :return: File path to the exe file
    """

    global exe_files

    if script_file.endswith("R"):
        program = "R"
        exe_name = "Rscript.exe"
        # Specify the path that should contain the R exe file
        path_to_search = r"C:\Users\Miriam\AppData\Local\Programs\R"
    else:
        exe_name = "python.exe"
        program = "Python"
        # Specify the path that should contain the Python exe file
        path_to_search = r"C:\Users\Miriam\PycharmProjects"

    if program not in exe_files:
        # Use the os.walk function to traverse the directory tree
        for dirpath, dirnames, filenames in os.walk(path_to_search):
            if exe_name in filenames:
                exe_files[program] = os.path.join(dirpath, exe_name)
                break
        else:
            raise FileNotFoundError(
                f"{exe_name} cannot be found on your computer. Please check if {program} is installed.")

    return exe_files[program]
