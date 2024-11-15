import pandas as pd
from pathlib import Path


def read_file(file_path):
    if isinstance(file_path, Path):
        if file_path.suffix == ".csv":
            return get_csv_file_content(file_path)
        elif file_path.suffix == ".txt":
            return get_txt_file_content(file_path)
    elif isinstance(file_path, str):
        if file_path.endswith(".csv"):
            return get_csv_file_content(file_path)
        elif file_path.endswith(".txt"):
            return get_txt_file_content(file_path)


def get_csv_file_content(file_path):
    return pd.read_csv(file_path, header=0)


def get_txt_file_content(file_path):
    return pd.read_csv(file_path, sep='\t', header=0)
