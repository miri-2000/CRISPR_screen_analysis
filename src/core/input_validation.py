import os
import pandas as pd
from abc import ABC, abstractmethod
from src.core.file_input import read_file


class InputValidator(ABC):
    def __init__(self):
        pass

    def validate_files(self, input_file, essential_genes, non_essential_genes, library_file):

        if not self.validate_file_path(input_file, essential_genes, non_essential_genes, library_file):
            return False

        input_data = read_file(input_file)

        if not self.validate_input_file(input_data):
            return False

        if not self.validate_gene_files(input_data, essential_genes, non_essential_genes):
            return False

        if not self.validate_library_file(library_file):
            return False

        return True

    def validate_file_path(self, input_file, essential_genes, non_essential_genes, library_file):
        file_paths = {"Screen Result File": input_file, "Essential Genes File": essential_genes,
                      "Non-Essential Genes File": non_essential_genes, "Library File": library_file}
        error_vars = []

        for file_path_name, file_path in file_paths.items():
            if not os.path.exists(file_path) or not (
                    str(file_path).endswith('.txt') or str(file_path).endswith('.csv')):
                error_vars += [f"'{file_path_name}'"]
        if error_vars:
            self.abort(
                f"The following path/s is/are invalid: {', '.join(error_vars)}. The "
                f"file path/s must be accurate and lead to a text or csv file.")
            return False

        return True

    def validate_input_file(self, input_data):

        if (len(input_data.iloc[:, 0].unique()) != len(input_data.iloc[:, 0]) or
                len(input_data.iloc[:, 1].unique()) != len(input_data.iloc[:, 1])):
            self.abort(
                f"The CRISPR screen input file requires the first column with the sgRNA names and the second column "
                f"with the sgRNA sequences to be unique and not contain missing values (except for nohit row with empty"
                f" sgRNA sequence.")
            return False

        if not self.check_invalid_dtypes(input_data):
            return False

        if not (input_data.iloc[:, 0] == "nohit_row").any():
            self.abort("The CRISPR screen input file requires a no-hit row.")
            return False

        # if not (input_data.iloc[:, 0] == re.compile(r"Non[-_. ]Targeting[-_. ]Control", re.IGNORECASE)).any():
        #     self.abort("The CRISPR screen input file requires a non-targeting control row.")

        if not input_data.columns.str.contains("guide_mm1").any():
            self.abort("The CRISPR screen input file requires a guide_mm1_ column.")
            return False

        return True

    def check_invalid_dtypes(self, input_data):
        # Check if either of the first two columns contain only strings
        check_strings = dict()
        for col in range(2):
            check_strings[col] = input_data.iloc[:, col].apply(lambda x: isinstance(x, str) or pd.isna(x)).all()

        if not check_strings[0] and not check_strings[1]:
            self.abort(
                "The sgRNA name and sequence columns (column 1 and 2) in the CRISPR screen input file should "
                "only contain strings.")
            return False
        elif not check_strings[0]:
            self.abort(
                "The sgRNA name column (column 1) in the CRISPR screen input file should only contain strings.")
            return False
        elif not check_strings[1]:
            self.abort(
                "The sgRNA sequence column (column 2) in the CRISPR screen input file should only contain strings.")
            return False

        return True

    def validate_gene_files(self, input_data, essential_genes, non_essential_genes):
        gene_files = {"Essential Genes File": essential_genes,
                      "Non-Essential Genes File": non_essential_genes}
        invalid_files = []

        for gene_file_name, gene_file in gene_files.items():
            genes = read_file(gene_file)

            input_data_split = input_data.iloc[:, 0].str.split("-", expand=True)[0]

            if not genes.iloc[:, 0].isin(input_data_split).all():
                invalid_files += [gene_file_name]
        if invalid_files:
            if len(invalid_files) == 1:
                message = (
                    f"All genes mentioned in the '{invalid_files[0]}' file need to be present in the CRISPR screen "
                    f"input file.")
            else:
                message = (
                    f"All genes mentioned in the '{invalid_files[0]}' and '{invalid_files[1]}' files need to be present"
                    f" in the CRISPR screen input file.")
            self.abort(message)
            return False

        return True

    def validate_library_file(self, library_file):
        library_data = read_file(library_file)

        if not pd.Index(["Target Gene Symbol", "sgRNA Target Sequence"]).isin(library_data.columns).all():
            self.abort("The library file requires the column 'Target Gene Symbol' (holding the gene names) and"
                       "'sgRNA Target Sequence' (holding the sgRNA sequence) to be present in the CRISPR screen.")
            return False

        return True

    # def validate_samples(self,input_data):
    #     sample_vars = ["target_samples","reference_samples"]
    #     for sample_var in sample_vars:
    #         sample = getattr(self, sample_var, None)
    #
    #         for sample_elem in sample.split(","):
    #             if sample_elem not in input_data.columns:
    #                 self.abort(f"The sample {sample_elem} from {sample_var} cannot be found in the CRISPR screen
    #                 input file. Please make sure that all mentioned samples represent columns in the input file.")

    def validate_int_fields(self, threshold_reads=None, top=None):
        int_field_vars = {"Minimum required sum of reads/guide": threshold_reads, "Number of hits per plot": top}
        fields_to_check = {var_name: var for var_name, var in int_field_vars.items() if var is not None}

        for field_name, field_var in fields_to_check.items():
            if isinstance(field_var, float):
                self.abort(f"The '{field_name}' needs to be a whole number.")
                return False
            try:
                field_var = int(field_var)
            except ValueError:
                self.abort(f"The '{field_name}' needs to be a whole number.")
                return False

        return True

    def validate_threshold_fdr(self, threshold_fdr):
        try:
            threshold_fdr = float(threshold_fdr)
        except ValueError:
            self.abort(f"The 'Significance threshold' needs to be a number between 0 and 1.")
            return False

        if not 0 <= threshold_fdr <= 1:
            self.abort(f"The 'Significance threshold' needs to be a number between 0 and 1.")
            return False

        return True

    def validate_choice_fields(self, x_axis, replicate_type):
        choice_fields = {"x_axis": [x_axis, "normZ", "log2 fold-change"],
                         "replicate_type": [replicate_type, "biological", "technical"]}
        for choice_field, choice_values in choice_fields.items():
            choice = choice_values[0]

            if choice not in choice_values[1:]:
                self.abort(f"The {choice_field} can only either be '{choice_values[1]}' or '{choice_values[2]}'.")
                return False

        return True

    def validate_working_dir(self, working_dir):
        if not os.path.isdir(working_dir):
            self.abort("The given working directory needs to be a directory. Please make sure you inserted a folder "
                       "location.")
            return False

        return True

    @abstractmethod
    def abort(self, message):
        pass
