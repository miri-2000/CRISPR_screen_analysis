import os
import pandas as pd
from abc import ABC, abstractmethod


class InputValidator(ABC):
    def __init__(self, args):
        self.working_dir = args.working_dir
        self.input_file = args.input_file
        self.dataset = os.path.basename(self.input_file).split('_')[0]
        self.target_samples = args.target_samples
        self.reference_samples = args.reference_samples
        self.essential_genes = args.essential_genes
        self.non_essential_genes = args.non_essential_genes
        self.library_file = args.library_file
        self.unwanted_columns = args.unwanted_columns
        self.unwanted_rows = args.unwanted_rows
        self.unwanted_row_substrings = args.unwanted_row_substrings
        self.threshold_reads = args.threshold_reads
        self.x_axis = args.x_axis
        self.threshold_fdr = args.threshold_fdr
        self.top = args.top
        self.distribution_condition1 = args.distribution_condition1
        self.distribution_condition2 = args.distribution_condition2
        self.replicate_type = args.replicate_type
        self.input_data = pd.read_csv(self.input_file, sep="\t")

    def validate_files(self):
        self.validate_file_path()
        self.validate_input_file()
        self.validate_gene_files()
        self.validate_library_file()

    def validate_file_path(self):
        file_path_vars = ["input_file", "essential_genes", "non_essential_genes", "library_file"]

        for file_path_var in file_path_vars:
            file_path = getattr(self, file_path_var, None)
            if not os.path.exists(file_path) or not (
                    str(file_path).endswith('.txt') or str(file_path).endswith('.csv')):
                self.abort(
                    f"The parameter '{file_path_var}' needs to contain a valid file path that leads to the specified "
                    f"file. The file must either be a CSV or TXT file.")

    def validate_input_file(self):

        if (len(self.input_data.iloc[:, 0].unique()) != len(self.input_data.iloc[:, 0]) or
                len(self.input_data.iloc[:, 1].unique()) != len(self.input_data.iloc[:, 1])):
            self.abort(
                f"The CRISPR screen input file requires the first column with the sgRNA names and the second column "
                f"with the sgRNA sequences to be unique.")

        if not (self.input_data.iloc[:, 0] == 'nohit_row').any():
            self.abort('The CRISPR screen input file requires a no-hit row.')

        # if not (self.input_data.iloc[:, 0] == re.compile(r'Non[-_. ]Targeting[-_. ]Control', re.IGNORECASE)).any():
        #     self.abort('The CRISPR screen input file requires a non-targeting control row.')

        if not self.input_data.columns.str.contains('guide_mm1').any():
            self.abort('The CRISPR screen input file requires a guide_mm1_ column.')

    def validate_gene_files(self):
        for gene_type in ["essential_genes", "non_essential_genes"]:
            gene_file = getattr(self, gene_type, None)
            genes = pd.read_csv(gene_file)

            input_data_split = self.input_data.iloc[:, 0].str.split('-', expand=True)[0]

            if not genes.iloc[:, 0].isin(input_data_split).all():
                self.abort(f'All genes mentioned in the "{gene_type}" file need to be present in the CRISPR screen '
                           f'input file.')

    def validate_library_file(self):
        library_data = pd.read_csv(self.library_file, sep="\t")

        if not pd.Index(['Target Gene Symbol', 'sgRNA Target Sequence']).isin(library_data.columns).all():
            self.abort('The library file requires the column "Target Gene Symbol" (holding the gene names) and'
                       '"sgRNA Target Sequence" (holding the sgRNA sequence) to be present in the CRISPR screen.')

    # def validate_samples(self):
    #     sample_vars = ['target_samples','reference_samples']
    #     for sample_var in sample_vars:
    #         sample = getattr(self, sample_var, None)
    #
    #         for sample_elem in sample.split(','):
    #             if sample_elem not in self.input_data.columns:
    #                 self.abort(f'The sample {sample_elem} from {sample_var} cannot be found in the CRISPR screen
    #                 input file. Please make sure that all mentioned samples represent columns in the input file.')

    def validate_int_fields(self, which='both'):
        int_field_vars = ['threshold_reads', 'top']
        fields_to_check = {
            "both": int_field_vars,
            "threshold_reads": [int_field_vars[0]],
            "top": [int_field_vars[1]]
        }  # Default to "both" if an invalid option is passed

        for int_field_var in fields_to_check[which]:
            int_field = getattr(self, int_field_var, None)

            if not isinstance(int_field, int):
                self.abort(f'The parameter "{int_field_var} needs to be a whole number.')

    def validate_threshold_fdr(self):
        if not isinstance(self.threshold_fdr, (int, float)) and not 0 <= self.threshold_fdr <= 1:
            self.abort(f'The parameter "threshold_fdr" needs to be a number between 0 and 1.')

    def validate_choice_fields(self):
        choice_fields = {'x_axis': ['normZ', 'log2 fold-change'], 'replicate_type': ['biological', 'technical']}
        for choice_field, choice_values in choice_fields.items():
            choice = getattr(self, choice_field, None)

            if choice not in choice_values:
                self.abort(f'The replicate type can only either be "{choice_values[0]}" or "{choice_values[1]}".')

    def validate_working_dir(self):
        if not os.path.isdir(self.working_dir):
            self.abort('The given working directory needs to be a directory. Please make sure you inserted a folder '
                       'location.')

    @abstractmethod
    def abort(self, message):
        pass
