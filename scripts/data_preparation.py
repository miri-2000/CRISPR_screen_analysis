# data_preparation: Script for preparing the given dataset for the hit identification
# Last modified 02.10.2024
# ------------------------------------
# Purpose of the program
# 1) remove unwanted rows (here, rows with 'multi:mismatch1'; cannot be mapped to specific Gene)
# 2) remove unwanted row substrings (here, mismatch indication (e.g. from USP17L7-1:mismatch1_USP17L19-1 to USP17L7-1))
# 3) remove columns that contain "guide_mm1_mismatch1" or "guide_mm1_nohit"
# 4) make a summary table showing the sum of the guides and guide_mm1s + the nohit_guide
# 5) store the cleaned table in a csv file
# 7) Call the R script "R_analysis_1" to normalize the read counts and create raw and normalized plots
# ------------------------------------
from pathlib import Path
import pandas as pd
import re
from analysis_tools import assign_type, run_script
import logging as log

log.basicConfig(level=log.DEBUG)
log_ = log.getLogger(__name__)


class DataPreparation:
    """
    Perform pre-processing to establish whether the quality of the screen output is
    good enough to continue with the hit identification.
    """

    def __init__(self, input_file, essential_genes, non_essential_genes, library_file, unwanted_columns, unwanted_rows,
                 unwanted_row_substrings, summary_file, threshold_reads, output_file, dataset, replicate_type,
                 target_samples, reference_samples, distribution_condition1, distribution_condition2):

        self.input_file = input_file
        self.essential_genes = essential_genes
        self.non_essential_genes = non_essential_genes
        self.library_file = library_file
        self.unwanted_columns = unwanted_columns
        self.unwanted_rows = unwanted_rows
        self.unwanted_row_substrings = unwanted_row_substrings
        self.summary_file = summary_file
        self.threshold_reads = threshold_reads
        self.output_file = output_file
        self.dataset = dataset
        self.replicate_type = replicate_type
        self.target_samples = target_samples
        self.reference_samples = reference_samples
        self.distribution_condition1 = distribution_condition1
        self.distribution_condition2 = distribution_condition2

    def prepare_data(self):

        log_.debug("Initiating data preparation")

        # Get the data
        reads, essential_genes, non_essential_genes, library = self.get_data()

        # Pre-process the data
        reads, conditions = self.preprocess_data(reads)

        # Compare the cleaned data file with the given library file
        self.compare_with_library(reads, library, essential_genes, non_essential_genes)

        # Store the cleaned data in a csv file
        log_.debug(f"\nStoring the resulting dataframe in {self.output_file}")
        reads.to_csv(self.output_file, sep=';', index=False)

        # Normalize data and creating quality control plots
        self.normalize_data(conditions)

        log_.debug("Finished Data preparation")

    def get_data(self):
        # Read the provided input files
        reads, essential_genes, non_essential_genes, library = self.read_files(self.input_file, self.essential_genes,
                                                                               self.non_essential_genes,
                                                                               self.library_file)

        # Check the input file for its correctness
        self.check_input(reads)

        # Clean up the input file, so it only contains required columns and rows
        reads = self.clean_up(reads, essential_genes, non_essential_genes, self.unwanted_columns, self.unwanted_rows,
                              self.unwanted_row_substrings)
        return reads, essential_genes, non_essential_genes, library

    @staticmethod
    def read_files(data_file, essential_genes_file, non_essential_genes_file, library_file):
        """
        Read the input files

        :param data_file: File with the read counts for each gRNA
        :param essential_genes_file: File with the essential genes data
        :param non_essential_genes_file: File with the non-essential genes data
        :param library_file: File with the library data
        :return: Tuple with the input files stored as dataframes
        """

        read_data = pd.read_csv(data_file, sep="\t", header=0)
        essential_genes = pd.read_csv(essential_genes_file)
        non_essential_genes = pd.read_csv(non_essential_genes_file)
        library = pd.read_csv(library_file, sep="\t")

        return read_data, essential_genes, non_essential_genes, library

    @staticmethod
    def check_input(read_data):
        """
        Check the input data for errors that prevent the data preparation from working.

        :param read_data: Dataset containing read counts for each gRNA
        :return: Error message in case the input is not correct
        """

        # Check if either of the first two columns contain only strings and do not have empty rows (except for the
        # nohit-row)
        for col in range(2):
            nan_count = read_data.iloc[:, col].isna().sum()
            if nan_count > 1:
                raise ValueError(
                    f"Column {col + 1} has multiple empty rows. Please make sure that every row has a value.")

            is_all_strings = read_data.iloc[:, col].apply(lambda x: isinstance(x, str) or pd.isna(x)).all()

            if not is_all_strings:
                raise ValueError(f"Column {col + 1} should contain only strings")

        # Control if the columns with the sgRNA names and sequences from the table are both unique
        log_.debug("Checking the uniqueness of the sgRNA names and sequences")
        if len(read_data) != len(set(read_data.iloc[:, 0])) or len(read_data) != len(set(read_data.iloc[:, 1])):
            raise ValueError(
                "Both columns with sgRNA names and sequences need to be unique. If this is true, please check whether "
                "these two columns are the first two columns of the input data file (as expected).")

        # Check if there is a nohit_row
        log_.debug("Identifying the nohit_row")
        nohit_row = read_data[read_data.iloc[:, 0] == 'nohit_row']
        if len(nohit_row) == 0:
            raise ValueError("There is no nohit_row. Please check your dataset for a row that is called 'nohit_row'.")

    def clean_up(self, read_data, essential_genes, non_essential_genes, unwanted_columns=None, unwanted_rows=None,
                 unwanted_row_substrings=None):
        """
        Remove unwanted rows and columns and prepare the dataset for
        sequential processing

        :param read_data: Dataset containing read counts for each gRNA
        :param essential_genes: list of genes that are considered essential
        :param non_essential_genes: list of genes that are considered non-essential
        :param unwanted_columns: List of substrings used to match and delete columns (Optional)
        :param unwanted_rows: List of substrings used to match and delete rows (Optional)
        :param unwanted_row_substrings: List of substrings used to match and remove from sgRNA name (Optional)

        :return: read_data: Dataset with the created changes
        """

        # Rename the first column to "sgRNA" to avoid inconsistent column names
        log_.debug("Renaming first column to 'sgRNA'")
        read_data.rename(columns={read_data.columns[0]: "sgRNA"}, inplace=True)

        if unwanted_columns:
            # Remove unwanted columns
            read_data = self.remove_unwanted_columns(read_data, unwanted_columns)

        if unwanted_rows:
            # Remove rows with multi-matches (contain "multi:mismatch1" in column "name")
            read_data = self.remove_unwanted_rows(read_data, unwanted_rows)

        if unwanted_row_substrings:
            # Remove unwanted substrings from the gRNA column
            read_data = self.remove_unwanted_row_substrings(read_data, unwanted_row_substrings)

        # The column names are converted to a common naming convention (stated in the documentation)
        # to avoid inconsistency in the output files between different screens.
        # e.g. T1_3D_to_2D_Rep3:TGGTCA is changed to t1_3d2d_r3
        read_data = self.standardize_column_names(read_data)

        # Create a new column called 'Gene' that contains the names of the genes
        # this is done by removing the extensions from the sgRNA names
        read_data = self.create_gene_column(read_data)

        # Assign a type to each gRNA
        read_data = self.assign_gene_type(read_data, essential_genes, non_essential_genes)

        # Return the cleaned dataframe
        return read_data

    @staticmethod
    def remove_unwanted_columns(read_data, unwanted_columns):
        """
        Remove unwanted columns from the input dataset

        :param read_data: Dataset containing read counts for each gRNA
        :param unwanted_columns: List of substrings used to match and delete columns (Optional)
        :return: Cleaned dataset without unwanted columns
        """

        log_.debug("Removing unwanted columns")
        log_.debug(f"Columns to be removed: {unwanted_columns}")
        # For this dataset: Remaining columns either contain the information for each condition or belong to "guide_mm1"
        length_before = len(read_data.columns)
        read_data = read_data.drop(
            columns=read_data.filter(regex=unwanted_columns.replace(",", "|")))
        if length_before != len(read_data.columns):
            log_.debug(
                f"Number of columns before: {length_before}, after: {len(read_data.columns)}, diff: {length_before - len(read_data.columns)}")
        else:
            log_.debug("No unwanted columns identified")

        return read_data

    @staticmethod
    def remove_unwanted_rows(read_data, unwanted_rows):
        """
        Remove unwanted rows from the input dataset

        :param read_data: Dataset containing read counts for each gRNA
        :param unwanted_rows: List of substrings used to match and delete rows (Optional)
        :return: Cleaned dataset without unwanted rows
        """

        log_.debug("Removing unwanted rows")
        log_.debug(
            f"Rows containing (one of) the following in the sgRNA name will be removed: {unwanted_rows}")
        rows_to_be_removed = read_data["sgRNA"].str.contains(unwanted_rows.replace(",", "|"))
        if any(rows_to_be_removed):
            log_.debug(
                f"Number of rows before: {len(read_data)}, after: {len(read_data[~rows_to_be_removed])}, diff: {len(read_data[rows_to_be_removed])}")
            read_data = read_data[~rows_to_be_removed]
        else:
            log_.debug("No unwanted rows identified")

        return read_data

    @staticmethod
    def remove_unwanted_row_substrings(read_data, unwanted_row_substrings):
        """
        Remove unwanted row substrings from the input dataset

        :param read_data: Dataset containing read counts for each gRNA
        :param unwanted_row_substrings: List of substrings used to match and remove from sgRNA name (Optional)
        :return: Cleaned dataset without unwanted row substrings
        """

        log_.debug("Removing unwanted substrings from rows of column sgRNA")
        log_.debug(
            f"The following substrings will be removed from the sgRNA name: {unwanted_row_substrings}")
        rows_to_be_edited = read_data[read_data["sgRNA"].str.contains(unwanted_row_substrings.replace(",", "|"))]
        if len(rows_to_be_edited) > 0:
            for substring in unwanted_row_substrings.split(','):
                read_data['sgRNA'] = read_data['sgRNA'].str.replace(substring + '.*', '', regex=True)
            log_.debug(f"Number of edited rows: {len(rows_to_be_edited)}")
        else:
            log_.debug("No unwanted row substrings identified")

        return read_data

    def preprocess_data(self, reads):
        # Store all column names representing the conditions and all columns that have "guide_mm1" in their name
        conditions, data_columns = self.get_conditions(reads)

        # Remove the nohit row
        reads, nohit_guide = self.remove_nohit(reads, conditions)

        # Create a summary file
        self.create_summary(reads, conditions, nohit_guide, self.summary_file)

        # Remove columns containing "guide_mm1"
        reads = self.remove_guide_mm1_columns(reads)

        # Remove the rows that have a total sum of counts below a certain threshold
        reads = self.remove_invalid_guides(reads, self.threshold_reads)

        return reads, conditions

    def standardize_column_names(self, read_data):
        """
        Convert column names to a common naming convention

        :param read_data: Dataset containing read counts for each gRNA
        :return: Cleaned dataset with normalized column names
        """
        log_.debug("Modifying column names")
        log_.debug("Changing names like 'T1_3D_to_2D_Rep3:TGGTCA' to 't1_3d2d_r3'")
        numeric_data = read_data.select_dtypes(include='number')
        annotations = read_data.select_dtypes(exclude='number')
        
        numeric_data.columns = (numeric_data.columns.str.replace(':.*', '', regex=True)
                                .str.lower()
                                .str.replace("untreated", "ut")
                                .str.replace("treated", "tr")
                                .str.replace("_rep", "_r")
                                .str.replace("3d_to_2d", "3d2d")
                                .str.replace(r"\(|\)", ""))

        # In case the column name structures are not equal across all columns:
        # e.g. change T2_Treated_(3D)_Rep2:ATTGGC to t2_3d_tr_r2
        new_columns = []
        for column in numeric_data.columns:
            if "(" in column:
                new_columns += [re.sub(r'_(tr|ut)_\((.*?)\)_(r\d)', r'_\2_\1_\3', column)]
            else:
                new_columns += [column]
        numeric_data.columns = new_columns

        read_data = pd.concat([annotations, numeric_data], axis=1)

        return read_data

    def create_gene_column(self, read_data):
        """
        Create a gene column from read_data

        :param read_data: Dataset containing read counts for each gRNA
        :return: Dataset with an added gene column
        """

        log_.debug("Creating a column for the gene names")
        pattern = r"-[0-9]*$"
        read_data.insert(1, "Gene", [re.sub(pattern, "", gene) for gene in read_data["sgRNA"]])

        return read_data

    def assign_gene_type(self, read_data, essential_genes, non_essential_genes):
        """
        Create a column with the assigned data type for each gene

        :param read_data: Dataset containing read counts for each gRNA
        :param essential_genes: list of genes that are considered essential
        :param non_essential_genes: list of genes that are considered non-essential

        :return: Dataset with an added column for the gene type
        """

        log_.debug("Creating a column with the assigned data type for each gene")

        read_data.insert(2, "type", assign_type(read_data['Gene'], essential_genes, non_essential_genes))

        log_.debug(
            f"\nDataset composition:\nNumber of essential genes (p) = {len(set(read_data['Gene'][read_data['type'] == 'p']))}, "
            f"Number of sgRNAs = {len(read_data['Gene'][read_data['type'] == 'p'])}\n"
            f"Number of non-essential genes (n) = {len(set(read_data['Gene'][read_data['type'] == 'n']))}, "
            f"Number of sgRNAs = {len(read_data['Gene'][read_data['type'] == 'n'])}\n"
            f"Number of non-targeting controls (o) = {len(read_data['Gene'][read_data['type'] == 'o'])}\n"
            f"Number of undefined genes (x) = {len(set(read_data['Gene'][read_data['type'] == 'x']))}, "
            f"Number of sgRNAs = {len(read_data['Gene'][read_data['type'] == 'x'])}\n")

        return read_data

    def get_conditions(self, read_data):
        """
        Select columns that are used for creating the counts_summary file
        :param read_data: Dataset containing read counts for each gRNA

        :return: Tuple with the selected conditions and the data columns
        """

        log_.debug("Selecting columns that are used for creating the counts_summary file")
        data_columns = read_data.select_dtypes(include='number').columns
        guide_mm1 = []
        conditions = [col for col in data_columns if not re.search("guide_mm1", col) or guide_mm1.append(col)]

        return conditions, data_columns

    def remove_nohit(self, read_data, conditions):
        """
        Remove the nohit row from the dataset

        :param read_data: Dataset containing read counts for each gRNA
        :param conditions: Selected conditions
        :return: Tuple containing dataset without nohit row and the nohit row
        """

        log_.debug("Selecting the nohit row and removing it from the table")
        nohit_guide = read_data[read_data.iloc[:, 0] == 'nohit_row'][conditions]
        read_data = read_data.drop(nohit_guide.index)

        return read_data, nohit_guide

    def remove_guide_mm1_columns(self, read_data):
        """
        Remove the guide_mm1 columns from the dataset

        :param read_data: Dataset containing read counts for each gRNA
        :return: Dataset with guide_mm1 columns removed
        """
        log_.debug(f"Removing guide_mm1 columns")
        read_data = read_data.drop(columns=read_data.filter(regex="guide_mm1"))

        return read_data

    def remove_invalid_guides(self, read_data, threshold_reads):
        """
        Remove invalid guides from the dataset

        :param read_data: Dataset containing read counts for each gRNA
        :param threshold_reads: Minimum number of total reads/guide so that the guide will not be discarded

        :return: Dataset with invalid guides removed
        """

        if threshold_reads == 0:
            log_.debug(f"Removing sgRNAs with a sum over all columns of {threshold_reads}")
        else:
            log_.debug(f"Removing sgRNAs with a sum over all columns below {threshold_reads}")
        # Calculate row sums for the data columns
        row_sums = read_data.sum(numeric_only=True, axis=1)
        # Check if the sum of the total and filtered rows are the same
        if len(read_data) != len(read_data[row_sums > threshold_reads]):
            log_.debug(
                f"Number of rows before filtering: {len(read_data)}, after filtering: {len(read_data[row_sums > threshold_reads])},"
                f" difference: {len(read_data[row_sums <= threshold_reads])}")
            read_data = read_data[row_sums > threshold_reads]
        else:
            log_.debug(f"No sgRNAs count sums below {threshold_reads} identified")

        return read_data

    def create_summary(self, read_data, conditions, nohit_guide, summary_file):
        """
        This function creates a summary file that show the sums of reads per condition.
        The program makes three different classes (condition, condition with 1 mismatch and
        no match to any condition) and calculates the sum of reads for each of these classes per condition.

        :param read_data: Cleaned data file with read counts per condition
        :param conditions: A list of all conditions from the screen
        :param nohit_guide: The row with the reads per condition that could not be mapped to any guide (nohit_row)
        :param summary_file: Name of the output file
        :return:
        """

        log_.debug("Creating a data summary file")

        # Calculate the sum for all columns
        log_.debug("Calculating sum of each column")
        data_sum = read_data.sum(numeric_only=True, axis=0)

        # Store the sum values of the condition columns
        sum_conditions = pd.DataFrame(data_sum.loc[data_sum.index.isin(conditions)])

        # Store the sum values of the guide_mm1 columns
        sum_guide_mm1 = pd.DataFrame(data_sum.loc[~data_sum.index.isin(conditions)]).T

        # Remove the "guide_mm1" notation from the column names
        log_.debug("Aligning column names from condition- and guide_mm1- columns")
        sum_guide_mm1.columns = sum_guide_mm1.columns.str.replace("guide_mm1_", '')

        # Combine the data frames
        # first combine the values vertically to change the column names
        # and then transpose (allows creating correct row indices)
        log_.debug("Creating one dataframe with column sums from conditions and guide_mm1 plus the nohit_row")
        sumDf = pd.concat([sum_conditions, sum_guide_mm1.T, nohit_guide.T], axis=1)
        sumDf.columns = ["sum_conditions", "sum_guide_mm1", "nohit_guide"]
        sumDf = sumDf.T

        # Calculate the sum for each row and column
        log_.debug("Calculating the sum for each row and column of the new dataframe")
        sumDf['row_sum'] = sumDf.sum(axis=1)
        sumDf.loc["col_sum", :] = sumDf.sum(axis=0)

        # Write the result to a CSV file
        log_.debug(f"Storing the final table in {summary_file}")
        sumDf.to_csv(summary_file, index=True, sep=";")

    def compare_with_library(self, data, library, essential_genes, non_essential_genes):
        """
        This function compares the library file with the read counts file and gives out
        a summary that states how many guides had reads from all the guides that are
        contained in the original screen library.

        :param data: Cleaned dataset with read counts per condition
        :param library: Library file stating all guides that are contained in the screen library
        :param essential_genes: File that lists all essential control genes from the screen
        :param non_essential_genes: File that lists all non-essential control genes from the screen

        :return:
        """

        # Select the required columns for the comparison from the library file
        lib = library[["Target Gene Symbol", "sgRNA Target Sequence"]]

        # Assign the type (target, essential, non-essential, non-targeting) to each guide in the library
        lib_copy = lib.copy()
        lib_copy.loc[:, "type"] = assign_type(lib["Target Gene Symbol"], essential_genes, non_essential_genes)

        # Calculate the ratio of the number of all sgRNAs contained in the read counts file compared to the number of
        # sgRNAs in the library for each type
        log_.debug("\nLibray composition")
        for gene_type in ["p", "n", "o", "x"]:
            if len(lib_copy[lib_copy['type'] == gene_type]) != 0:
                ratio = round(100 * len(data[data['type'] == gene_type]) / len(lib_copy[lib_copy['type'] == gene_type]),
                              1)
            else:
                ratio = "-Inf"
            log_.debug(
                f"Number of type {gene_type} sgRNAs in library: {len(lib_copy[lib_copy['type'] == gene_type])}, in dataset: {len(data[data['type'] == gene_type])}, in percentage: {ratio} %")

    def normalize_data(self, conditions):
        # Remove the replicate notation ("_r1") from the condition names
        conditions_new = [condition.rsplit('_', 1)[0] for condition in conditions]

        log_.debug("Normalizing data and creating quality control plots\n")
        run_script(Path(__file__).parent / "R_analysis_1.R",
                   additional_args=[self.output_file, self.dataset, ",".join(conditions_new), self.replicate_type,
                                    self.target_samples, self.reference_samples, self.distribution_condition1,
                                    self.distribution_condition2])
