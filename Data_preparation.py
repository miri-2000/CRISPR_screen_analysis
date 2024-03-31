# Data_preparation: Script for preparing the given dataset for the hit identification
# Last modified 19.11.2023
# ------------------------------------
# Purpose of the program
# 1) remove unwanted rows (here, rows with 'multi:mismatch1'; cannot be mapped to specific Gene)
# 2) remove unwanted row substrings (here, mismatch indication (e.g. from USP17L7-1:mismatch1_USP17L19-1 to USP17L7-1))
# 3) remove columns that contain "guide_mm1_mismatch1" or "guide_mm1_nohit"
# 4) make a summary table showing the sum of the guides and guide_mm1s + the nohit_guide
# 5) store the cleaned table in a csv file
# 7) Call the R script "R_analysis_1" to normalize the read counts and create raw and normalized plots
# ------------------------------------

import pandas as pd
import re
from analysis_tools import assign_type, run_script
import logging as log

log.basicConfig(level=log.INFO)
log_ = log.getLogger(__name__)


def check_input(data):
    """
    Checks the input data file for errors that prevent the data preparation from working.

    :param data: input data file
    :return: Error message in case the input is not correct.
    """

    # Check if either of the first two columns contain only strings and do not have empty rows (except for the nohit-row)
    for col in range(2):
        nan_count = data.iloc[:, col].isna().sum()
        if nan_count > 1:
            raise ValueError(f"Column {col + 1} has multiple empty rows. Please make sure that every row has a value.")

        is_all_strings = data.iloc[:, col].apply(lambda x: isinstance(x, str) or pd.isna(x)).all()

        if not is_all_strings:
            raise ValueError(f"Column {col + 1} should contain only strings")

    # Control if the columns with the sgRNA names and sequences from the table are both unique
    log_.debug("Checking the uniqueness of the sgRNA names and sequences")
    if len(data) != len(set(data.iloc[:, 0])) or len(data) != len(set(data.iloc[:, 1])):
        raise ValueError(
            "Both columns with sgRNA names and sequences need to be unique. If this is true, please check whether "
            "these two columns are the first two columns of the input data file (as expected).")

    # Check if there is a nohit_row
    log_.debug("Identifying the nohit_row")
    nohit_row = data[data.iloc[:, 0] == 'nohit_row']
    if len(nohit_row) == 0:
        raise ValueError("There is no nohit_row. Please check your dataset for a row that is called 'nohit_row'.")


def clean_up(data, essential_genes, non_essential_genes, unwanted_columns=None, unwanted_rows=None,
             unwanted_row_substrings=None):
    """
    This function removes unwanted rows and columns and makes the dataset ready for
    sequential processing.

    :param data: Dataset containing read counts for each gRNA.
    :param essential_genes: list of genes that are considered essential
    :param non_essential_genes: list of genes that are considered non-essential
    :param unwanted_columns: List of substrings used to match and delete columns. (Optional)
    :param unwanted_rows: List of substrings used to match and delete rows. (Optional)
    :param unwanted_row_substrings: List of substrings used to match and remove from sgRNA name. (Optional)

    :return: cleaned_data: Dataset with the created changes
    """

    # Create a copy of the current dataframe
    cleaned_data = data

    # Rename the first column to "sgRNA" to avoid inconsistent column names
    log_.debug("Renaming first column to 'sgRNA'")
    cleaned_data.rename(columns={cleaned_data.columns[0]: "sgRNA"}, inplace=True)

    if unwanted_columns:
        # Remove unwanted columns
        log_.debug("Removing unwanted columns")
        log_.debug(f"Columns to be removed: {', '.join(unwanted_columns)}")
        # For this dataset: Remaining columns either contain the information for each condition or belong to "guide_mm1"
        length_before = len(cleaned_data.columns)
        cleaned_data = cleaned_data.drop(
            columns=cleaned_data.filter(regex="|".join(unwanted_columns)))
        if length_before != len(cleaned_data.columns):
            log_.debug(
                f"Number of columns before: {length_before}, after: {len(cleaned_data.columns)}, diff: {length_before - len(cleaned_data.columns)}")
        else:
            log_.debug("No unwanted columns identified")

    if unwanted_rows:
        # Remove rows with multi-matches (contain "multi:mismatch1" in column "name")
        log_.debug("Removing unwanted rows")
        log_.debug(
            f"Rows containing (one of) the following in the sgRNA name will be removed: {', '.join(unwanted_rows)}")
        rows_to_be_removed = cleaned_data["sgRNA"].str.contains('|'.join(unwanted_rows))
        if any(rows_to_be_removed):
            log_.debug(
                f"Number of rows before: {len(cleaned_data)}, after: {len(cleaned_data[~rows_to_be_removed])}, diff: {len(cleaned_data[rows_to_be_removed])}")
            cleaned_data = cleaned_data[~rows_to_be_removed]
        else:
            log_.debug("No unwanted rows identified")

    if unwanted_row_substrings:
        # Remove unwanted substrings from the gRNA column
        log_.debug("Removing unwanted substrings from rows of column sgRNA")
        log_.debug(
            f"The following substrings will be removed from the sgRNA name: {', '.join(unwanted_row_substrings)}")
        rows_to_be_edited = cleaned_data[cleaned_data["sgRNA"].str.contains('|'.join(unwanted_row_substrings))]
        if len(rows_to_be_edited) > 0:
            for substring in unwanted_row_substrings:
                cleaned_data['sgRNA'] = cleaned_data['sgRNA'].str.replace(substring + '.*', '', regex=True)
            log_.debug(f"Number of edited rows: {len(rows_to_be_edited)}")
        else:
            log_.debug("No unwanted row substrings identified")

    # The column names are converted to a common naming convention (stated in the documentation)
    # to avoid inconsistency in the output files between different screens.
    # e.g. T1_3D_to_2D_Rep3:TGGTCA is changed to t1_3d2d_r3
    log_.debug("Modifying column names")
    log_.debug("Changing names like 'T1_3D_to_2D_Rep3:TGGTCA' to 't1_3d2d_r3'")
    data_columns = cleaned_data.select_dtypes(include='number')
    annotations = cleaned_data.select_dtypes(exclude='number')
    data_columns.columns = (data_columns.columns.str.replace(':.*', '', regex=True)
                            .str.lower()
                            .str.replace("untreated", "ut")
                            .str.replace("treated", "tr")
                            .str.replace("_rep", "_r")
                            .str.replace("3d_to_2d", "3d2d")
                            .str.replace(r"\(|\)", ""))

    # In case the column name structures are not equal across all columns:
    # e.g. change T2_Treated_(3D)_Rep2:ATTGGC to t2_3d_tr_r2
    new_columns = []
    for column in data_columns.columns:
        if "(" in column:
            new_columns += [re.sub(r'_(tr|ut)_\((.*?)\)_(r\d)', r'_\2_\1_\3', column)]
        else:
            new_columns += [column]
    data_columns.columns = new_columns
    cleaned_data = pd.concat([annotations, data_columns], axis=1)

    # Create a new column called 'Gene' that contains the names of the genes
    # this is done by removing the extensions from the sgRNA names
    log_.debug("Creating a column for the gene names")
    pattern = r"-[0-9]*$"
    cleaned_data.insert(1, "Gene", [re.sub(pattern, "", gene) for gene in cleaned_data["sgRNA"]])

    # Assign a type to each gRNA
    log_.debug("Creating a column with the assigned data type for each gene")

    cleaned_data.insert(2, "type", assign_type(cleaned_data['Gene'], essential_genes, non_essential_genes))

    log_.debug(
        f"Dataset composition:\nNumber of essential genes (p) = {len(set(cleaned_data['Gene'][cleaned_data['type'] == 'p']))}, "
        f"Number of sgRNAs = {len(cleaned_data['Gene'][cleaned_data['type'] == 'p'])}\n"
        f"Number of non-essential genes (n) = {len(set(cleaned_data['Gene'][cleaned_data['type'] == 'n']))}, "
        f"Number of sgRNAs = {len(cleaned_data['Gene'][cleaned_data['type'] == 'n'])}\n"
        f"Number of non-targeting controls (o) = {len(cleaned_data['Gene'][cleaned_data['type'] == 'o'])}\n"
        f"Number of undefined genes (x) = {len(set(cleaned_data['Gene'][cleaned_data['type'] == 'x']))}, "
        f"Number of sgRNAs = {len(cleaned_data['Gene'][cleaned_data['type'] == 'x'])}\n")

    # Return the cleaned dataframe
    return cleaned_data


def create_summary(data, conditions, nohit_guide, summary_file):
    """
    This function creates a summary file that show the sums of reads per condition.
    The program makes three different classes (condition, condition with 1 mismatch and
    no match to any condition) and calculates the sum of reads for each of these classes per condition.

    :param data: Cleaned data file with read counts per condition
    :param conditions: A list of all conditions from the screen
    :param nohit_guide: The row with the reads per condition that could not be mapped to any guide (nohit_row)
    :param summary_file: Name of the output file
    :return:
    """

    log_.debug("Creating a data summary file")

    # Calculate the sum for all columns
    log_.debug("Calculating sum of each column")
    data_sum = data.sum(numeric_only=True, axis=0)

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


def compare_with_library(data, library, essential_genes, non_essential_genes):
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
    for gene_type in ["p", "n", "o", "x"]:
        if len(lib_copy[lib_copy['type'] == gene_type]) != 0:
            ratio = round(100 * len(data[data['type'] == gene_type]) / len(lib_copy[lib_copy['type'] == gene_type]), 1)
        else:
            ratio = "-Inf"
        log_.info(
            f"Number of type {gene_type} sgRNAs in library: {len(lib_copy[lib_copy['type'] == gene_type])}, in dataset: {len(data[data['type'] == gene_type])}, in percentage: {ratio} %")


def data_preparation(args):
    """
    Perform pre-processing to establish whether the quality of the screen output is
    good enough to continue with the hit identification.

    :param args: All input parameters stated in the run_analysis.py file

    :return:
    """
    log_.debug("Initiating data preparation")

    # Read the provided input files
    reads = pd.read_csv(args.input_file, sep="\t", header=0)
    essential_genes = pd.read_csv(args.essential_genes)
    non_essential_genes = pd.read_csv(args.non_essential_genes)
    library = pd.read_csv(args.library_file, sep="\t")

    # Check the input file for its correctness
    check_input(reads)

    # Clean up the input file, so it only contains required columns and rows
    reads = clean_up(reads, essential_genes, non_essential_genes, args.unwanted_columns, args.unwanted_rows,
                     args.unwanted_row_substrings)

    # Store all column names representing the conditions and all columns that have "guide_mm1" in their name
    log_.debug("Selecting columns that are used for creating the counts_summary file")
    data_columns = reads.select_dtypes(include='number').columns
    guide_mm1 = []
    conditions = [col for col in data_columns if not re.search("guide_mm1", col) or guide_mm1.append(col)]

    # Create a summary file if wanted by the user
    log_.debug("Selecting the nohit_row and removing it from the table")
    nohit_guide = reads[reads.iloc[:, 0] == 'nohit_row'][conditions]
    reads = reads.drop(nohit_guide.index)
    create_summary(data=reads, conditions=conditions, nohit_guide=nohit_guide, summary_file=args.summary_file)

    # Remove columns containing "guide_mm1"
    log_.debug(f"Removing guide_mm1 columns")
    reads = reads.drop(columns=reads.filter(regex="guide_mm1"))

    # Remove the rows that have a total sum of counts below a certain threshold
    if args.threshold == 0:
        log_.debug(f"Removing sgRNAs with a sum over all columns of {args.threshold}")
    else:
        log_.debug(f"Removing sgRNAs with a sum over all columns below {args.threshold}")
    # Calculate row sums for the data columns
    row_sums = reads.sum(numeric_only=True, axis=1)
    # Check if the sum of the total and filtered rows are the same
    if len(reads) != len(reads[row_sums > args.threshold]):
        log_.debug(
            f"Number of rows before filtering: {len(reads)}, after filtering: {len(reads[row_sums > args.threshold])}, "
            f"difference: {len(reads[row_sums <= args.threshold])}")
        reads = reads[row_sums > args.threshold]
    else:
        log_.debug(f"No sgRNAs count sums below {args.threshold} identified")

    # Compare the cleaned data file with the given library file
    compare_with_library(reads, library, essential_genes, non_essential_genes)

    # Store the cleaned data in a csv file
    log_.debug(f"Storing the resulting dataframe in {args.output_file}")
    reads.to_csv(args.output_file, sep=';', index=False)

    log_.info("Normalizing data and creating quality control plots\n")
    # Remove the replicate notation ("_r1") from the condition names
    conditions_new = [condition.rsplit('_', 1)[0] for condition in conditions]
    run_script(r"D:\D\Ausbildung\Master\1st year\Internships\NKI\Report\Program test\Program files\R_analysis_1.R",
               additional_args=[args.output_file, args.dataset, ",".join(conditions_new), "biological",
                                args.target_samples, args.reference_samples, args.distribution_condition1,
                                args.distribution_condition2])
    log_.debug("Finished Data preparation")
