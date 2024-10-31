from pathlib import Path
import pytest
from file_based_launch import Args
from input_validation_cl import InputValidatorCL


@pytest.fixture
def args_instance():
    """Fixture to create a valid Args instance for each tests."""
    return Args()


class TestInputValidator:

    def validate_and_assert(self, expected_message, args_instance):
        """Helper method to validate input_file and assert the error message."""
        validator = InputValidatorCL(args_instance)

        with pytest.raises(SystemExit) as excinfo:
            validator.validate()  # Assuming validate_file_path checks this

        # Check the output message
        assert str(excinfo.value) == expected_message

    def test_missing_variable(self, args_instance):
        """Test with an empty input_file."""
        args_instance.input_file = None

        expected_message = ("The parameter 'input_file' needs to be specified. If the parameter should remain empty set"
                            f"it by saying e.g. 'unwanted_rows = ""'")
        self.validate_and_assert(expected_message, args_instance)

    def test_empty_input_file(self, args_instance):
        """Test with an empty input_file."""
        args_instance.input_file = ""

        expected_message = ("The parameter 'input_file' needs to contain a valid file path that leads to the "
                            "specified file. The file must either be a CSV or TXT file.")
        self.validate_and_assert(expected_message, args_instance)

    def test_input_file_with_invalid_path(self, args_instance):
        """Test with an invalid input_file path."""
        args_instance.input_file = r"/path/that/does/not/exist.csv"

        expected_message = ("The parameter 'input_file' needs to contain a valid file path that leads to the "
                            "specified file. The file must either be a CSV or TXT file.")
        self.validate_and_assert(expected_message, args_instance)

    def test_input_file_with_wrong_file_type(self, args_instance):
        """Test with an input file with the wrong file type."""
        args_instance.input_file = Path(__file__).parents[1] / "tests" / "test_file.json"

        expected_message = ("The parameter 'input_file' needs to contain a valid file path that leads to the "
                            "specified file. The file must either be a CSV or TXT file.")
        self.validate_and_assert(expected_message, args_instance)

    def test_input_file_with_non_unique_first_columns(self, args_instance):
        """Test input_file with non-unique first columns."""
        args_instance.input_file = Path(__file__).parents[1] / "tests" / "read_count_file_with_duplicated_rows.txt"

        expected_message = (
            "The CRISPR screen input file requires the first column with the sgRNA names and the second column "
            f"with the sgRNA sequences to be unique.")

        self.validate_and_assert(expected_message, args_instance)

    def test_input_file_without_nohit_row(self, args_instance):
        """Test input file without a no-hit row."""
        args_instance.input_file = Path(__file__).parents[1] / "tests" / "read_count_file_without_nohit_row.txt"

        expected_message = (
            "The CRISPR screen input file requires a no-hit row.")

        self.validate_and_assert(expected_message, args_instance)

    def test_input_file_without_guide_mm1_column(self, args_instance):
        """Test input_file without a guide mm1 column."""
        args_instance.input_file = Path(__file__).parents[1] / "tests" / "read_count_file_without_guide_mm1_columns.txt"

        expected_message = (
            "The CRISPR screen input file requires a guide_mm1_ column.")

        self.validate_and_assert(expected_message, args_instance)

    def test_essential_genes_not_in_input_file(self, args_instance):
        """Test essential_genes file with genes that are not in input file."""
        args_instance.essential_genes = Path(__file__).parents[
                                            1] / "tests" / "essential_genes_with_genes_not_in_input_file.csv"

        expected_message = (
            "All genes mentioned in the 'essential_genes' file need to be present in the CRISPR screen "
            f"input file.")

        self.validate_and_assert(expected_message, args_instance)

    def test_library_file_with_renamed_required_columns(self, args_instance):
        """Test library file with renamed required columns."""
        args_instance.library_file = Path(__file__).parents[
                                         1] / "tests" / "library_file_with_renamed_required_columns.txt"

        expected_message = (
            "The library file requires the column 'Target Gene Symbol' (holding the gene names) and"
            "'sgRNA Target Sequence' (holding the sgRNA sequence) to be present in the CRISPR screen.")

        self.validate_and_assert(expected_message, args_instance)

    def test_library_file_with_missing_required_columns(self, args_instance):
        """Test library_file with missing required columns."""
        args_instance.library_file = Path(__file__).parents[
                                         1] / "tests" / "library_file_with_missing_required_columns.txt"

        expected_message = (
            "The library file requires the column 'Target Gene Symbol' (holding the gene names) and"
            "'sgRNA Target Sequence' (holding the sgRNA sequence) to be present in the CRISPR screen.")

        self.validate_and_assert(expected_message, args_instance)

    def test_int_field_empty(self, args_instance):
        """Test int_files with an empty string."""
        args_instance.threshold_reads = ""

        expected_message = (
            "The parameter 'threshold_reads' needs to be a whole number.")

        self.validate_and_assert(expected_message, args_instance)

    def test_int_field_wrong_data_type(self, args_instance):
        """Test int_file with a wrong data type."""
        args_instance.threshold_reads = 10.5

        expected_message = (
            "The parameter 'threshold_reads' needs to be a whole number.")

        self.validate_and_assert(expected_message, args_instance)

    def test_threshold_reads_empty(self, args_instance):
        """Test threshold_reads with an empty string."""
        args_instance.threshold_fdr = ""

        expected_message = (
            "The parameter 'threshold_fdr' needs to be a number between 0 and 1.")

        self.validate_and_assert(expected_message, args_instance)

    def test_threshold_reads_number_out_of_limit(self, args_instance):
        """Test threshold_reads with a number that is not between 0 and 1."""
        args_instance.threshold_fdr = 10

        expected_message = (
            "The parameter 'threshold_fdr' needs to be a number between 0 and 1.")

        self.validate_and_assert(expected_message, args_instance)

    def test_choice_fields_invalid_choice(self, args_instance):
        """Test a choice_field with an invalid choice."""
        args_instance.x_axis = 10

        expected_message = (
            "The x_axis can only either be 'normZ' or 'log2 fold-change'.")

        self.validate_and_assert(expected_message, args_instance)

    def test_working_dir_wrong_data_type(self, args_instance):
        """Test working_dir with the wrong data type."""
        args_instance.working_dir = 10

        expected_message = (
            "The given working directory needs to be a directory. Please make sure you inserted a folder "
            "location.")

        self.validate_and_assert(expected_message, args_instance)

    def test_working_dir_invalid_path(self, args_instance):
        """Test working_dir with an invalid path."""
        args_instance.working_dir = Path(__file__).parents[
                                        1] / "tests" / "library_file_with_missing_required_columns.txt"

        expected_message = (
            "The given working directory needs to be a directory. Please make sure you inserted a folder "
            "location.")

        self.validate_and_assert(expected_message, args_instance)
