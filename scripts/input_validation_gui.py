from tkinter import messagebox
from input_validation import InputValidator


class InputValidatorGUI(InputValidator):
    def __init__(self):
        super().__init__()

    def validate_page_one(self, file_paths):
        input_file = file_paths['Screen Result File'].get()
        essential_genes = file_paths['Essential Genes File'].get()
        non_essential_genes = file_paths['Non-Essential Genes File'].get()
        library_file = file_paths['Library File'].get()

        if not self.validate_files(input_file, essential_genes, non_essential_genes, library_file):
            return False

        return True

    def validate_page_two(self, threshold_reads):
        if not self.validate_int_fields(threshold_reads=threshold_reads):
            return False

        return True

    def validate_page_three(self, top, threshold_fdr, working_dir):
        if not self.validate_int_fields(top=top):
            return False

        if not self.validate_threshold_fdr(threshold_fdr):
            return False

        if not self.validate_working_dir(working_dir):
            return False

        return True

    def abort(self, message):
        messagebox.showinfo('Input Error', message)
