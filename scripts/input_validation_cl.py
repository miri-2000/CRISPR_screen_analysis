import sys
from input_validation import InputValidator


class InputValidatorCL(InputValidator):
    def __init__(self):
        super().__init__()

    def validate(self, input_file, essential_genes, non_essential_genes, library_file, threashold_reads, top,
                 threshold_fdr, x_axis, replicate_type, working_dir):

        self.validate_files(input_file, essential_genes, non_essential_genes, library_file)
        # self.validate_samples()
        self.validate_int_fields(threshold_reads=threashold_reads, top=top)
        self.validate_threshold_fdr(threshold_fdr)
        self.validate_choice_fields(x_axis, replicate_type)
        self.validate_working_dir(working_dir)

    def abort(self, message):
        sys.exit(message)
