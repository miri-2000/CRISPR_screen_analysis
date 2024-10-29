import sys
from constants import INPUT_FIELDS
from input_validation import InputValidator


class InputValidatorCL(InputValidator):
    def __init__(self, args):
        super().__init__(args)

    def validate(self):
        input_fields = INPUT_FIELDS
        for field in input_fields:
            value = getattr(self, field, None)
            if value is None:
                sys.exit(
                    f"The parameter '{field}' needs to be specified. If the parameter should remain empty set"
                    f"it by saying e.g. 'unwanted_rows = ""'")

        self.validate_files()
        # self.validate_samples()
        self.validate_int_fields()
        self.validate_threshold_fdr()
        self.validate_choice_fields()
        self.validate_working_dir()

    def abort(self, message):
        sys.exit(message)
