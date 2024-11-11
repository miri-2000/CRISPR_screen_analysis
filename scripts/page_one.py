import tkinter.ttk as ttk
from base_frame import BaseFrame
from tkinter import messagebox
# from input_validation_gui_stateless import InputValidatorGUI

class StartPage(BaseFrame):
    """Start page of the CRISPR screen analysis tool."""

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        self.controller = controller

        # Define the file paths and associated labels
        self.file_paths = {
            "Screen Result File": controller.shared_data["input_file"],
            "Essential Genes File": controller.shared_data["essential_genes"],
            "Non-Essential Genes File": controller.shared_data["non_essential_genes"],
            "Library File": controller.shared_data["library_file"]
        }

        self.create_header("Welcome to the CRISPR screen analysis tool")
        self.create_description("Welcome to the CRISPR screen analysis tool for the identification\n"
                                "of drug-gene interactions using DrugZ! This tool is designed to analyze\n"
                                "CRISPR screens and identify potential drug-gene interactions. To get started,\n"
                                "please provide the following parameters:")
        self.create_widgets()

        # Button to navigate to the next page
        ttk.Button(self, text="Next", command=self.validate_and_proceed).grid(row=9, column=1, pady=20)

    def create_widgets(self):
        """Create input fields and buttons for user parameters."""
        # Targets entry
        self.create_labeled_entry(label_text="Name of treated conditions/timepoints:", command_action=lambda: self.show_info(
            "Treated conditions (e.g. with exposure to drugs,...) to be compared against baseline conditions"),
                                  textvariable=self.controller.shared_data["target_samples"], row=2)

        # Add a label or star indicator next to each entry
        self.add_indicator_label(row=2, label_text="Name of treated conditions/timepoints:")

        # References entry
        self.create_labeled_entry("Name of untreated/baseline conditions/timepoints:",
                                  lambda: self.show_info("Baseline conditions (e.g. with no added drugs)"),
                                  self.controller.shared_data["reference_samples"], 3)

        # Add a label or star indicator next to each entry
        self.add_indicator_label(row=3, label_text="Name of untreated/baseline conditions/timepoints:")

        for i, (label_text, var) in enumerate(self.file_paths.items(), start=4):
            self.create_labeled_entry(f"{label_text}:",
                                      lambda l=label_text: self.show_info(f"File path to the {l}"),
                                      var, i)

            # Create a browse button
            self.create_browse_button(var, row=i)

            # Add a label or star indicator next to each entry
            self.add_indicator_label(row=i, label_text=label_text)

    def validate_and_proceed(self):
        # validator = InputValidatorGUI()
        # if validator.validate_page_one(self.file_paths):
        #     self.controller.show_frame(self.get_next_class())
        """Validate input fields and navigate to the next page if valid."""
        self.invalid_file_types.clear()
        self.invalid_values.clear()

        self.add_trace(self.controller.shared_data["target_samples"], self.indicator_labels["Name of target samples:"],
                       "non_empty")
        self.add_trace(self.controller.shared_data["reference_samples"],
                       self.indicator_labels["Name of reference samples:"], "non_empty")

        invalid_value_flag = self.check_invalid_values()
        invalid_file_type_flag = self.check_invalid_file_types()

        if not invalid_value_flag and not invalid_file_type_flag:
            self.controller.show_frame(self.get_next_class())

    def check_invalid_file_types(self):
        """Check for any invalid file types."""
        for label_text, var in self.file_paths.items():
            self.add_trace(var, self.indicator_labels[label_text], "file_type")

        invalid_value_labels = [k for k, v in self.invalid_file_types.items() if v is True]
        if invalid_value_labels:
            messagebox.showwarning("Invalid File Paths/Types",
                                   f"The following paths are invalid: {', '.join(invalid_value_labels)}. The "
                                   f"file path must be accurate and lead to a text or csv file.")
            return True

    @staticmethod
    def get_next_class():
        """Get the next class for navigation."""
        from page_two import PageTwo

        # Return the next class
        return PageTwo
