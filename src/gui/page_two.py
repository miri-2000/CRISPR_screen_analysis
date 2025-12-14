import tkinter.ttk as ttk
from src.gui.base_frame import BaseFrame
from src.core.input_validation_gui import InputValidatorGUI


class PageTwo(BaseFrame):
    """Page two of the CRISPR screen core tool."""

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Define class attributes
        self.threshold_reads_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_reads"])

        self.create_header("Data Preparation Settings")
        self.create_description("Optional data preparation settings (default values shown in line)")
        self.create_widgets()

        # Button to navigate to the previous page
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(self.get_previous_class())).grid(row=9,
                                                                                                             column=1,
                                                                                                             pady=10)
        # Button to navigate to the next page
        ttk.Button(self, text="Next", command=self.validate_and_proceed).grid(row=9,
                                                                              column=2,
                                                                              pady=10)

    def create_widgets(self):
        """Create input fields and buttons for user parameters."""
        # Unwanted columns entry
        self.create_labeled_entry("Unwanted Columns:", lambda: self.show_info("Columns that should be removed "
                                                                              "from the input file"),
                                  self.controller.shared_data["unwanted_columns"],
                                  2)

        # Unwanted rows entry
        self.create_labeled_entry("Unwanted Rows:", lambda: self.show_info("Rows that should be removed from "
                                                                           "the input file"),
                                  self.controller.shared_data["unwanted_rows"], 3)

        # Unwanted row substrings entry
        self.create_labeled_entry("Unwanted Row Substrings:", lambda: self.show_info("Row substrings that "
                                                                                     "should be removed from gRNA "
                                                                                     "names"),
                                  self.controller.shared_data["unwanted_row_substrings"], 4)

        # Min number of reads entry
        self.create_labeled_entry("Minimum required sum of reads/guide:",
                                  lambda: self.show_info("Minimum number of total reads/guide so that "
                                                         "the guide will not be discarded"),
                                  self.threshold_reads_entry, 5, "0")

    def validate_and_proceed(self):
        validator = InputValidatorGUI()
        if validator.validate_page_two(self.threshold_reads_entry.get()):
            self.controller.show_frame(self.get_next_class())

    @staticmethod
    def get_previous_class():
        from page_one import StartPage

        # Return the previous class
        return StartPage

    @staticmethod
    def get_next_class():
        from src.gui.page_three import PageThree

        # Return the next class
        return PageThree
