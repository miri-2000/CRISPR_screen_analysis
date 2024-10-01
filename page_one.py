import tkinter as tk
import tkinter.ttk as ttk
from base_frame import BaseFrame
from tkinter import filedialog, messagebox
import os


class StartPage(BaseFrame):
    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Welcome to the CRISPR screen analysis tool",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=3, padx=20)

        # Text label
        text_label = ttk.Label(self,
                               text="Welcome to the CRISPR screen analysis tool for the identification\n"
                                    "of drug-gene interactions using DrugZ! This tool is designed to analyze\n"
                                    "CRISPR screens and identify potential drug-gene interactions. To get started,\n"
                                    "please provide the following parameters:", justify="center")
        text_label.grid(row=1, column=0, columnspan=3, pady=(5, 15), padx=20)

        self.controller = controller

        # Define the file paths and associated labels
        self.file_paths = {
            "Screen Result File": controller.shared_data["input_file"],
            "Essential Genes File": controller.shared_data["essential_genes"],
            "Non-Essential Genes File": controller.shared_data["non_essential_genes"],
            "Library File": controller.shared_data["library_file"]
        }

        self.tracers_added = False  # Flag to check if tracers are added

        self.create_widgets()

        # Button to navigate to the next page
        ttk.Button(self, text="Next Page", command=self.validate_and_proceed).grid(row=9, column=1, pady=20)

    def create_widgets(self):
        # Targets entry
        self.create_labeled_entry(label_text="Name of target samples:", command_action=lambda: self.show_info(
            "Treated conditions (e.g. with exposure to drugs,...) to be compared against baseline conditions"),
                                  textvariable=self.controller.shared_data["target_samples"], row=2)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=2, column=3, sticky="w")
        self.indicator_labels["Name of target samples:"] = indicator_label

        # References entry
        self.create_labeled_entry("Name of reference samples:",
                                  lambda: self.show_info("Baseline conditions (e.g. with no added drugs)"),
                                  self.controller.shared_data["reference_samples"], 3)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=3, column=3, sticky="w")
        self.indicator_labels["Name of reference samples:"] = indicator_label

        for i, (label_text, var) in enumerate(self.file_paths.items(), start=4):
            self.create_labeled_entry(f"{label_text}:",
                                      lambda l=label_text: self.show_info(f"File path to the {l}"),
                                      var, i)

            ttk.Button(self, text="Browse",
                       command=lambda var=self.file_paths[label_text]: self.browse_files(var)).grid(row=i, column=2,
                                                                                                    sticky="w")

            # Add a label or star indicator next to each entry
            indicator_label = ttk.Label(self, text="")
            indicator_label.grid(row=i, column=3, sticky="w")
            self.indicator_labels[label_text] = indicator_label

    def validate_and_proceed(self):
        self.invalid_file_types.clear()
        self.invalid_values.clear()

        self.add_trace(self.controller.shared_data["target_samples"], self.indicator_labels["Name of target samples:"],
                       "non_empty")
        self.add_trace(self.controller.shared_data["reference_samples"],
                       self.indicator_labels["Name of reference samples:"], "non_empty")

        invalid_value_labels = [k for k, v in self.invalid_values.items() if v is True]
        if invalid_value_labels:
            messagebox.showwarning("Missing information",
                                   f"The following fields cannot remain empty: {', '.join(invalid_value_labels)}. "
                                   f"Please specify the parameters.")

        for label_text, var in self.file_paths.items():
            self.add_trace(var, self.indicator_labels[label_text], "file_type")

        invalid_value_labels = [k for k, v in self.invalid_file_types.items() if v is True]
        if invalid_value_labels:
            messagebox.showwarning("Invalid Paths",
                                   f"The following paths are invalid: {', '.join(invalid_value_labels)}. The "
                                   f"file path must be accurate and lead to a text or csv file.")

        if not invalid_value_labels and not invalid_value_labels:
            self.controller.show_frame(self.get_next_class())

    @staticmethod
    def get_next_class():
        from page_two import PageTwo
        # Return the next class
        return PageTwo
