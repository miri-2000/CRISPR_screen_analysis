import tkinter.ttk as ttk
from base_frame import BaseFrame
from tkinter import messagebox
from start_program import CRISPRScreenAnalysis


class PageThree(BaseFrame):
    """Page three of the CRISPR screen analysis tool."""

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Define class attributes
        self.top_entry = ttk.Entry(self, textvariable=self.controller.shared_data["top"])
        self.threshold_fdr_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_fdr"])
        self.condition = {"Positive Control": self.controller.shared_data["distribution_condition1"],
                          "Negative Control": self.controller.shared_data["distribution_condition2"]}

        self.create_header("Visualization and Storage Settings")
        self.create_description("Visualization settings (default values shown in line) and storage settings\n")
        self.create_widgets()

        # Button to navigate to the previous page
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(self.get_previous_class())).grid(row=13,
                                                                                                             column=1,
                                                                                                             pady=10)
        # Button to navigate to start the program
        ttk.Button(self, text="Start Computation", command=self.validate_and_proceed).grid(row=13, column=2, pady=10)

    def create_widgets(self):
        """Create input fields and buttons for user parameters."""
        # Title label
        self.create_subheader("Gene significance plots", 2)

        # Set the initial value of the selected option
        self.controller.shared_data["x_axis"].set("normZ")

        # Create radiobuttons with different text and values
        self.create_radio_buttons("X-axis value:", 3,
                                  lambda: self.show_info("Metric that should be taken for the x axis of the "
                                                         "distribution plot"), "normZ", "log2 fold-change",
                                  self.controller.shared_data["x_axis"])

        # Gene significance threshold
        self.create_labeled_entry("Significance threshold:",
                                  lambda: self.show_info("Maximum FDR for a gene to be considered significant"),
                                  self.threshold_fdr_entry, 4, "0.25")

        # Top hits entry
        self.create_labeled_entry("Number of hits per plot:",
                                  lambda: self.show_info("Maximum number of significant genes "
                                                         "that will be displayed"),
                                  self.top_entry, 5, "15")

        # Title label
        self.create_subheader("Distribution Plot", 6)

        # Positive Control entry
        self.create_labeled_entry("Positive Control:", lambda: self.show_info("Sample that should be taken as the "
                                                                              "positive control of the screen"),
                                  self.controller.shared_data["distribution_condition1"], 7)

        # Add a label or star indicator next to each entry
        self.add_indicator_label(row=7, label_text="Positive Control")

        # Negative Control entry
        self.create_labeled_entry("Negative Control:", lambda: self.show_info("Sample that should be taken as the "
                                                                              "negative control of the screen"),
                                  self.controller.shared_data["distribution_condition2"], 8)

        # Add a label or star indicator next to each entry
        self.add_indicator_label(row=8, label_text="Negative Control")

        # Title label
        self.create_subheader("Correlation Plot", 9)

        # Set the initial value of the selected option
        self.controller.shared_data["replicate_type"].set("biological")

        # Create radiobuttons with different text and values
        self.create_radio_buttons("Replicate type:", 10,
                                  lambda: self.show_info(
                                      "Defines whether replicate samples are biological or technical replicates"),
                                  "biological", "technical",
                                  self.controller.shared_data["replicate_type"])

        # Title label
        self.create_subheader("Storage settings", 11)

        self.create_labeled_entry("Result storage location:",
                                  lambda: self.show_info("Directory where the results should be stored"),
                                  self.controller.shared_data["working_dir"], 12)

        self.create_browse_button(self.controller.shared_data["working_dir"], row=12)

        # Add a label or star indicator next to each entry
        self.add_indicator_label(12, label_text="Result storage location")

    def create_radio_buttons(self, label_text, row, command_action, button1_text, button2_text, variable):
        """Create radio buttons for user parameters."""
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=20)
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0,
                                                                                            sticky="e")
        ttk.Radiobutton(self, text=button1_text, variable=variable, value=button1_text).grid(
            row=row, column=1, sticky="w")
        ttk.Radiobutton(self, text=button2_text, variable=variable,
                        value=button2_text).grid(row=row,
                                                 column=2,
                                                 sticky="w")

    def validate_and_proceed(self):
        """Validate input fields and start program if valid."""
        self.invalid_values.clear()

        for label_text, var in self.indicator_labels.items():
            if label_text != "Result storage location":
                self.add_trace(self.condition[label_text], var, "non_empty")
            else:
                self.add_trace(self.controller.shared_data["working_dir"], var, "folder_location")

        invalid_value_flag = self.check_invalid_values()
        invalid_folder_location_flag = self.check_invalid_folder_location()

        if not invalid_value_flag and not invalid_folder_location_flag:
            self.start_computation()

    def check_invalid_folder_location(self):
        """Check for any invalid folder locations."""
        invalid_folder_location_labels = [k for k, v in self.invalid_folder_locations.items() if v is True]
        if invalid_folder_location_labels:
            messagebox.showwarning("Not a directory",
                                   f"The following field/s is/are not (a) directory/ies: {', '.join(
                                       invalid_folder_location_labels)}. Please specify a valid folder location "
                                   f"to store the results.")
            return True

    def start_computation(self):
        """Start the CRISPR screen analysis program."""

        # Prepare user input
        self.save_shared_data()

        # Call the CRISPR screen analysis program
        CRISPRScreenAnalysis(self)

    def save_shared_data(self):
        """Store the user input data for the start of the CRISPR screen analysis program."""
        # Iterate through shared_data and store values as strings
        for key, variable in self.controller.shared_data.items():
            # Get the value of the variable and convert it to string/int
            try:
                value = int(variable.get())
            except ValueError:
                value = variable.get()

            # Dynamically create instance variables using setattr
            setattr(self, key, value)

    @staticmethod
    def get_previous_class():
        """Get the next class for navigation."""
        from page_two import PageTwo

        # Return the previous class
        return PageTwo
