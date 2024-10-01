import tkinter.ttk as ttk
from base_frame import BaseFrame
from tkinter import messagebox
from start_program import CRISPR_screen_analysis


class PageThree(BaseFrame):
    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Visualization and Storage Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=3)

        # Text label
        text_label = ttk.Label(self,
                               text="Visualization settings (default values shown in line) and storage settings\n",
                               justify="center")
        text_label.grid(row=1, column=0, columnspan=3, pady=(5, 0))

        # Define class attributes
        self.top_entry = ttk.Entry(self, textvariable=self.controller.shared_data["top"])
        self.threshold_fdr_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_fdr"])
        self.condition = {"Positive Control": self.controller.shared_data["distribution_condition1"],
                          "Negative Control": self.controller.shared_data["distribution_condition2"]}

        self.create_widgets()

        # Button to navigate to the previous page
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(self.get_previous_class())).grid(row=11, column=1, pady=10)

        # Button to navigate to start the program
        ttk.Button(self, text="Start Computation", command=self.validate_and_proceed).grid(row=11, column=2, pady=10)

    def create_widgets(self):
        # Title label
        text_label = ttk.Label(self,
                               text="Gene significance plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=2, column=0, sticky="w", pady=(3, 5))

        # Set the initial value of the selected option
        self.controller.shared_data["x_axis"].set("normZ")  # Make Option 1 the default

        # Create radiobuttons with different text and values
        ttk.Label(self, text="X-axis value:").grid(row=3, column=0, sticky="w")
        ttk.Button(self, text="?", command=lambda: self.show_info("Metric that should be taken for the x axis of the "
                                                                  "distribution plot"), width=1).grid(row=3, column=0,
                                                                                                      sticky="e")
        ttk.Radiobutton(self, text="normZ", variable=self.controller.shared_data["x_axis"], value="normZ").grid(
            row=3, column=1, sticky="w")
        ttk.Radiobutton(self, text="log2 fold-change", variable=self.controller.shared_data["x_axis"],
                        value="log2 fold-change").grid(row=3,
                                                       column=2,
                                                       sticky="w")

        # Gene significance threshold
        self.create_labeled_entry("Significance threshold:",
                                                 lambda: self.show_info("Genes with a FDR above the defined threshold "
                                                                        "will not be considered significant"),
                                                 self.threshold_fdr_entry, 4, "0.25")

        # Top hits entry
        self.create_labeled_entry("Number of hits per plot:",
                                                 lambda: self.show_info("Number defines the maximum number of genes "
                                                                        "that will be displayed"),
                                                 self.top_entry, 5, "15")

        # Title label
        text_label = ttk.Label(self,
                               text="Distribution Plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=6, column=0, sticky="w", pady=(10, 5))

        # Positive Control entry
        self.create_labeled_entry("Positive Control:", lambda: self.show_info("Sample that should be taken as the "
                                                                              "positive control of the screen"),
                                  self.controller.shared_data["distribution_condition1"], 7)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=7, column=2, sticky="w")
        self.indicator_labels["Positive Control"] = indicator_label

        # Negative Control entry
        self.create_labeled_entry("Negative Control:", lambda: self.show_info("Sample that should be taken as the "
                                                                              "negative control of the screen"),
                                  self.controller.shared_data["distribution_condition2"], 8)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=8, column=2, sticky="w")
        self.indicator_labels["Negative Control"] = indicator_label

        # Title label
        text_label = ttk.Label(self,
                               text="Storage settings",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=9, column=0, sticky="w", pady=(10, 5))

        self.create_labeled_entry("Result storage location:",
                                  lambda: self.show_info("Directory where the results should be stored"),
                                  self.controller.shared_data["working_dir"], 10)

        ttk.Button(self, text="Browse",
                   command=lambda: self.browse_files(self.controller.shared_data["working_dir"])).grid(row=10, column=2,
                                                                                                       sticky="e")

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=10, column=3, sticky="w")
        self.indicator_labels["Result storage location"] = indicator_label

    def validate_and_proceed(self):
        self.invalid_values.clear()

        for label_text, var in self.indicator_labels.items():
            if label_text != "Result storage location":
                self.add_trace(self.condition[label_text], var, "non_empty")
            else:
                self.add_trace(self.controller.shared_data["working_dir"], var, "folder_location")

        invalid_value_labels = [k for k, v in self.invalid_values.items() if v is True]
        invalid_folder_location_labels = [k for k, v in self.invalid_folder_locations.items() if v is True]
        if invalid_value_labels:
            messagebox.showwarning("Missing information",
                                   f"The following fields cannot remain empty: {', '.join(invalid_value_labels)}. "
                                   f"Please specify the parameters.")
        elif invalid_folder_location_labels:
            messagebox.showwarning("Not a directory",
                                   f"The following field/s is/are not (a) directory/ies: {', '.join(invalid_value_labels)}. "
                                   f"Please specify a valid folder location to store the results.")
        else:
            self.start_computation()

    def save_shared_data(self):
        # Iterate through shared_data and store values as strings
        for key, variable in self.controller.shared_data.items():
            # Get the value of the variable and convert it to string/int
            try:
                value = int(variable.get())
            except ValueError:
                value = variable.get()
            # value = str(variable.get())

            # Dynamically create instance variables using setattr
            setattr(self, key, value)

            # Print for demonstration (you can save it to a file or use it as needed)
            print(f"{key} = {getattr(self, key)}")

    def start_computation(self):

        print(self.controller.shared_data)
        # Perform conversion
        self.save_shared_data()

        # Call the function from another file with self as the argument
        CRISPR_screen_analysis(self)

    @staticmethod
    def get_previous_class():
        from page_two import PageTwo
        # Return the previous class
        return PageTwo
