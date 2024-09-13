import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox
# from gui_base_frame import BaseFrame
# from gui_page_one import StartPage
# from gui_page_two import PageTwo
# from gui_page_three import PageThree
from run_analysis import CRISPR_screen_analysis
import os


class BaseFrame(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        self.invalid_values = {}  # Dictionary to hold invalid values
        self.invalid_file_types = {}  # Dictionary to hold invalid file types
        self.invalid_folder_locations = {}  # Dictionary to hold invalid folder locations
        self.indicator_labels = {}

    def create_labeled_entry(self, label_text, command_action, textvariable, row):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=(20, 0))
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0, sticky="e")
        ttk.Entry(self, textvariable=textvariable).grid(row=row, column=1, sticky="w")

    def create_labeled_entry_with_threshold(self, label_text, command_action, threshold_entry, default_value, row):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=20)
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0, sticky="e")
        threshold_entry.insert(0, default_value)
        threshold_entry.bind("<FocusIn>", lambda event: self.on_entry_click(event, threshold_entry, default_value))
        threshold_entry.bind("<FocusOut>", lambda event: self.on_entry_leave(event, threshold_entry, default_value))
        threshold_entry.grid(row=row, column=1)

    def show_info(self, label_text):
        messagebox.showinfo("Information", label_text)

    def on_entry_click(self, event, entry, default_text):
        if entry.get() == default_text:
            entry.delete(0, tk.END)

    def on_entry_leave(self, event, entry, default_text):
        if entry.get() == "":
            entry.insert(0, default_text)

    def add_trace(self, value, label, check):
        value.trace_add("write", self.update_indicator(value, label, check))

    def update_indicator(self, var, indicator_label, check):
        var_value = var.get()
        label_text = [k for k, v in self.indicator_labels.items() if v == indicator_label][0]

        if check == "non_empty":
            self.check_empty_file(indicator_label, label_text, var_value)
        elif check == "file_type":
            self.check_file_type(indicator_label, label_text, var_value)
        elif check == "folder_location":
            self.check_folder_location(indicator_label, label_text, var_value)

    def check_empty_file(self, indicator_label, label_text, var_value):
        if var_value.strip() == "":
            indicator_label.config(text=" * Invalid")
            self.invalid_values[label_text] = True
        else:
            indicator_label.config(text="")
            self.invalid_values[label_text] = False

    def check_file_type(self, indicator_label, label_text, path):
        print(f"Checking path: {path}")  # Debug print
        if not os.path.exists(path) or not (path.endswith('.txt') or path.endswith('.csv')):
            indicator_label.config(text=" * Invalid")
            print(f"Path is invalid: {path}")  # Debug print
            self.invalid_file_types[label_text] = True
        else:
            indicator_label.config(text="")
            print(f"Path is valid: {path}")  # Debug print
            self.invalid_file_types[label_text] = False

    def check_folder_location(self, indicator_label, label_text, path):
        print(f"Checking path: {path}")  # Debug print
        if not os.path.isdir(path):
            indicator_label.config(text=" * Invalid")
            print(f"Path is invalid: {path}")  # Debug print
            self.invalid_folder_locations[label_text] = True
        else:
            indicator_label.config(text="")
            print(f"Path is valid: {path}")  # Debug print
            self.invalid_folder_locations[label_text] = False

    def browse_files(self, output_variable):
        file_path = filedialog.askopenfilename()
        output_variable.set(file_path)


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
            self.controller.show_frame(PageTwo)


class PageTwo(BaseFrame):

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Data Preparation Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=3, padx=20)

        # Title label
        text_label = ttk.Label(self,
                               text="Optional data preparation settings (default values shown in line)")
        # ,justify="center")
        text_label.grid(row=1, column=0, columnspan=3, pady=(5, 15), padx=20)

        # Define class attributes
        self.threshold_reads_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_reads"])

        self.create_widgets()

        # Button to navigate to the previous page
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(StartPage)).grid(row=9, column=1, pady=10)

        # Button to navigate to the next page
        ttk.Button(self, text="Next", command=lambda: controller.show_frame(PageThree)).grid(row=9, column=2, pady=10)

    def create_widgets(self):
        # Unwanted columns entry
        self.create_labeled_entry("Unwanted Columns:", lambda: self.show_info("Name of the unwanted columns"),
                                  self.controller.shared_data["unwanted_columns"],
                                  2)

        # Unwanted rows entry
        self.create_labeled_entry("Unwanted Rows:", lambda: self.show_info("Name of the unwanted rows"),
                                  self.controller.shared_data["unwanted_rows"], 3)

        # Unwanted row substrings entry
        self.create_labeled_entry("Unwanted Row Substrings:", lambda: self.show_info("Name of the unwanted "
                                                                                     "row substrings"),
                                  self.controller.shared_data["unwanted_row_substrings"], 4)

        # Min number of reads entry
        self.create_labeled_entry_with_threshold("Minimum required sum of reads/guide:",
                                                 lambda: self.show_info("Minimum number of total reads/guide so that "
                                                                        "the guide will not be discarded"),
                                                 self.threshold_reads_entry, "0", 5)


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
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(PageTwo)).grid(row=11, column=1, pady=10)

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
        self.create_labeled_entry_with_threshold("Significance threshold:",
                                                 lambda: self.show_info("Genes with a FDR above the defined threshold "
                                                                        "will not be considered significant"),
                                                 self.threshold_fdr_entry, "0.25", 4)

        # Top hits entry
        self.create_labeled_entry_with_threshold("Number of hits per plot:",
                                                 lambda: self.show_info("Number defines the maximum number of genes "
                                                                        "that will be displayed"),
                                                 self.top_entry, "15", 5)

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


class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # Create a container to hold multiple frames
        container = tk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")

        # Allow the container to resize with the window
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        self.shared_data = {'input_file': tk.StringVar(),
                            'essential_genes': tk.StringVar(),
                            'non_essential_genes': tk.StringVar(),
                            'library_file': tk.StringVar(),
                            'target_samples': tk.StringVar(),
                            'reference_samples': tk.StringVar(),
                            'unwanted_columns': tk.StringVar(),
                            'unwanted_rows': tk.StringVar(),
                            'unwanted_row_substrings': tk.StringVar(),
                            'threshold_reads': tk.StringVar(),
                            'x_axis': tk.StringVar(),
                            'threshold_fdr': tk.StringVar(),
                            'top': tk.StringVar(),
                            'distribution_condition1': tk.StringVar(),
                            'distribution_condition2': tk.StringVar(),
                            'working_dir': tk.StringVar()
                            }

        # Create and add pages to the application
        for F in (StartPage, PageTwo, PageThree):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

        # Show the start page initially
        self.show_frame(StartPage)

    def show_frame(self, cont):
        # Show the specified frame
        frame = self.frames[cont]
        frame.tkraise()


# Run the application
if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
