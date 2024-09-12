import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox
from run_analysis import CRISPR_screen_analysis
import os


class BaseFrame(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

        self.invalid_paths = {}  # Dictionary to hold invalid paths
        self.invalid_values = {}  # Dictionary to hold invalid paths
        self.indicator_labels = {}
        # self.threshold = ttk.Entry(self)

    def create_labeled_entry(self, label_text, command_action, textvariable, row):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w")
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0, sticky="e")
        ttk.Entry(self, textvariable=textvariable).grid(row=row, column=1, sticky="w")

    def create_labeled_entry_with_threshold(self, label_text, row, threshold_entry, default_value):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=(0, 20))
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=row, column=0, sticky="e")
        threshold_entry.insert(0, default_value)
        threshold_entry.bind("<FocusIn>", lambda event: self.on_entry_click(event, threshold_entry, default_value))
        threshold_entry.bind("<FocusOut>", lambda event: self.on_entry_leave(event, threshold_entry, default_value))
        threshold_entry.grid(row=row, column=1)

    def show_info(self, label_text):
        messagebox.showinfo("Information", f"Information about {label_text}.")

    def on_entry_click(self, event, entry, default_text):
        if entry.get() == default_text:
            entry.delete(0, tk.END)

    def on_entry_leave(self, event, entry, default_text):
        if entry.get() == "":
            entry.insert(0, default_text)

    def add_trace(self, value, label, check_path=True):
        value.trace_add("write", self.update_indicator(value, label, check_path))

    def update_indicator(self, var, indicator_label, check_path):
        var_value = var.get()
        label_text = [k for k, v in self.indicator_labels.items() if v == indicator_label][0]

        if check_path:
            self.check_path_file(indicator_label, label_text, var_value)
        else:
            self.check_empty_file(indicator_label, label_text, var_value)

    def check_path_file(self, indicator_label, label_text, path):
        print(f"Checking path: {path}")  # Debug print
        if not os.path.exists(path) or not (path.endswith('.txt') or path.endswith('.csv')):
            indicator_label.config(text=" * Invalid")
            print(f"Path is invalid: {path}")  # Debug print
            self.invalid_paths[label_text] = True
        else:
            indicator_label.config(text="")
            print(f"Path is valid: {path}")  # Debug print
            self.invalid_paths[label_text] = False

    def check_empty_file(self, indicator_label, label_text, var_value):
        if var_value.strip() == "":
            indicator_label.config(text=" * Invalid")
            self.invalid_values[label_text] = True
        else:
            indicator_label.config(text="")
            self.invalid_values[label_text] = False

    # def validate_int_input(value_if_allowed):
    #     if value_if_allowed == "":
    #         try:
    #             test = int(value_if_allowed).isdigit()  # Allow empty string or digits only
    #             return True
    #         except ValueError:
    #             pass
    #     return False


class StartPage(BaseFrame):
    # class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Welcome to the CRISPR screen analysis tool",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=4)

        # Text label
        text_label = ttk.Label(self,
                               text="Welcome to the CRISPR screen analysis tool for the identification\n"
                                    "of drug-gene interactions using DrugZ! This tool is designed to analyze\n"
                                    "CRISPR screens and identify potential drug-gene interactions. To get started,\n"
                                    "please provide the following parameters:", justify="center")
        text_label.grid(row=1, column=0, columnspan=4, pady=(5, 15))

        self.controller = controller
        # self.input_file = tk.StringVar()
        # self.essential_genes = tk.StringVar()
        # self.non_essential_genes = tk.StringVar()
        # self.library_file = tk.StringVar()
        # self.type = tk.StringVar()
        # self.target_samples = tk.StringVar()
        # self.reference_samples = tk.StringVar()

        # Define the file paths and associated labels
        self.file_paths = {
            # "Screen Result File": self.input_file,
            # "Essential Genes File": self.essential_genes,
            # "Non-Essential Genes File": self.non_essential_genes,
            # "Library File": self.library_file
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
                                      lambda l=label_text: self.show_info(l),
                                      var, i)

            ttk.Button(self, text="Browse",
                       command=lambda var=self.file_paths[label_text]: self.browse_files(var)).grid(row=i, column=2,
                                                                                                    sticky="w")

            # Add a label or star indicator next to each entry
            indicator_label = ttk.Label(self, text="")
            indicator_label.grid(row=i, column=3, sticky="w")
            self.indicator_labels[label_text] = indicator_label

    def browse_files(self, output_variable):
        file_path = filedialog.askopenfilename()
        output_variable.set(file_path)

    def validate_and_proceed(self):
        self.invalid_paths.clear()
        self.invalid_values.clear()

        self.add_trace(self.controller.shared_data["target_samples"], self.indicator_labels["Name of target samples:"],
                       False)
        self.add_trace(self.controller.shared_data["reference_samples"],
                       self.indicator_labels["Name of reference samples:"], False)

        invalid_value_labels = [k for k, v in self.invalid_values.items() if v is True]
        if invalid_value_labels:
            messagebox.showwarning("Missing information",
                                   f"The following fields cannot remain empty: {', '.join(invalid_value_labels)}. "
                                   f"Please specify the parameters.")

        for label_text, var in self.file_paths.items():
            self.add_trace(var, self.indicator_labels[label_text])

        invalid_path_labels = [k for k, v in self.invalid_paths.items() if v is True]
        if invalid_path_labels:
            messagebox.showwarning("Invalid Paths",
                                   f"The following paths are invalid: {', '.join(invalid_path_labels)}. The "
                                   f"file path must be accurate and lead to a text or csv file.")

        if not invalid_value_labels and not invalid_path_labels:
            self.controller.show_frame(PageTwo)


class PageTwo(BaseFrame):

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Data Preparation Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0)

        # Title label
        text_label = ttk.Label(self,
                               text="Optional data preparation settings (default values shown in line)")
        # ,justify="center")
        text_label.grid(row=1, column=0, columnspan=2, pady=(5, 15))

        # Define class attributes
        # self.unwanted_columns = tk.StringVar()
        # self.unwanted_rows = tk.StringVar()
        # self.unwanted_row_substrings = tk.StringVar()
        self.threshold_reads_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_reads"])

        self.create_widgets()

        # Button to navigate to the previous page
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(StartPage)).grid(row=9, column=1, pady=10)

        # Button to navigate to the next page
        ttk.Button(self, text="Next", command=lambda: controller.show_frame(PageThree)).grid(row=9, column=2, pady=10)

    def create_widgets(self):
        # Unwanted columns entry
        self.create_labeled_entry("Unwanted Columns:", self.show_info, self.controller.shared_data["unwanted_columns"],
                                  2)

        # Unwanted rows entry
        self.create_labeled_entry("Unwanted Rows:", self.show_info, self.controller.shared_data["unwanted_rows"], 3)

        # Unwanted row substrings entry
        self.create_labeled_entry("Unwanted Row Substrings:", self.show_info,
                                  self.controller.shared_data["unwanted_row_substrings"], 4)

        # Min number of reads entry
        self.create_labeled_entry_with_threshold("Minimum required sum of reads/guide:", 5,
                                                 self.threshold_reads_entry, "0")


class PageThree(BaseFrame):
    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Visualization Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=3)

        # Text label
        text_label = ttk.Label(self,
                               text="Optional result visualization settings (default values shown in line)\n",
                               justify="center")
        text_label.grid(row=1, column=0, columnspan=3, pady=(5, 0))

        # Define class attributes
        # self.x_axis = tk.StringVar()
        # self.top = ttk.Entry(self)
        self.top_entry = ttk.Entry(self, textvariable=self.controller.shared_data["top"])
        self.threshold_fdr_entry = ttk.Entry(self, textvariable=self.controller.shared_data["threshold_fdr"])
        # self.distribution_condition1 = tk.StringVar()
        # self.distribution_condition2 = tk.StringVar()
        self.condition = {"Positive Control": self.controller.shared_data["distribution_condition1"],
                          "Negative Control": self.controller.shared_data["distribution_condition2"]}
        # self.working_dir = tk.StringVar()

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
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=3, column=0, sticky="e")
        ttk.Radiobutton(self, text="normZ", variable=self.controller.shared_data["x_axis"], value="normZ").grid(
            row=3, column=1, sticky="w")
        ttk.Radiobutton(self, text="log2 fold-change", variable=self.controller.shared_data["x_axis"],
                        value="log2 fold-change").grid(row=3,
                                                       column=2,
                                                       sticky="w")

        # Gene significance threshold
        self.create_labeled_entry_with_threshold("Significance threshold:", 4, self.threshold_fdr_entry, "0.25")

        # Top hits entry
        # self.create_labeled_entry_with_threshold("Number of hits per plot:", 5, self.top, "15")
        self.create_labeled_entry_with_threshold("Number of hits per plot:", 5, self.top_entry, "15")

        # Title label
        text_label = ttk.Label(self,
                               text="Distribution Plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=6, column=0, sticky="w", pady=(10, 5))

        # Positive Control entry
        self.create_labeled_entry("Positive Control:", self.show_info,
                                  self.controller.shared_data["distribution_condition1"], 7)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=7, column=2, sticky="w")
        self.indicator_labels["Positive Control"] = indicator_label

        # Negative Control entry
        self.create_labeled_entry("Negative Control:", self.show_info,
                                  self.controller.shared_data["distribution_condition2"], 8)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=8, column=2, sticky="w")
        self.indicator_labels["Negative Control"] = indicator_label

        self.create_labeled_entry("Result storage location:", self.show_info,
                                  self.controller.shared_data["working_dir"], 10)
        # self.working_dir = self.working_dir.get()

    def validate_and_proceed(self):
        self.invalid_values.clear()

        for label_text, var in self.indicator_labels.items():
            self.add_trace(self.condition[label_text], var, False)

        invalid_path_labels = [k for k, v in self.invalid_values.items() if v is True]
        if invalid_path_labels:
            messagebox.showwarning("Missing information",
                                   f"The following fields cannot remain empty: {', '.join(invalid_path_labels)}. "
                                   f"Please specify the parameters.")
        else:
            self.start_computation()

    # def convert_stringvars_to_str(self):
    #     for attr_name, attr_value in self.__dict__.items():
    #         if isinstance(attr_value, tk.StringVar):
    #             self.__dict__[attr_name] = attr_value.get()
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


# Run the application
class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # Create a container to hold multiple frames
        container = tk.Frame(self)
        container.grid(row=0, column=0, sticky="nsew")
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
