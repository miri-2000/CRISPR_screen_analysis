import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox

import run_analysis
import os


class BaseFrame(ttk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)
        self.controller = controller

    def create_labeled_entry(self, label_text, command_action, textvariable, row):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w")
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0, sticky="e")
        ttk.Entry(self, textvariable=textvariable).grid(row=row, column=1, sticky="w")

    def create_labeled_entry_with_threshold(self, label_text, row, threshold_entry):
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=(0, 20))
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=row, column=0, sticky="e")
        threshold_entry.insert(0, "0")
        threshold_entry.bind("<FocusIn>", lambda event: self.on_entry_click(event, threshold_entry, "0"))
        threshold_entry.bind("<FocusOut>", lambda event: self.on_entry_leave(event, threshold_entry, "0"))
        threshold_entry.grid(row=row, column=1)

    def show_info(self, label_text):
        messagebox.showinfo("Information", f"Information about {label_text}.")

    def on_entry_click(self, event, entry, default_text):
        if entry.get() == default_text:
            entry.delete(0, tk.END)

    def on_entry_leave(self, event, entry, default_text):
        if entry.get() == "":
            entry.insert(0, default_text)


class StartPage(BaseFrame):
    # class StartPage(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent, controller)
        # def __init__(self, root, controller):
        #     tk.Frame.__init__(self, root)

        # Header label
        header_label = ttk.Label(self, text="Welcome to the CRISPR screen analysis tool",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=4)
        # header_label.pack()

        # Text label
        text_label = ttk.Label(self,
                               text="Welcome to the CRISPR screen analysis tool for the identification\n"
                                    "of drug-gene interactions using DrugZ! This tool is designed to analyze\n"
                                    "CRISPR screens and identify potential drug-gene interactions. To get started,\n"
                                    "please provide the following parameters:", justify="center")
        text_label.grid(row=1, column=0, columnspan=4, pady=(5, 15))
        # text_label.pack()

        self.controller = controller
        self.input_file = tk.StringVar()
        self.essential_genes = tk.StringVar()
        self.non_essential_genes = tk.StringVar()
        self.library_file = tk.StringVar()
        self.type = tk.StringVar()
        self.target_samples = tk.StringVar()
        self.reference_samples = tk.StringVar()

        # Define the file paths and associated labels
        self.file_paths = {
            "Screen Result File": self.input_file,
            "Essential Genes File": self.essential_genes,
            "Non-Essential Genes File": self.non_essential_genes,
            "Library File": self.library_file
        }

        self.invalid_paths = {}  # Dictionary to hold invalid paths
        self.indicator_labels = {}  # Dictionary to hold indicator labels
        self.tracers_added = False  # Flag to check if tracers are added

        self.create_widgets()

        # Button to navigate to PreperationPage
        button1 = ttk.Button(self, text="Next Page", command=self.validate_and_proceed)
        button1.grid(row=9, column=1, pady=20)

    def create_widgets(self):
        # Targets entry
        # ttk.Label(self, text="Name of target samples:").grid(row=2, column=0, sticky="w")
        # ttk.Button(self, text="?", command=lambda: self.show_info(
        #     "Treated conditions (e.g. with exposure to drugs,...) to be compared against baseline conditions"),
        #            width=1).grid(row=2, column=0, sticky="e")
        # ttk.Entry(self, textvariable=self.target_samples).grid(row=2, column=1, sticky="w")
        self.create_labeled_entry(label_text="Name of target samples:", command_action=lambda: self.show_info(
            "Treated conditions (e.g. with exposure to drugs,...) to be compared against baseline conditions"),
                                  textvariable=self.target_samples, row=2)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=2, column=1, sticky="w")
        self.indicator_labels["Name of target samples:"] = indicator_label

        # References entry
        # ttk.Label(self, text="Name of references samples:").grid(row=3, column=0, sticky="w")
        # ttk.Button(self, text="?", command=lambda: self.show_info("Baseline conditions (e.g. with no added drugs)"),
        #            width=1).grid(row=3, column=0, sticky="e")
        # ttk.Entry(self, textvariable=self.reference_samples).grid(row=3, column=1, sticky="w")
        self.create_labeled_entry("Name of reference samples:",
                                  lambda: self.show_info("Baseline conditions (e.g. with no added drugs)"),
                                  self.reference_samples, 3)

        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=3, column=1, sticky="w")
        self.indicator_labels["Name of reference samples:"] = indicator_label

        for i, (label_text, var) in enumerate(self.file_paths.items(), start=4):
            self.create_labeled_entry(f"{label_text}:",
                                      lambda l=label_text: self.show_info(l),
                                      var, i)
            # ttk.Label(self, text=f"{label_text}:").grid(row=i, column=0, sticky="w")
            # ttk.Button(self, text="?", command=lambda l=label_text: self.show_info(l), width=1).grid(row=i, column=0,
            #                                                                                          sticky="e")

            # ttk.Entry(self, textvariable=var).grid(row=i, column=1, sticky="w")

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

    def add_trace(self, var, indicator_label):
        var.trace_add("write", self.update_indicator(indicator_label, var))
        # var.trace_add("write", lambda *args: self.update_indicator(indicator_label, var))

    def update_indicator(self, indicator_label, var):
        path = var.get()
        label_text = [k for k, v in self.indicator_labels.items() if v == indicator_label][0]
        print(f"Checking path: {path}")  # Debug print
        if not os.path.exists(path) or not (path.endswith('.txt') or path.endswith('.csv')):
            indicator_label.config(text=" * Invalid")
            print(f"Path is invalid: {path}")  # Debug print
            self.invalid_paths[label_text] = True
        else:
            indicator_label.config(text="")
            print(f"Path is valid: {path}")  # Debug print
            self.invalid_paths[label_text] = False

    def validate_and_proceed(self):
        self.invalid_paths.clear()
        # invalid_paths = []

        for label_text, var in self.file_paths.items():
            path = var.get()
            self.add_trace(self.file_paths[label_text], self.indicator_labels[label_text])

        invalid_path_labels = [k for k, v in self.invalid_paths.items() if v is True]
        if invalid_path_labels:
            messagebox.showwarning("Invalid Paths",
                                   f"The following paths are invalid: {', '.join(invalid_path_labels)}. The "
                                   f"file path must be accurate and lead to a text or csv file.")
        else:
            self.controller.show_frame(PageTwo)


class PageTwo(BaseFrame):
    # class PageTwo(tk.Frame):
    # def __init__(self, root, controller):
    #     tk.Frame.__init__(self, root)

    def __init__(self, parent, controller):
        super().__init__(parent, controller)

        # Header label
        header_label = ttk.Label(self, text="Data Preparation Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0)

        # Text label
        text_label = ttk.Label(self,
                               text="Optional data preparation settings (default values shown in line)",
                               justify="center")
        text_label.grid(row=1, column=0, columnspan=2, pady=(5, 15))

        # Define class attributes
        self.unwanted_columns = tk.StringVar()
        self.unwanted_rows = tk.StringVar()
        self.unwanted_row_substrings = tk.StringVar()
        self.threshold = ttk.Entry(self)
        # self.threshold = tk.IntVar(value=0)

        self.create_widgets()

        # Button to navigate to the previous page
        button1 = ttk.Button(self, text="Back", command=lambda: controller.show_frame(StartPage))
        button1.grid(row=9, column=1, pady=10)

        # Button to navigate to the next page
        button2 = ttk.Button(self, text="Next", command=lambda: controller.show_frame(PageThree))
        button2.grid(row=9, column=2, pady=10)

    def create_widgets(self):
        # Unwanted columns entry
        ttk.Label(self, text="Unwanted Columns:").grid(row=2, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=2, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.unwanted_columns).grid(row=2, column=1)

        # Unwanted rows entry
        ttk.Label(self, text="Unwanted Rows:").grid(row=3, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=3, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.unwanted_rows).grid(row=3, column=1)

        # Unwanted row substrings entry
        ttk.Label(self, text="Unwanted Row Substrings:").grid(row=4, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=4, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.unwanted_row_substrings).grid(row=4, column=1, sticky="w")

        # Min number of reads entry
        ttk.Label(self, text="Minimum required sum of reads/guide:").grid(row=5, column=0, sticky="w", padx=(0, 20))
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=5, column=0, sticky="e")
        self.threshold.insert(0, "0")
        self.threshold.bind("<FocusIn>", lambda event: self.on_entry_click(event, self.threshold, "0"))
        self.threshold.bind("<FocusOut>",
                            lambda event: self.on_entry_leave(event, self.threshold, "0"))
        self.threshold.grid(row=5, column=1)


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
        self.x_axis_new = tk.StringVar()
        self.x_axis = tk.StringVar()
        self.threshold_fdr = ttk.Entry(self)
        self.top = ttk.Entry(self)
        self.distribution_condition1 = tk.StringVar()
        self.distribution_condition2 = tk.StringVar()
        self.condition = {"Positive Control": self.distribution_condition1,
                          "Negative Control": self.distribution_condition2}
        self.invalid_paths = {}
        self.indicator_labels = {}

        self.create_widgets()

        # Button to navigate to the previous page
        button1 = ttk.Button(self, text="Back", command=lambda: controller.show_frame(PageTwo))
        button1.grid(row=9, column=1, pady=10)

        # Button to navigate to the next page
        button2 = ttk.Button(self, text="Start Computation", command=self.validate_and_proceed)
        # button2 = ttk.Button(self, text="Start Computation", command=self.start_computation)
        button2.grid(row=9, column=2, pady=10)

        # Button to navigate to ThirdPag
        # button1 = ttk.Button(self, text="Start computation", command=self.start_computation)
        # button1.grid()

    def create_widgets(self):
        # Text label
        text_label = ttk.Label(self,
                               text="Gene significance plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=2, column=0, sticky="w", pady=(3, 5))

        # Set the initial value of the selected option
        self.x_axis_new.set("normZ")  # Make Option 1 the default

        # Create radiobuttons with different text and values
        ttk.Label(self, text="X-axis value:").grid(row=3, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=3, column=0, sticky="e")
        ttk.Radiobutton(self, text="normZ", variable=self.x_axis_new, value="normZ").grid(row=3, column=1, sticky="w")
        ttk.Radiobutton(self, text="log2 fold-change", variable=self.x_axis_new, value="log2 fold-change").grid(row=3,
                                                                                                                column=2,
                                                                                                                sticky="w")
        # ttk.Checkbutton(self, text="normZ value", variable=self.x_axis, onvalue="normZ", offvalue="").grid(row=2, column=1)
        # ttk.Checkbutton(self, text="log2 fold-change", variable=self.x_axis, onvalue="log2fc", offvalue="").grid(row=2, column=2)

        # Gene significance threshold
        ttk.Label(self, text="Significance threshold:").grid(row=4, column=0, sticky="w")
        self.threshold_fdr.insert(0, "0.25")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=4, column=0, sticky="e")
        self.threshold_fdr.bind("<FocusIn>", lambda event: self.on_entry_click(event, self.threshold_fdr, "0.25"))
        self.threshold_fdr.bind("<FocusOut>",
                                lambda event: self.on_entry_leave(event, self.threshold_fdr, "0.25"))
        self.threshold_fdr.grid(row=4, column=1)

        # Top hits entry
        ttk.Label(self, text="Number of hits per plot:").grid(row=5, column=0, sticky="w")
        self.top.insert(0, "15")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=5, column=0, sticky="e")
        self.top.bind("<FocusIn>", lambda event: self.on_entry_click(event, self.top, "15"))
        self.top.bind("<FocusOut>",
                      lambda event: self.on_entry_leave(event, self.top, "15"))
        self.top.grid(row=5, column=1)

        # Text label
        text_label = ttk.Label(self,
                               text="Distribution Plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=6, column=0, sticky="w", pady=(10, 5))

        # Positive Control entry
        ttk.Label(self, text="Positive Control:").grid(row=7, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=7, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.distribution_condition1).grid(row=7, column=1)
        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=7, column=2, sticky="w")
        self.indicator_labels["Positive Control"] = indicator_label

        # Negative Control entry
        ttk.Label(self, text="Negative Control:").grid(row=8, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=8, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.distribution_condition2).grid(row=8, column=1)
        # Add a label or star indicator next to each entry
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=8, column=2, sticky="w")
        self.indicator_labels["Negative Control"] = indicator_label

    def add_trace(self, value, label):
        value.trace_add("write", self.update_indicator(label, value))
        # var.trace_add("write", lambda *args: self.update_indicator(indicator_label, var))

    def update_indicator(self, indicator_label, var):
        label_text = [k for k, v in self.indicator_labels.items() if v == indicator_label][0]
        if var.get().strip() == "":
            indicator_label.config(text=" * Invalid")
            self.invalid_paths[label_text] = True
        else:
            indicator_label.config(text="")
            self.invalid_paths[label_text] = False

    def validate_and_proceed(self):
        for label_text, var in self.indicator_labels.items():
            self.add_trace(self.condition[label_text], var)

        invalid_path_labels = [k for k, v in self.invalid_paths.items() if v is True]
        if invalid_path_labels:
            messagebox.showwarning("Missing information",
                                   f"The following fields cannot remain empty: {', '.join(invalid_path_labels)}. "
                                   f"Please specify the parameters.")
        else:
            self.start_computation()

    def start_computation(self):
        # Call the function from another file with self as the argument
        run_analysis.CRISPR_screen_analysis(self)


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
