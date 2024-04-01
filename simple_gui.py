import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox

import run_analysis
from run_analysis import CRISPR_screen_analysis


class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # Create a container to hold multiple frames
        container = tk.Frame(self)
        # container.pack(side="top", fill="both", expand=True)
        container.grid(row=0, column=0, sticky="nsew")
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        # Create and add pages to the application
        for F in (StartPage, PreperationPage, VisualizationPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
            # frame.pack(side="top", fill="both", expand=True)
            # frame.grid_rowconfigure(0, weight=1)
            # frame.grid_columnconfigure(0, weight=1)

        # Show the start page initially
        self.show_frame(StartPage)

    def show_frame(self, cont):
        # Show the specified frame
        frame = self.frames[cont]
        frame.tkraise()


# Define the pages of the application
class StartPage(tk.Frame):
    def __init__(self, root, controller):
        tk.Frame.__init__(self, root)

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

        # Define class attributes
        self.input_file = tk.StringVar()
        self.essential_genes = tk.StringVar()
        self.non_essential_genes = tk.StringVar()
        self.library_file = tk.StringVar()
        self.type = tk.StringVar()
        self.target_samples = tk.StringVar()
        self.reference_samples = tk.StringVar()

        self.create_widgets()

        # Button to navigate to PageOne
        button1 = ttk.Button(self, text="Next Page", command=lambda: controller.show_frame(PreperationPage))
        button1.grid(row=9, column=1, pady=20)

    def create_widgets(self):
        ttk.Label(self, text="Screen Result File:").grid(row=2, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=2, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.input_file).grid(row=2, column=1)
        ttk.Button(self, text="Browse", command=self.browse_input_file).grid(row=2, column=2)

        # Targets entry
        ttk.Label(self, text="Name of target samples:").grid(row=3, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=3, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.target_samples).grid(row=3, column=1)

        # References entry
        ttk.Label(self, text="Name of references samples:").grid(row=4, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=4, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.reference_samples).grid(row=4, column=1)

        ttk.Label(self, text="File with list of essential genes:").grid(row=5, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=5, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.essential_genes).grid(row=5, column=1)
        ttk.Button(self, text="Browse", command=self.browse_essential_genes).grid(row=5, column=2)

        ttk.Label(self, text="File with list of non-essential genes:").grid(row=6, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=6, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.non_essential_genes).grid(row=6, column=1)
        ttk.Button(self, text="Browse", command=self.browse_non_essential_genes).grid(row=6, column=2)

        ttk.Label(self, text="Library File:").grid(row=7, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=7, column=0, sticky="e")
        ttk.Entry(self, textvariable=self.library_file).grid(row=7, column=1)
        ttk.Button(self, text="Browse", command=self.browse_library_file).grid(row=7, column=2)

        # Replicate type checkboxes
        ttk.Label(self, text="Replicate type:").grid(row=8, column=0, sticky="w")
        ttk.Button(self, text="?", command=self.show_info, width=1).grid(row=8, column=0, sticky="e")
        ttk.Checkbutton(self, text="Biological", variable=self.type, onvalue="biological", offvalue="").grid(row=8,
                                                                                                             column=1)
        ttk.Checkbutton(self, text="Technical", variable=self.type, onvalue="technical", offvalue="").grid(row=8,
                                                                                                           column=2)

    def browse_essential_genes(self):
        file_path = filedialog.askopenfilename()
        self.essential_genes.set(file_path)

    def browse_non_essential_genes(self):
        file_path = filedialog.askopenfilename()
        self.non_essential_genes.set(file_path)

    def browse_library_file(self):
        file_path = filedialog.askopenfilename()
        self.library_file.set(file_path)

    def browse_input_file(self):
        file_path = filedialog.askopenfilename()
        self.input_file.set(file_path)

    def show_info(self):
        messagebox.showinfo("Information", "This is a tooltip explaining what the user needs to enter.")


class PreperationPage(tk.Frame):
    def __init__(self, root, controller):
        tk.Frame.__init__(self, root)

        # Header label
        header_label = ttk.Label(self, text="Data Preparation Settings",
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=2)

        # Text label
        text_label = ttk.Label(self,
                               text="Optional data preparation settings (default values shown in line)",
                               justify="center")
        text_label.grid(row=1, column=0, columnspan=2, pady=(5, 15))

        # Define class attributes
        self.unwanted_columns = tk.StringVar()
        self.unwanted_rows = tk.StringVar()
        self.unwanted_row_substrings = tk.StringVar()
        self.threshold = tk.IntVar(value=0)

        self.create_widgets()

        # Button to navigate to PageOne
        button1 = ttk.Button(self, text="Next Page", command=lambda: controller.show_frame(VisualizationPage))
        button1.grid(row=9, column=1, pady=10)

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
        ttk.Entry(self, textvariable=self.threshold).grid(row=5, column=1, sticky="w")

    def show_info(self):
        messagebox.showinfo("Information", "This is a tooltip explaining what the user needs to enter.")


class VisualizationPage(tk.Frame):
    def __init__(self, root, controller):
        tk.Frame.__init__(self, root)

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

        self.create_widgets()

        # Button to navigate to ThirdPag
        button1 = ttk.Button(self, text="Start computation", command=self.start_computation)
        button1.grid()

    def create_widgets(self):
        # Text label
        text_label = ttk.Label(self,
                               text="Gene significance plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=2, column=0,sticky="w",pady=(3,5))


        # Set the initial value of the selected option
        self.x_axis_new.set("normZ")  # Make Option 1 the default

        # Create radiobuttons with different text and values
        ttk.Label(self, text="X-axis value for gene significance plots:").grid(row=3, column=0, sticky="w")
        ttk.Radiobutton(self, text="normZ", variable=self.x_axis_new, value="normZ").grid(row=3, column=1, sticky="w")
        ttk.Radiobutton(self, text="log2 fold-change", variable=self.x_axis_new, value="log2 fold-change").grid(row=3,
                                                                                                                column=2, sticky="w")
        # ttk.Checkbutton(self, text="normZ value", variable=self.x_axis, onvalue="normZ", offvalue="").grid(row=2, column=1)
        # ttk.Checkbutton(self, text="log2 fold-change", variable=self.x_axis, onvalue="log2fc", offvalue="").grid(row=2, column=2)

        # Gene significance threshold
        ttk.Label(self, text="Gene significance threshold:").grid(row=4, column=0, sticky="w")
        self.threshold_fdr.insert(0, "0.25")
        self.threshold_fdr.bind("<FocusIn>", lambda event: self.on_entry_click(event, self.threshold_fdr, "0.25"))
        self.threshold_fdr.bind("<FocusOut>",
                                lambda event: self.on_entry_leave(event, self.threshold_fdr, "0.25"))
        self.threshold_fdr.grid(row=4, column=1)

        # Top hits entry
        ttk.Label(self, text="Maximum number of hits:").grid(row=5, column=0, sticky="w")
        self.top.insert(0, "15")
        self.top.bind("<FocusIn>", lambda event: self.on_entry_click(event, self.top, "15"))
        self.top.bind("<FocusOut>",
                                lambda event: self.on_entry_leave(event, self.top, "15"))
        self.top.grid(row=5, column=1)

        # Text label
        text_label = ttk.Label(self,
                               text="Distribution Plots",
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=6, column=0,sticky="w",pady=(10,5))

        # Positive Control entry
        ttk.Label(self, text="Positive Control:").grid(row=7, column=0, sticky="w")
        ttk.Entry(self, textvariable=self.distribution_condition1).grid(row=7, column=1)

        # Negative Control entry
        ttk.Label(self, text="Negative Control:").grid(row=8, column=0, sticky="w")
        ttk.Entry(self, textvariable=self.distribution_condition2).grid(row=8, column=1)

    def on_entry_click(self, event, entry_widget, default_text):
        if entry_widget.get() == default_text:
            entry_widget.delete(0, tk.END)

    def on_entry_leave(self, event, entry_widget, default_text):
        if entry_widget.get() == "":
            entry_widget.insert(0, default_text)

    def start_computation(self):
        # Call the function from another file with self as the argument
        run_analysis.CRISPR_screen_analysis(self)


# Run the application
if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
