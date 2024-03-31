import tkinter as tk
from tkinter import filedialog

class ParameterGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Parameter Input")

        # Define class attributes
        self.essential_genes = tk.StringVar()
        self.non_essential_genes = tk.StringVar()
        self.library_file = tk.StringVar()
        self.input_file = tk.StringVar()
        self.type = tk.StringVar()
        self.unwanted_columns = tk.StringVar()
        self.unwanted_rows = tk.StringVar()
        self.unwanted_row_substrings = tk.StringVar()
        self.threshold = tk.IntVar(value=0)
        self.target_samples = tk.StringVar()
        self.reference_samples = tk.StringVar()
        self.x_axis = tk.StringVar()
        self.threshold_fdr = tk.IntVar(value=0.25)
        self.top = tk.IntVar(value=15)
        self.distribution_condition1 = tk.StringVar()
        self.distribution_condition2 = tk.StringVar()

        # Create and place widgets
        self.create_widgets()

    def create_widgets(self):
        # Essential genes entry
        tk.Label(self.root, text="File with list of essential genes:").pack()
        tk.Entry(self.root, textvariable=self.essential_genes).pack()
        tk.Button(self.root, text="Browse", command=self.browse_essential_genes).pack()

        # Non-essential genes entry
        tk.Label(self.root, text="File with list of non-essential genes:").pack()
        tk.Entry(self.root, textvariable=self.non_essential_genes).pack()
        tk.Button(self.root, text="Browse", command=self.browse_non_essential_genes).pack()

        # Library file entry
        tk.Label(self.root, text="Library File:").pack()
        tk.Entry(self.root, textvariable=self.library_file).pack()
        tk.Button(self.root, text="Browse", command=self.browse_library_file).pack()

        # Input file entry
        tk.Label(self.root, text="Input File:").pack()
        tk.Entry(self.root, textvariable=self.input_file).pack()
        tk.Button(self.root, text="Browse", command=self.browse_input_file).pack()

        # Replicate type checkboxes
        tk.Label(self.root, text="Replicate type:").pack()
        tk.Checkbutton(self.root, text="Biological", variable=self.type, onvalue="biological", offvalue="").pack()
        tk.Checkbutton(self.root, text="Technical", variable=self.type, onvalue="technical", offvalue="").pack(),

        # Unwanted columns entry
        tk.Label(self.root, text="Unwanted Columns (optional):").pack()
        tk.Entry(self.root, textvariable=self.unwanted_columns).pack()

        # Unwanted rows entry
        tk.Label(self.root, text="Unwanted Rows (optional):").pack()
        tk.Entry(self.root, textvariable=self.unwanted_rows).pack()

        # Unwanted row substrings entry
        tk.Label(self.root, text="Unwanted Row Substrings (optional):").pack()
        tk.Entry(self.root, textvariable=self.unwanted_row_substrings).pack()

        # Min number of reads entry
        tk.Label(self.root, text="Minimum required sum of reads per guide:").pack()
        tk.Entry(self.root, textvariable=self.threshold).pack()

        # Targets entry
        tk.Label(self.root, text="Name of target samples:").pack()
        tk.Entry(self.root, textvariable=self.target_samples).pack()

        # References entry
        tk.Label(self.root, text="Name of references samples:").pack()
        tk.Entry(self.root, textvariable=self.reference_samples).pack()

        # X-axis entry
        tk.Label(self.root, text="X-axis value for significant plots:").pack()
        tk.Checkbutton(self.root, text="normZ value", variable=self.type, onvalue="normZ", offvalue="").pack()
        tk.Checkbutton(self.root, text="log2 fold-change", variable=self.type, onvalue="log2fc", offvalue="").pack(),

        # Gene significance threshold
        tk.Label(self.root, text="Gene significance threshold:").pack()
        tk.Entry(self.root, textvariable=self.threshold).pack()

        # Top hits entry
        tk.Label(self.root, text="Maximum number of hits:").pack()
        tk.Entry(self.root, textvariable=self.top).pack()

        # Positive Control entry
        tk.Label(self.root, text="Positive Control:").pack()
        tk.Entry(self.root, textvariable=self.distribution_condition1).pack()

        # Negative Control entry
        tk.Label(self.root, text="Negative Control:").pack()
        tk.Entry(self.root, textvariable=self.distribution_condition2).pack()

        # Submit button
        tk.Button(self.root, text="Submit", command=self.submit).pack()

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

    # Add similar methods for other parameters...

    def submit(self):
        # Access and print class attributes
        print("Essential Genes:", self.essential_genes.get())
        print("Non-Essential Genes:", self.non_essential_genes.get())
        print("Library File:", self.library_file.get())
        # Print other parameters...

if __name__ == "__main__":
    root = tk.Tk()
    app = ParameterGUI(root)
    root.mainloop()
