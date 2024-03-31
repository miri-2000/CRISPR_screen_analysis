import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog

class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        # Create a container to hold multiple frames
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        # Create and add pages to the application
        for F in (StartPage, PageOne, PageTwo):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")

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
        header_label = ttk.Label(self, text="Welcome to the CRISPR screen analysis tool", font=("Helvetica", 18, "bold"))
        header_label.pack(pady=20)

        # Text label
        text_label = ttk.Label(self, text="Script for the identification of drug-gene interactions in CRISPR screens using DrugZ")
        text_label.pack(pady=10)

        # Define class attributes
        self.essential_genes = tk.StringVar()
        self.non_essential_genes = tk.StringVar()
        self.library_file = tk.StringVar()
        self.input_file = tk.StringVar()
        self.type = tk.StringVar()
        self.unwanted_columns = tk.StringVar()
        self.unwanted_rows = tk.StringVar()

        self.create_widgets()

        # Button to navigate to PageOne
        button1 = tk.Button(self, text="Next Page", command=lambda: controller.show_frame(PageOne))
        button1.pack()

    def create_widgets(self):
        # Essential genes entry
        ttk.Label(self, text="File with list of essential genes:").pack()
        tk.Entry(self, textvariable=self.essential_genes).pack()
        tk.Button(self, text="Browse", command=self.browse_essential_genes).pack()

        # Non-essential genes entry
        ttk.Label(self, text="File with list of non-essential genes:").pack()
        tk.Entry(self, textvariable=self.non_essential_genes).pack()
        tk.Button(self, text="Browse", command=self.browse_non_essential_genes).pack()

        # Library file entry
        ttk.Label(self, text="Library File:").pack()
        tk.Entry(self, textvariable=self.library_file).pack()
        tk.Button(self, text="Browse", command=self.browse_library_file).pack()

        # Input file entry
        ttk.Label(self, text="Input File:").pack()
        tk.Entry(self, textvariable=self.input_file).pack()
        tk.Button(self, text="Browse", command=self.browse_input_file).pack()

        # Replicate type checkboxes
        ttk.Label(self, text="Replicate type:").pack()
        tk.Checkbutton(self, text="Biological", variable=self.type, onvalue="biological", offvalue="").pack()
        tk.Checkbutton(self, text="Technical", variable=self.type, onvalue="technical", offvalue="").pack(),

        # Unwanted columns entry
        ttk.Label(self, text="Unwanted Columns (optional):").pack()
        tk.Entry(self, textvariable=self.unwanted_columns).pack()

        # Unwanted rows entry
        ttk.Label(self, text="Unwanted Rows (optional):").pack()
        tk.Entry(self, textvariable=self.unwanted_rows).pack()

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


class PageOne(tk.Frame):
    def __init__(self, root, controller):
        tk.Frame.__init__(self, root)

        # Define class attributes
        self.type = tk.StringVar()
        self.unwanted_row_substrings = tk.StringVar()
        self.threshold = tk.IntVar(value=0)
        self.target_samples = tk.StringVar()
        self.reference_samples = tk.StringVar()
        self.x_axis = tk.StringVar()
        self.threshold_fdr = tk.IntVar(value=0.25)
        self.top = tk.IntVar(value=15)
        self.distribution_condition1 = tk.StringVar()
        self.distribution_condition2 = tk.StringVar()

        self.create_widgets()

        # Button to navigate to ThirdPag
        button2 = tk.Button(self, text="Next Page", command=lambda: controller.show_frame(PageTwo))
        button2.pack()

    def create_widgets(self):
        # Essential genes entry
        # Unwanted row substrings entry
        ttk.Label(self, text="Unwanted Row Substrings (optional):").pack()
        tk.Entry(self, textvariable=self.unwanted_row_substrings).pack()

        # Min number of reads entry
        ttk.Label(self, text="Minimum required sum of reads per guide:").pack()
        tk.Entry(self, textvariable=self.threshold).pack()

        # Targets entry
        ttk.Label(self, text="Name of target samples:").pack()
        tk.Entry(self, textvariable=self.target_samples).pack()

        # References entry
        ttk.Label(self, text="Name of references samples:").pack()
        tk.Entry(self, textvariable=self.reference_samples).pack()

        # X-axis entry
        ttk.Label(self, text="X-axis value for significant plots:").pack()
        tk.Checkbutton(self, text="normZ value", variable=self.type, onvalue="normZ", offvalue="").pack()
        tk.Checkbutton(self, text="log2 fold-change", variable=self.type, onvalue="log2fc", offvalue="").pack(),

        # Gene significance threshold
        ttk.Label(self, text="Gene significance threshold:").pack()
        tk.Entry(self, textvariable=self.threshold).pack()

        # Top hits entry
        ttk.Label(self, text="Maximum number of hits:").pack()
        tk.Entry(self, textvariable=self.top).pack()

        # Positive Control entry
        ttk.Label(self, text="Positive Control:").pack()
        tk.Entry(self, textvariable=self.distribution_condition1).pack()

        # Negative Control entry
        ttk.Label(self, text="Negative Control:").pack()
        tk.Entry(self, textvariable=self.distribution_condition2).pack()


class PageTwo(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text="Page Two")
        label.pack(pady=10, padx=10)

        # Button to navigate to StartPage
        button1 = tk.Button(self, text="Go to Start Page",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

# Run the application
if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
