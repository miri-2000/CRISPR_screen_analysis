import tkinter as tk
from page_one import StartPage
from page_two import PageTwo
from page_three import PageThree


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
