import tkinter as tk

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
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        # Header label
        header_label = tk.Label(self, text="Welcome to the CRISPR screen analysis tool", font=("Helvetica", 18, "bold"))
        header_label.pack(pady=20)

        # Text label
        text_label = tk.Label(self,
                              text="This is the start page.\nYou can write some instructions or information here.")
        text_label.pack(pady=10)


        # Button to navigate to PageOne
        button1 = tk.Button(self, text="Go to Page One",
                            command=lambda: controller.show_frame(PageOne))

        button1.pack()


class PageOne(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Page One")
        label.pack(pady=10, padx=10)

        # Button to navigate to StartPage
        button1 = tk.Button(self, text="Go to Start Page",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        # Button to navigate to PageTwo
        button2 = tk.Button(self, text="Go to Page Two",
                            command=lambda: controller.show_frame(PageTwo))
        button2.pack()

class PageTwo(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Page Two")
        label.pack(pady=10, padx=10)

        # Button to navigate to StartPage
        button1 = tk.Button(self, text="Go to Start Page",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

# Run the application
if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()
