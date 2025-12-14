import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox
import logging as log

log.basicConfig(level=log.DEBUG)
log_ = log.getLogger(__name__)


class BaseFrame(ttk.Frame):
    """
    Base Class specifying common elements for all frames.
    This class provides common UI components and validation methods for derived frames.

    Attributes:
        controller: The main application controller that manages navigation.
        invalid_values (dict): Dictionary to hold invalid values.
        invalid_file_types (dict): Dictionary to hold invalid file types.
        invalid_folder_locations (dict): Dictionary to hold invalid folder locations.
        indicator_labels (dict): Dictionary to hold indicator label references.
    """

    def __init__(self, parent, controller):
        """
        Initializes the BaseFrame.

        Args:
            parent (tk.Widget): The parent widget for the frame.
            controller: The main application controller that manages navigation.
        """
        super().__init__(parent)
        self.controller = controller

        self.invalid_values = {}
        self.invalid_file_types = {}
        self.invalid_folder_locations = {}
        self.indicator_labels = {}

    def create_header(self, header_text):
        """Create header label for the page."""
        header_label = ttk.Label(self, text=header_text,
                                 font=("Helvetica", 18, "bold"))
        header_label.grid(row=0, column=0, columnspan=3, padx=20)

    def create_description(self, description_text):
        """Create description label for the page."""
        text_label = ttk.Label(self,
                               text=description_text, justify="center")
        text_label.grid(row=1, column=0, columnspan=3, pady=(5, 15), padx=20)

    def create_subheader(self, subtitle_text, row):
        text_label = ttk.Label(self,
                               text=subtitle_text,
                               font=("Helvetica", 10, "bold"))
        text_label.grid(row=row, column=0, sticky="w", pady=(3, 5), padx=20)

    def create_labeled_entry(self, label_text, command_action, textvariable, row, default_value=None):
        """
        Creates a labeled entry with an associated button and an optional threshold entry.

        Args:
            label_text (str): The text for the label.
            command_action (callable): The function to call when the button is pressed.
            textvariable: The variable to link with the entry.
            row (int): The row index in the grid layout.
            default_value (str, optional): The default value to insert into the threshold entry. Defaults to None.
        """
        ttk.Label(self, text=label_text).grid(row=row, column=0, sticky="w", padx=20)
        ttk.Button(self, text="?", command=command_action, width=1).grid(row=row, column=0, sticky="e")

        # Create the main entry widget
        if default_value is None:
            main_entry = ttk.Entry(self, textvariable=textvariable)
            main_entry.grid(row=row, column=1, sticky="w")

        # If a threshold entry is provided, configure it
        else:
            textvariable.insert(0, default_value)
            self.bind_entry_focus(textvariable, default_value)
            textvariable.grid(row=row, column=1, sticky="w")

    def create_browse_button(self, var, row):
        """Create a browse button for file selection."""
        ttk.Button(self, text="Browse", command=lambda: self.browse_files(var)).grid(row=row, column=2, sticky="w")

    def add_indicator_label(self, row, label_text):
        """Add an indicator label for validation feedback."""
        indicator_label = ttk.Label(self, text="")
        indicator_label.grid(row=row, column=3, sticky="w")
        self.indicator_labels[label_text] = indicator_label

    @staticmethod
    def show_info(label_text):
        """
        Displays an information message box.

        Args:
            label_text (str): The message to display in the info box.
        """
        messagebox.showinfo("Information", label_text)

    def bind_entry_focus(self, entry, default_value):
        entry.bind("<FocusIn>", lambda event: self.on_entry_click(entry, default_value))
        entry.bind("<FocusOut>", lambda event: self.on_entry_leave(entry, default_value))

    @staticmethod
    def on_entry_click(entry, default_text):
        """
        Clears the entry if it contains the default text.

        Args:
            entry (tk.Entry): The entry widget to check.
            default_text (str): The default text to clear.
        """
        if entry.get() == default_text:
            entry.delete(0, tk.END)

    @staticmethod
    def on_entry_leave(entry, default_text):
        """
        Restores the default text if the entry is empty.

        Args:
            entry (tk.Entry): The entry widget to check.
            default_text (str): The default text to restore.
        """
        if entry.get() == "":
            entry.insert(0, default_text)

    @staticmethod
    def browse_files(output_variable):
        """
        Opens a file dialog to select a file and set the output variable.

        Args:
            output_variable (tk.StringVar): The variable to set with the selected file path.
        """
        file_path = filedialog.askopenfilename()
        output_variable.set(file_path)
