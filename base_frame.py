import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog, messagebox
import os
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

    def add_trace(self, variable, label, check):
        """
        Adds a trace to monitor changes in a variable.

        Args:
            variable (tk.StringVar()): The variable to monitor.
            label (tk.Label): The label associated with the variable.
            check (str): The type of check to perform.
        """
        variable.trace_add("write", self.update_indicator(variable, label, check))

    def update_indicator(self, variable, indicator_label, check):
        """
        Updates the indicator based on the value of the monitored variable.

        Args:
            variable (tk.StringVar()): The variable being monitored.
            indicator_label (tk.Label): The label to update based on validation.
            check (str): The type of validation to perform.
        """
        var_value = variable.get()
        label_text = [k for k, v in self.indicator_labels.items() if v == indicator_label][0]

        if check == "non_empty":
            self.check_empty_file(indicator_label, label_text, var_value)
        elif check == "file_type":
            self.check_file_type(indicator_label, label_text, var_value)
        elif check == "folder_location":
            self.check_folder_location(indicator_label, label_text, var_value)

    def check_empty_file(self, indicator_label, label_text, var_value):
        """
        Checks if the provided value is empty and updates the indicator label accordingly.

        Args:
            indicator_label (tk.Label): The label to update.
            label_text (str): The text associated with the indicator.
            var_value (str): The value to check.
        """
        if var_value.strip() == "":
            indicator_label.config(text=" * Invalid")
            self.invalid_values[label_text] = True
        else:
            indicator_label.config(text="")
            self.invalid_values[label_text] = False

    def check_file_type(self, indicator_label, label_text, path):
        """
        Checks if the file type of the provided path is valid.

        Args:
            indicator_label (tk.Label): The label to update.
            label_text (str): The text associated with the indicator.
            path (str): The file path to check.
        """
        log_.debug(f"Checking path: {path}")
        if not os.path.exists(path) or not (path.endswith('.txt') or path.endswith('.csv')):
            indicator_label.config(text=" * Invalid")
            log_.debug(f"Path is invalid: {path}")
            self.invalid_file_types[label_text] = True
        else:
            indicator_label.config(text="")
            log_.debug(f"Path is valid: {path}")
            self.invalid_file_types[label_text] = False

    def check_folder_location(self, indicator_label, label_text, path):
        """
        Checks if the provided path is a valid folder location.

        Args:
            indicator_label (tk.Label): The label to update.
            label_text (str): The text associated with the indicator.
            path (str): The folder path to check.
        """
        log_.debug(f"Checking path: {path}")
        if not os.path.isdir(path):
            indicator_label.config(text=" * Invalid")
            log_.debug(f"Path is invalid: {path}")
            self.invalid_folder_locations[label_text] = True
        else:
            indicator_label.config(text="")
            log_.debug(f"Path is valid: {path}")
            self.invalid_folder_locations[label_text] = False

    @staticmethod
    def browse_files(output_variable):
        """
        Opens a file dialog to select a file and set the output variable.

        Args:
            output_variable (tk.StringVar): The variable to set with the selected file path.
        """
        file_path = filedialog.askopenfilename()
        output_variable.set(file_path)
