import tkinter as tk
import tkinter.ttk as ttk
from base_frame import BaseFrame
from tkinter import filedialog, messagebox
import os


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
        ttk.Button(self, text="Back", command=lambda: controller.show_frame(self.get_previous_class())).grid(row=9,
                                                                                                             column=1,
                                                                                                             pady=10)

        # Button to navigate to the next page
        ttk.Button(self, text="Next", command=lambda: controller.show_frame(self.get_next_class())).grid(row=9,
                                                                                                         column=2,
                                                                                                         pady=10)

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
        self.create_labeled_entry("Minimum required sum of reads/guide:",
                                  lambda: self.show_info("Minimum number of total reads/guide so that "
                                                         "the guide will not be discarded"),
                                  self.threshold_reads_entry, 5, "0")

    @staticmethod
    def get_previous_class():
        from page_one import StartPage
        # Return the previous class
        return StartPage

    @staticmethod
    def get_next_class():
        from page_three import PageThree
        # Return the next class
        return PageThree
