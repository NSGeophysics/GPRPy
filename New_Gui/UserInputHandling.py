import tkinter as tk

class UserInputHandling:
    def __init__(self):
        self.delimiter = None

    def getDelimiter(self):
        # Method to prompt the user to select a delimiter (comma or tab)
        commaQuery = tk.Toplevel()
        commaQuery.title("Comma or tab separated?")
        text = tk.Label(commaQuery, text="Is this a comma- or tab-separated file?", fg='red')
        text.pack(padx=10, pady=10)
        commaButton = tk.Button(commaQuery, text="comma", width=10,
                                command=lambda: [self.setComma(), commaQuery.destroy()])
        commaButton.pack(side="left")
        tabButton = tk.Button(commaQuery, text="tab", width=10,
                              command=lambda: [self.setTab(), commaQuery.destroy()])
        tabButton.pack(side="right")
        commaQuery.wait_window()

    def setComma(self):
        self.delimiter = ','
        print("Delimiter set to comma")

    def setTab(self):
        self.delimiter = '\t'
        print("Delimiter set to tab")