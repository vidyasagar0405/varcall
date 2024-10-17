import tkinter as tk
from tkinter import filedialog

def open_file_dialog():
    files = filedialog.askopenfilenames()
    if files:
        print("OK")
    return files

if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()

    files = open_file_dialog()
    if files:
        print("Selected files:")
        file_string = []
        for file in files:
            file_string.append(file)
        print(file_string)
        print()
        print(" ".join(file_string))
    else:
        print("FAILED")
