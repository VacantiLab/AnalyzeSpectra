from tkinter import *
from tkinter import filedialog

root = Tk()
root.withdraw()

directory = filedialog.askdirectory() + '/'
root.update() #required so the directory request dialog box disappears and does not freeze
