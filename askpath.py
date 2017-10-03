from tkinter import *
from tkinter import filedialog

root = Tk()
root.withdraw()

root.update()
directory = filedialog.askdirectory() + '/'
