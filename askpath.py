from tkinter import *
from tkinter import filedialog

root = Tk()
root.withdraw()

directory = filedialog.askdirectory() + '/'
root.update()
