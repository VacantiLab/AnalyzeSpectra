from tkinter import Tk
from tkinter.filedialog import askdirectory

root = Tk()
root.withdraw()

directory = askdirectory() + '/'
root.update() #required so the directory request dialog box disappears and does not freeze
