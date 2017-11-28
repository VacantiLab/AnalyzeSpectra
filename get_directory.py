def get_directory(retrieve_directory_method):
    from tkinter import Tk #allows for asking for a directory through a GUI
    from tkinter.filedialog import askdirectory #allows for asking for a directory through a GUI

    if retrieve_directory_method == 'gui':
        #ask for the directory where the netCDF and library.txt files are
        root = Tk()
        root.withdraw() #closes the tkinter GUI window because the rest of the program is not run through the GUI
        file_directory = askdirectory() + '/'
        root.update() #required so the directory request dialog box disappears and does not freeze

    if retrieve_directory_method == 'manual':
        #manually enter file directory because tkinter is giving troubles when used with pdb.set_trace()
        file_directory = '/Users/nate/Dropbox/Research/Lehtio_Laboratory/Projects/GCMS_raw_data/2017_11_03_sknbe2_nmyckd_bt20_phen05uM.AIA/polar_metabolites_mox_tbdms/'

    return(file_directory)
