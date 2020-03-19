def GetFileBatch(file_directory):
    # opens the batch_to_alkane.txt in the data location directory

    import pandas as pd
    import pdb
    import copy
    import numpy as np

    file_to_batch_map_directory = file_directory + 'file_to_batch.txt'
    file_name_to_batch = pd.read_table(file_to_batch_map_directory, sep="\t", header=0)
    batches = pd.unique(file_name_to_batch['batch'])

    batch_to_alkane_map_directory = file_directory + 'batch_to_alkane.txt'
    batch_to_alkane = pd.read_table(batch_to_alkane_map_directory, sep="\t", header=0)
    alkane_files = pd.unique(batch_to_alkane['alkane_file'])

    return(file_name_to_batch,batch_to_alkane,alkane_files,batches)
