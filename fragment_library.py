def fragment_library():

    import numpy as np #this is numpy, allows for data frame and matrix handling
    import pandas

    atom_abundances = dict()
    atom_abundances['c'] = pandas.Series([0.9893,0.0107,0],index=['M0','M1','M2'])
    atom_abundances['h'] = pandas.Series([0.999885,0.000115,0],index=['M0','M1','M2'])
    atom_abundances['n'] = pandas.Series([0.99632,0.00368,0],index=['M0','M1','M2'])
    atom_abundances['o'] = pandas.Series([0.99757,0.00038,0.00205],index=['M0','M1','M2'])
    atom_abundances['si'] = pandas.Series([0.922297,0.046832,0.030872],index=['M0','M1','M2'])

    atom_abundances_df = pandas.DataFrame(atom_abundances)
    atom_abundances_df = pandas.DataFrame.transpose(atom_abundances_df)

    fragment_dict = dict()
    fragment_dict['pyr174'] = {'formula':{'c':6,'h':12,'n':1,'o':3,'si':1},
                               'rt':450,
                               'mz':np.array([174,175,176,177,178]),
                               'areas':dict(),
                               'mid':dict(),
                               'tot_area':0
                              }


    #atom_masses = {'c':12.0107,'h':1.00794,'n':14.0067,'o':15.999,'si':28.0855}
    #atom_quantities = np.array(list(fragment_dict['pyr174']['formula'].values()))
    #atom_masses = np.array(list(atom_masses.values()))

    #fragment_dict['pyr174']['mass_values'] = atom_quantities*atom_masses
    #fragment_dict['pyr174']['tot_mass'] = sum(fragment_dict['pyr174']['mass_values'])
    #fragment_dict['pyr174']['mass_fractions'] = fragment_dict['pyr174']['mass_values']/fragment_dict['pyr174']['tot_mass']
    #mass_fractions_dict = dict()
    #mass_fractions_dict['pyr174'] = pandas.Series(fragment_dict['pyr174']['mass_fractions'],index=['c','h','n','o','si'])
    #mass_fractions_df = pandas.DataFrame(mass_fractions_dict)
    #mass_fractions_df = pandas.DataFrame.transpose(mass_fractions_df)

    #fragment_dict['pyr174']['unlabeled_mid'] = mass_fractions_df.dot(atom_abundances_df)

    return(fragment_dict)
