def calc_natural_mid(formula):

    import pdb
    import importlib #allows fresh importing of modules
    import numpy as np #this is numpy
    import pandas
    import expand_polynomial
    import re
    importlib.reload(expand_polynomial) #reload the custom function in case it was changed

    atom_abundances = dict()
    atom_abundances['C'] = pandas.Series([0.9893,0.0107,0],index=['M0','M1','M2'])
    atom_abundances['H'] = pandas.Series([0.999885,0.000115,0],index=['M0','M1','M2'])
    atom_abundances['N'] = pandas.Series([0.99632,0.00368,0],index=['M0','M1','M2'])
    atom_abundances['O'] = pandas.Series([0.99757,0.00038,0.00205],index=['M0','M1','M2'])
    atom_abundances['Si'] = pandas.Series([0.922297,0.046832,0.030872],index=['M0','M1','M2'])

    n_rel_abuns = len(atom_abundances['C'])

    atom_abundances_df = pandas.DataFrame(atom_abundances)
    atom_abundances_df = pandas.DataFrame.transpose(atom_abundances_df)

    broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', formula))
    # A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'

    odd_indices = np.array(range(1,len(broken_formula),2))
    even_indices = np.array(range(0,len(broken_formula),2))
    formula_numbers = broken_formula[odd_indices].astype(np.int)
    formula_atoms = broken_formula[even_indices]

    n_row = np.sum(formula_numbers)

    atom_mids = np.zeros([n_row,n_rel_abuns])

    atom_index = 0
    for i in range(0,n_row):
        if i < sum(formula_numbers[0:atom_index+1]):
            atom_mids[i,] = atom_abundances[formula_atoms[atom_index]]
        if i >= sum(formula_numbers[0:atom_index+1]):
            atom_mids[i,] = atom_abundances[formula_atoms[atom_index+1]]
            atom_index = atom_index + 1

    expanded_placeholder = [1]
    for i in range(0,n_row):
        expanded = expand_polynomial.expand_polynomial(expanded_placeholder,atom_mids[i,])
        expanded_placeholder = np.array(expanded['prob'])

    natural_mid = np.array(expanded.iloc[:,0])
    return(natural_mid)
