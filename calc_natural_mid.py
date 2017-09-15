def calc_natural_mid(formula):
    #formula: the formula for the fragment whose natural isotopic abundance is being calculated
    #    each atomic symbol must start with a capital letter, if there is a second letter it must be lowercase
    #    a number needs to follow each atomic symbol, even if that number is 1

    #import required modules
    import pdb
    import importlib #allows fresh importing of modules
    import numpy as np #this is numpy
    import pandas
    import expand_polynomial
    import re
    importlib.reload(expand_polynomial) #reload the custom function in case it was changed

    #define the atomic isotopic abundances, there must be an equal number for each atom
    atom_abundances = dict()
    atom_abundances['C'] = pandas.Series([0.9893,0.0107,0],index=['M0','M1','M2'])
    atom_abundances['H'] = pandas.Series([0.999885,0.000115,0],index=['M0','M1','M2'])
    atom_abundances['N'] = pandas.Series([0.99632,0.00368,0],index=['M0','M1','M2'])
    atom_abundances['O'] = pandas.Series([0.99757,0.00038,0.00205],index=['M0','M1','M2'])
    atom_abundances['Si'] = pandas.Series([0.922297,0.046832,0.030872],index=['M0','M1','M2'])

    #determine how many relative abundances are being considered
    n_rel_abuns = len(atom_abundances['C'])

    #make a data frame out of the atom_abundances library
    #atom_abundances_df = pandas.DataFrame(atom_abundances)
    #atom_abundances_df = pandas.DataFrame.transpose(atom_abundances_df)

    #break the fragment formula up into its atomic symbol and letter components
    broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', formula))
    #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
    #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
    #        all components are strings

    #get separate arrays containing the atomic symbols and corresponding formula numbers
    odd_indices = np.array(range(1,len(broken_formula),2)) #all of the odd indices of the array broken_formula
    even_indices = np.array(range(0,len(broken_formula),2)) #all of the even indices of the array broken_formula
    formula_numbers = broken_formula[odd_indices].astype(np.int) #the atomic numbers array
    formula_atoms = broken_formula[even_indices] #the atomic symbols array

    #determine the number of atoms in the fragment
    n_row = np.sum(formula_numbers)

    #create a matrix containing all of the the atom mids
    #    each row corresponds to an atom in the fragment
    #    if an atom appears more than once in the fragment, the matrix will contain idenical rows for each appearance
    #    the isotope mass varies across the columns
    #    example:
    #        row corresponds to:
    #                     carbon MID:   M0 M1 M2 M3...
    #                     carbon MID:   M0 M1 M2 M3...
    #                                        .
    #                                        .
    #                     hydrogen MID: M0 M1 M2 M3...
    #                     hydrogen MID: M0 M1 M2 M3...
    #                                        .
    #                                        .
    atom_mids = np.zeros([n_row,n_rel_abuns]) #initialize the matrix
    atom_index = 0 #initialize the atom counter, this tracks the identity of the atom
                   #    if you go from the first to the second C, this will not changed
                   #    if you go from the last C to the first H this will increase by one
    for i in range(0,n_row): #iterate through the rows of the atom mid matrix
        new_atom_id_index = sum(formula_numbers[0:atom_index+1]) #'+1' is necessary because the end of the range is not included in python
        #corresponds to the next row that the identity of the atom will change
        if i < new_atom_id_index: #if the atom identity has not changed
            current_atom = formula_atoms[atom_index]
            atom_mids[i,] = atom_abundances[current_atom]
        if i >= new_atom_id_index: #if the atom identity has changed
            current_atom = formula_atoms[atom_index+1]
            atom_mids[i,] = atom_abundances[current_atom]
            atom_index = atom_index + 1 #only iterate the atom_index if the atom identity has changed

    expanded_placeholder = [1]
    for i in range(0,n_row):
        expanded = expand_polynomial.expand_polynomial(expanded_placeholder,atom_mids[i,])
        expanded_placeholder = np.array(expanded['prob'])

    natural_mid = np.array(expanded.iloc[:,0])
    return(natural_mid)
