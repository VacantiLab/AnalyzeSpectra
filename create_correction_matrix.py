def create_correction_matrix(formula):

    import importlib
    import pdb
    import copy
    import numpy as np
    import pandas
    import re
    import calc_natural_mid
    importlib.reload(calc_natural_mid)

    #break the formula up so atoms and quantities are consecutive entries in a numpy array
    broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', formula))
    #    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
    #    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
    #        all components are strings
    n_formula_entries = len(broken_formula)

    #the number of rows of the correction matrix is equal to the quantity of the atom being corrected for
    #    this should really be the number of that atom in the metabolite, will shorten run time and give more accurate results
    #    metabolite atoms need to be recorded in the library
    atom_index = np.where(broken_formula=='C')[0][0]
    atom_quantity_index = atom_index+1
    atom_quantity = broken_formula[atom_quantity_index]
    atom_quantity = atom_quantity.astype(np.int)
    n_rows = atom_quantity

    #temporary for citrate!!!
    #need to replace this with the number of metabolite atoms in the formula
    atom_quantity = 6

    #add the "heavy atom to the end of the broken formula array", initially its quantity is 0
    broken_formula = np.append(broken_formula,np.array(['Hv','0']))
    n_formula_entries_heavy = len(broken_formula)

    #replace each atom of interest with a heavy atom and get the natural mid of the result
    #    these mids fill the rows of the correction matrix
    broken_formula_correct = copy.copy(broken_formula) #initialize the array to carry the formula with a heavy atom
    correction_matrix_dict = dict() #initialize a dictionary to hold the rows of the correction matrix
    n_theoretical_mid_entries = atom_quantity + 4 #the theoretical length is all possible atoms to label plus 4 (after this the relative abundances are assumed negligible)
    for i in range(0,atom_quantity+1):
        #subtract an atom of interest from the formula
        broken_formula_correct[atom_quantity_index] = broken_formula[atom_quantity_index].astype(np.int) - i
        broken_formula_correct[atom_quantity_index] = broken_formula_correct[atom_quantity_index].astype(np.str)

        #replace that atom with a heavy atom
        broken_formula_correct[n_formula_entries+1] = broken_formula[n_formula_entries+1].astype(np.int) + i
        broken_formula_correct[n_formula_entries+1] = broken_formula_correct[n_formula_entries+1].astype(np.str)

        #update the string version of the formula from the array version
        new_formula = ''
        for j in range(0,n_formula_entries_heavy):
            new_formula = new_formula + broken_formula_correct[j]

        #get the mid due to natural abundances of the updated formula (with one or more heavy atoms)
        correction_matrix_dict[i] = calc_natural_mid.calc_natural_mid(new_formula)

        #shorten each theoretical MID with given quantities of heavy atoms to the specified length of the theoretical MIDs
        correction_matrix_dict[i] = calc_natural_mid.calc_natural_mid(new_formula)[0:n_theoretical_mid_entries]
        correction_row_normalization_factor = sum(correction_matrix_dict[i])
        correction_matrix_dict[i] = correction_matrix_dict[i]/correction_row_normalization_factor

    #make the correction matrix dictionary into a matrix
    CM = pandas.DataFrame(correction_matrix_dict)
    CM = pandas.DataFrame.transpose(CM)
    CM = pandas.DataFrame.as_matrix(CM)

    #find the right inverse (pseudo-inverse in numpy jargon) of the correction matrix
    CM_i = np.linalg.pinv(CM)

    return(CM_i)
