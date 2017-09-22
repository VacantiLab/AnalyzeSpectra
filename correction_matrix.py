import importlib
import pdb
import copy
import numpy as np
import re
import calc_natural_mid
importlib.reload(calc_natural_mid)

formula = 'C6H12N1O3Si1'
mid_u = np.array([ 0.8444115,0.10784203,0.04012439,0.00762208,0.0])
n_uncor = len(mid_u)
n_col = n_uncor

#break the formula up so atoms and quantities are consecutive entries in a numpy array
broken_formula = np.array(re.findall('[A-Z][a-z]?|[0-9]+', formula))
#    '[A-Z][a-z]?|[0-9]+': A capital letter [A-Z], followed by a lowercase letter [a-z], which is optional '?', or '|' a number '[0-9]', and possibly more numbers '+'
#    example: this command will take formula = C6H12N1O3Si1 and return broken_formula = array(['C','6','H','12','N','1','O','3','Si','1'])
#        all components are strings
n_formula_entries = len(broken_formula)

#the number of rows of the correction matrix is equal to the quantity of the atom being corrected for
atom_index = np.where(broken_formula=='C')[0][0]
atom_quantity_index = atom_index+1
atom_quantity = broken_formula[atom_quantity_index]
atom_quantity = atom_quantity.astype(np.int)
n_rows = atom_quantity

#add the "heavy atom to the end of the broken formula array", initially its quantity is 0
broken_formula = np.append(broken_formula,np.array(['Hv','0']))
n_formula_entries_heavy = len(broken_formula)

#replace each atom of interest with a heavy atom and get the natural mid of the result
#    these mids fill the rows of the correction matrix
broken_formula_correct = copy.copy(broken_formula) #initialize the array to carry the formula with a heavy atom
correction_matrix_dict = dict() #initialize a dictionary to hold the rows of the correction matrix
for i in range(0,atom_quantity+1):
    #subtract an atom of interest from the formula
    broken_formula_correct[atom_quantity_index] = broken_formula[atom_quantity_index].astype(np.int) - i
    broken_formula_correct[atom_quantity_index] = broken_formula_correct[atom_quantity_index].astype(np.str)

    #replace that atom with a heavy atom
    broken_formula_correct[n_formula_entries+1] = broken_formula[n_formula_entries+1].astype(np.int) + i
    broken_formula_correct[n_formula_entries+1] = broken_formula_correct[n_formula_entries+1].astype(np.str)

    #update the string version of the formula from the array version
    print(broken_formula_correct)
    new_formula = ''
    for j in range(0,n_formula_entries_heavy):
        new_formula = new_formula + broken_formula_correct[j]

    #get the mid due to natural abundances of the updated formula (with one or more heavy atoms)
    correction_matrix_dict[i] = calc_natural_mid.calc_natural_mid(new_formula)
