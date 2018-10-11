def expand_polynomial(a,b):
# a: a python array representing the mass isotopomer abundances of atom a
#    the first entry is the abundance of M0, the second M1, the third M2, ...
# b: a python array representing the mass isotopomer abundances of atom b
#    the first entry is the abundance of M0, the second M1, the third M2, ...
# a pandas DataFrame is returned where the rows are the mass isotopomers and the column contains the mass isotopomer abundances as 'prob'
#    This retunred mass isotopomer distribution (MID) is for the molecule comprised of the atoms ab
#    Note this can be used in series to get a multi-atom molecule MID
#        Say you have the molecule abc
#        First run this function with the MIDs for the atoms a and b as inputs; the output is then the MID for the combined ab
#        Run the function again with the inputs being the previous output and the MID for c; the output is then the MID for abc


    #import the module that will allow all other modules to be imported in a single line
    #import ImportModules #the module that contains import statements for the other required modules
    import numpy as np
    import pandas
    import importlib #allows fresh importing of modules
    import pdb

    #redefine the input polynomial python arrays as numpy row matrices
    mid_a = np.matrix(a) #mass isotopomer distribution for a
    mid_b = np.matrix(b) #mass isotopomer distribution for b

    #find the number of recorded relative mass isotopomer abundances for each atom (the number of columns)
    n_a = mid_a.shape[1] #the second entry of the returned tuple
    n_b = mid_b.shape[1] #the second entry of the returned tuple

    #arrange information for atom a into a pandas DataFrame
    mid_a_m = np.zeros([n_a,2]) #initialize the matrix the size of data frame that will be used to hold the information for atom a
    mid_a_df = pandas.DataFrame(mid_a_m,columns=['prob','m_isotopomer']) #initialize the data frame containing the information for atom a
    mid_a_df[['prob']] = np.transpose(mid_a) #enter the relative abundances of each mass isotopomer
    mid_a_df[['m_isotopomer']] = np.transpose(np.matrix(range(0,n_a))) #enter the mz mass isotopomer values corresponding to each mass isotopomer

    #arrange information for atom b into a pandas DataFrame
    mid_b_m = np.zeros([n_b,2])
    mid_b_df = pandas.DataFrame(mid_b_m,columns=['prob','m_isotopomer'])
    mid_b_df[['prob']] = np.transpose(mid_b)
    mid_b_df[['m_isotopomer']] = np.transpose(np.matrix(range(0,n_b)))

    #initialize a DataFrame to contain information on all possible combinations of atoms a and b in the molecule ab
    factored_init_m = np.zeros([n_a*n_b,2])
    factored = pandas.DataFrame(factored_init_m,columns=['prob','m_isotopomer'])

    #fill in the factored_df by expanding the polynomial and keeping track of the mass isotopomer mz value (M0,M1,M2,...) for each term
    #    This is expanding the polynomial and storing each term
    #    The mass isotopomer mz values of each factor are summed for each term to get a new mass isotopomer mz value for each term
    iterator = 0
    for i in range(0,n_a):
        for k in range(0,n_b):
            factored.loc[iterator,'prob'] = mid_a_df.loc[i,'prob']*mid_b_df.loc[k,'prob']
            factored.loc[iterator,'m_isotopomer'] = mid_a_df.loc[i,'m_isotopomer'] + mid_b_df.loc[k,'m_isotopomer']
            iterator = iterator + 1

    #combine the relative mass isotopmer abundances for ab combinations with the same mass isotopomer mz value
    grouped = factored.groupby('m_isotopomer')
    factored = grouped.aggregate(np.sum)

    #shorten factored to remove negligible values
    n_factored = len(factored)
    if n_factored > 22:
        factored = factored.iloc[0:22]

    #return the DataFrame containing the mass isotopomer distrubution of the molecule ab
    return(factored)
