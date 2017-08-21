def remove_repeats(mz,ic,si):
    import numpy as np
    import pdb
    n_scns = len(si)
    si = np.array(si)
    n_mz = len(mz)
    mz = list(mz[0:n_mz])
    ic = list(ic[0:n_mz])
    mz_r = [round(element,0) for element in mz]
    mz_dif = np.subtract(mz,mz_r)
    mz_dif = [abs(element) for element in mz_dif]
    loop_indices = range(0,n_mz-1)
    del_dict = {}
    j=0
    k=0
    for i in loop_indices:
        if mz_r[i] == mz_r[i+1]:
            pair_list = [mz_dif[i],mz_dif[i+1]]
            max_pair_list_index = pair_list.index(max(pair_list))
            del_dict[j] = i + max_pair_list_index
            j = j+1
    to_remove = list(dict.values(del_dict))
    to_remove_desc = sorted(to_remove,reverse=True)
    n_to_rem = len(to_remove_desc)
    for i in range(0,n_to_rem):
        del mz_r[to_remove_desc[i]]
        del ic[to_remove_desc[i]]
        si = np.where(si>to_remove_desc[i],si-1,si)
    si = list(si)
    return mz_r, ic, si
