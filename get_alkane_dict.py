def get_alkane_dict():
    import numpy as np
    import pdb

    alkane_names = np.array([])
    for nc in alkane_nc:
        to_append = np.array_str(nc) + '_carbon_alkane'
        alkane_names = np.append(alkane_names,to_append )
        #should look something like the following, but in the reverse order...
        alkane_names = np.array(['08_carbon_alkane','09_carbon_alkane','10_carbon_alkane','11_carbon_alkane','12_carbon_alkane',])


    #will then be able to make the dictionary as in the following using a loop...
    #    the first group is the parent mass, there is also a single tail peak
    #    the second is the parent mass less 30, there is also the parent mass less 29 which is larger
    alkane_dict = dict()
    alkane_dict['8_carbon_alkane']['peak_profile'] =
    alkane_dict['9_carbon_alkane']['peak_profile'] =
    alkane_dict['10_carbon_alkane']['peak_profile'] =
    alkane_dict['11_carbon_alkane']['peak_profile'] =
    alkane_dict['12_carbon_alkane']['peak_profile'] =
    alkane_dict['13_carbon_alkane']['peak_profile'] = 184 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['14_carbon_alkane']['peak_profile'] = 198 100000 0.7 0.3 168 90000 0.3 0.6
    alkane_dict['15_carbon_alkane']['peak_profile'] = 212 100000 0.7 0.3 182 90000 0.3 0.6
    alkane_dict['16_carbon_alkane']['peak_profile'] = 226 100000 0.7 0.3 196 90000 0.3 0.6
    alkane_dict['17_carbon_alkane']['peak_profile'] = 240 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['18_carbon_alkane']['peak_profile'] = 254 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['19_carbon_alkane']['peak_profile'] = 268 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['20_carbon_alkane']['peak_profile'] = 282 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['21_carbon_alkane']['peak_profile'] = 296 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['22_carbon_alkane']['peak_profile'] = 310 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['23_carbon_alkane']['peak_profile'] = 324 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['24_carbon_alkane']['peak_profile'] = 338 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['25_carbon_alkane']['peak_profile'] = 352 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['26_carbon_alkane']['peak_profile'] = 366 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['27_carbon_alkane']['peak_profile'] = 380 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['28_carbon_alkane']['peak_profile'] = 394 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['29_carbon_alkane']['peak_profile'] = 408 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['30_carbon_alkane']['peak_profile'] = 422 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['31_carbon_alkane']['peak_profile'] = 436 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['32_carbon_alkane']['peak_profile'] = 450 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['33_carbon_alkane']['peak_profile'] = 464 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['34_carbon_alkane']['peak_profile'] = 478 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['35_carbon_alkane']['peak_profile'] = 492 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['36_carbon_alkane']['peak_profile'] = 506 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['37_carbon_alkane']['peak_profile'] = 520 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['38_carbon_alkane']['peak_profile'] = 534 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['39_carbon_alkane']['peak_profile'] = 548 100000 0.7 0.3 154 90000 0.3 0.6
    alkane_dict['40_carbon_alkane']['peak_profile'] = 562 100000 0.7 0.3 154 90000 0.3 0.6
