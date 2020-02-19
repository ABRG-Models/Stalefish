# Example following https://docs.h5py.org/en/latest/quick.html
import h5py
import numpy as np

filename = '../../ucr/V_Id2_1.h5'

with h5py.File (filename, 'r') as f:
    print("Keys: {0}".format(list(f.keys())))

    for group_key in list(f.keys()):
        print ('The data in object {0} is:'.format (group_key))

        # Get the data
        data = list(f[group_key])
        for d in data:
            print (' {0}'.format(d))
