#
# Depth plot example
#
# This is a Stalefish HDF5 data extraction example, using h5py (See
# https://docs.h5py.org/en/latest/quick.html)
#
# It loads data from the .h5 file and make some depth plots for a single,
# selected slice. The plots are saved into a subdirectory ./plots/
#
# Author: Seb James
# Date: October 2020
#

import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
from matplotlib.patches import Polygon
import matplotlib as mpl

# Choose your filename:
filename = '../../data/testimg.h5'

# Set which frame you'd like to get depth plots for (counting from 1):
selected_frame = 3

# Set plotting defaults
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
mpl.rc('font', **fnt)

# This incantation is important to ensure that, in an SVG output file, text is
# written as 'objects that can be edited in inkscape'.
plt.rcParams['svg.fonttype'] = 'none'

# Make sure ./plots exists
if not os.path.isdir('./plots'):
    os.makedirs('./plots')

try:
    with h5py.File(filename, 'r') as f:

        # To see all the 'keys' in the h5 file:
        # print("Keys: {0}".format(list(f.keys())))

        # In my opinion, the way you get the value of a scalar object with h5py
        # is a bit ugly.  You have to convert the hdf5 object into a list, then
        # get the 0th value of the list (which always has one element)
        nf = list(f['/nframes'])[0]
        print('Number of frames: {0}'.format(nf))

        # Create the 'frame tag' (e.g. 'Frame000')
        frame = 'Frame{0:03}'.format(selected_frame)

        # How many boxes are there?
        key = '{0}/nboxes'.format(frame)
        nboxes = list(f[key])[0]

        for j in range(0,nboxes):

            # The means of the sampled boxes are in means or means_autoscaled
            key = '{0}/lmalign/box_depth/box{1}'.format(frame,j)
            box_depth = list(f[key])
            key = '{0}/signal/postproc/boxes/box{1}'.format(frame,j)
            box_signal = list(f[key])

            # Can now plot box_signal vs. box_depth
            plt.clf()
            plt.plot(box_depth, box_signal, color='k', markersize=5, marker='.', linestyle='none')
            plt.margins()
            plt.xlim((0, 1.1))
            plt.ylim((0, 1))
            plt.xlabel('Depth (mm)')
            plt.ylabel('ISH expression')
            plt.tight_layout()

            fname = './plots/Frame{0:03}_box{1:03}.svg'.format(selected_frame, j)
            print('Saving depth plot: {0}'.format(fname))
            plt.savefig(fname, transparent=True)
except:
    print("The file '{0}' couldn't be opened, or it didn't contain a Frame{1:03}. Please open this script and edit the variable 'filename' so that it contains a path to a valid Stalefish .h5 file. You will also need to set 'selected_frame' to one that exists within your file (selected_frame is currently {1}).".format(filename, selected_frame))
