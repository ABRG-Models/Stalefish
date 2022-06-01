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
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib import cm
from scipy.ndimage.interpolation import shift
from matplotlib.patches import Polygon
import matplotlib as mpl

# Colours
gray70=(0.7019607843,0.7019607843,0.7019607843)
deepskyblue2=(0,0.6980392157,0.9333333333)

# Choose your filename:
filename = 'example_data/Fig_depth.h5'

# Set which frame you'd like to get depth plots for (counting from 1):
selected_frame = 1

# Set plotting defaults
fs = 18
fnt = {'family' : 'Arial',
       'weight' : 'regular',
       'size'   : fs}
mpl.rc('font', **fnt)

# This incantation is important to ensure that, in an SVG output file, text is
# written as 'objects that can be edited in inkscape'.
plt.rcParams['svg.fonttype'] = 'none'
# We're going to open quite a few figures...
plt.rcParams['figure.max_open_warning'] = 120

# Make sure ./svg dir exists
if not os.path.isdir('./svg'):
    os.makedirs('./svg')

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

        # Box centers give linear distances
        key = '{0}/autoalign/sbox_centers'.format(frame)
        sbox_centers = list(f[key])
        first=1
        Y=np.zeros(len(sbox_centers))
        ypos=0
        ii=int(0)
        for sbc in sbox_centers:
            if first:
                first=0
            else:
                dy = sbc[1]-sbc_last[1]
                dz = sbc[2]-sbc_last[2]
                dl = np.sqrt (dy*dy + dz*dz)
                ypos += dl
                Y[ii] = ypos
                #print ('Box {0} is at linear distance position {1}'.format(ii, ypos))
                ii = ii+1
            sbc_last=sbc

        print ('Generating SVG plots...')
        do_dplots=1
        for j in range(0,nboxes): # or 0,nboxes

            # The means of the sampled boxes are in means or means_autoscaled
            key = '{0}/autoalign/box_depth/box{1}'.format(frame,j)
            box_depth = list(f[key])
            key = '{0}/signal/postproc/boxes/box{1}'.format(frame,j)
            box_signal = list(f[key])

            # Can now plot box_signal vs. box_depth
            if do_dplots:
                plt.clf()
                fig = plt.figure (figsize=(6,3))
                axs = fig.add_subplot(111)
                axs.plot(box_depth, box_signal, color=gray70, markersize=2, marker='.', linestyle='none')

            # Plot histogrammed. Histogram first with signal as weights
            hist, bedges = np.histogram (box_depth, bins=100, weights=box_signal)
            # And second with no weights to get the counts
            hcounts, bedges = np.histogram (box_depth, bins=100)
            # then compute mean
            mean_signal = hist/hcounts
            # Plot on lower subplot
            if do_dplots:
                axs.plot(bedges[1:], mean_signal, color=deepskyblue2, linewidth=2, marker='None', linestyle='-')
            # Also add histogrammed data onto a container for plotting
            X=bedges[1:]
            if j>0:
                msigs = np.vstack((msigs,mean_signal))
            else:
                msigs = mean_signal

            if do_dplots:
                axs.margins()
                axs.set_xlim((0, 0.6))
                axs.set_ylim((0.2, 0.8))
                axs.set(ylabel='ISH')
                axs.set(xlabel='Depth (mm)', ylabel='Id2 ({0:.1f} mm)'.format(Y[j]))

                plt.tight_layout()

                fname = './svg/Frame{0:03}_box{1:03}.svg'.format(selected_frame, j)
                #print('Saving depth plot: {0}'.format(fname))
                plt.savefig(fname, transparent=True)

        X, Y = np.meshgrid(X, Y)
        plt.clf()
        fig = plt.figure()
        ax = fig.gca()
        flatmap = ax.imshow (msigs.T, cmap=cm.viridis, extent=[np.min(Y),np.max(Y),np.max(X),np.min(X)], aspect='auto')
        ax.set(xlabel='Lin dist (mm)', ylabel='Depth (mm)')

        fname = './svg/Frame{0:03}_surf.svg'.format(selected_frame)
        plt.savefig(fname, transparent=True)

        print ('Done. See ./svg/Frame*.svg for the output of this script.')

except Exception as e:
    print("Exception: {0}".format(str(e)))
