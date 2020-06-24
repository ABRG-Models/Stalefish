# Utility developed from read_example.py
#
# Load data and make some plots


import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
from matplotlib.patches import Polygon
import matplotlib as mpl
import sys

# Get a file name from the cmd line
if len(sys.argv) < 2:
    print ('Usage: {0} path/to/experiment.h5'.format (sys.argv[0]))
    exit(0)

filename = sys.argv[1]

with h5py.File (filename, 'r') as f:
    #print("Keys: {0}".format(list(f.keys())))

    # It's a bit ugly, the way you get the value of the nframes
    # object; you have to convert the hdf5 object into a list, then
    # get the 0th value of the list (which always has one element)
    nf = list(f['/nframes'])[0]
    print ('Number of frames: {0}'.format(nf))

    # Empty lists to hold data. These'll be lists of lists
    x = []
    y = []
    means = []
    thicks = []

    for i in range(1,nf+1):

        # Create the 'frame tag' (e.g. 'Frame000')
        frame = 'Frame{0:03}'.format(i)
        #print ('frame: {0}'.format(frame))

        # The means of the sampled boxes are in means or means_autoscaled
        key = '{0}/means_autoscaled'.format(frame)
        means_ = list(f[key])

        # The position on the 3D x-axis is in class/layer_x:
        key = '{0}/class/layer_x'.format(frame)

        # This creates a list with as many elements as there are
        # means, all with the same value (class/layer_x)
        x_ = list(f[key]) * len(means_)

        # For y values, use sbox_linear_distance
        key = '{0}/sbox_linear_distance'.format(frame)
        y_ = list(f[key])

        # Use the per-frame slice thickness when plotting rectangular polygons
        key = '{0}/class/thickness'.format(frame)
        thicks_ = list(f[key]) * len(means_)

        # Now have four lists, x_, y_, means_ and thicks_. These can
        # be appended onto the 'lists of lists' for all slices.
        x.append(x_)
        y.append(y_)
        means.append(means_)
        thicks.append(thicks_)

    # We now have (as lists of lists) x, y, means and the thicknesses
    # (in dimension x) of all the boxes. A regular scatter graph gives
    # awful results. Instead, plot each x/y as a rectangle with colour
    # given by means.
    plt.title('Flattened heat plot')
    for ii in range(0,len(x)): # len(x) is the number of frames.
        # Within each frame, make polygons
        x_l = np.array(x[ii])
        x_r = np.add (np.array(x[ii]), np.array(thicks[ii]))
        y_ar = np.array(y[ii]) # convert the list y[ii] into a np array
        y_t = y_ar + (shift(y_ar, -1, cval=y_ar[-1]) - y_ar) / 2.0 # half way between y and the next y
        y_b = y_ar - (y_ar - shift(y_ar, 1, cval=y_ar[0])) / 2. # half way between y and the previous y

        for jj in range(0,len(x[ii])): # len(x[ii]) is the number of bins in a frame
            # Boxes have vertices: (x_l,y_b) -> (x_r,y_b) -> (x_r,y_t) -> (x_l,y_t)
            v = np.array([ [ x_l[jj], y_b[jj] ], # v stands for 'vertices'
                           [ x_r[jj],  y_b[jj] ],
                           [ x_r[jj], y_t[jj] ],
                           [ x_l[jj], y_t[jj] ] ])
            poly = Polygon (v, facecolor=mpl.cm.jet(means[ii][jj]), edgecolor='None')
            plt.gca().add_patch(poly)


    plt.xlabel('x')
    plt.ylabel('y')
    plt.autoscale()
    plt.show()
