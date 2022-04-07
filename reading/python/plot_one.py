#
# Opens a resampled (and possibly transformed) 2D map from a file that was
# written out by sfview. Makes a plot and shows coordinate axes.
#

import numpy as np
import pylab as pl
import h5py
from PIL import Image

def getZ(fname):
    f = h5py.File(fname,'r')
    x=f['output_map']['twod']['coordinates_resampled'][:,0]
    y=f['output_map']['twod']['coordinates_resampled'][:,1]
    z=f['output_map']['twod']['expression_resampled'][:]
    dims = f['output_map']['twod']['widthheight_resampled'][:]
    X=np.reshape(x, [dims[1],dims[0]])
    Y=np.reshape(y, [dims[1],dims[0]])
    Z=np.reshape(z, [dims[1],dims[0]])
    f.close()
    return Z

def getMask(Z, T=0.2):
    mask = np.where(Z<T)
    return mask

def combi(fnames, thresh=0.2):
    Z = getZ(fnames[0])*0.
    MaskTemp = np.ones(Z.shape)
    M = MaskTemp.copy()
    for i,j in enumerate(fnames):
        z = getZ(j)
        m = MaskTemp.copy()
        m2 = getMask(z, T=thresh)
        m[m2]=np.nan
        Z = np.dstack([Z,z])
        M = np.dstack([M,m])
    return Z[:,:,1:], M[:,:,1:]


D = ['example_data/P9_Mouse_2_Litter_4_RZRB.TF.602630314_annot.h5']

Z, M = combi(D, 0.2)

# PLOTTING
_cmap = 'gray'
ax1 = Z.shape[1]
ax2 = Z.shape[0]
ax = True

F = pl.figure(figsize=(6,6))
f = F.add_subplot(111)
f.imshow(Z[:,:,0]*M[:,:,0],cmap=_cmap)

_img = Image.fromarray(np.uint8(Z[:,:,0]*255))
_img.save ('plot_one.png')

f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

F.tight_layout()

pl.show()
