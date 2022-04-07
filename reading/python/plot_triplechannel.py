#
# Opens a resampled (and possibly transformed) 2D map from a file that was
# written out by sfview. Makes a plot and shows coordinate axes.
#
# This version opens 3 channel colour info, as generated when making a 2D map
# from coloured, atlas slices.
#

import numpy as np
import pylab as pl
import h5py
from PIL import Image

def getZColour(fname):
    f = h5py.File(fname,'r')
    z = f['output_map']['twod']['boxcolours_resampled'][:]
    dims = f['output_map']['twod']['widthheight_resampled'][:]
    Z = np.reshape(z, [dims[1],dims[0],3])
    f.close()
    return Z

# Read data
Z = getZColour ('example_data/602630314_annot.TF.602630314_annot.h5')

ax1 = Z.shape[1]
ax2 = Z.shape[0]
print ('Image size: {0} by {1}'.format (ax1, ax2))

print ('Z0 max: {0}'.format (np.max(Z[:,:,0])))
print ('Z1 max: {0}'.format (np.max(Z[:,:,0])))
print ('Z2 max: {0}'.format (np.max(Z[:,:,0])))

#Z = Z/np.max(Z)

F = pl.figure(figsize=(8,18))
f = F.add_subplot(111)
f.imshow(Z)

_img = Image.fromarray(np.uint8(Z*255), mode="RGB")
_img.save ('triplechannel.png')

f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

F.tight_layout()

pl.savefig('triplechannel_fig.png')
pl.show()
