import numpy as np
import pylab as pl
import h5py
import cv2 as cv
import sys

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

def getP(Z, T=0.0008, k1=43, k2=93):

    # Flip on both axes
    Z = np.flipud(Z)
    Z = np.fliplr(Z)

    # Set outside pixels to mean inside (min. edge effects on blurring)
    mask = np.where(Z<T)
    avg = np.mean(Z[np.where(Z>=T)])
    Z[mask]=avg

    # Difference of Gaussians blurring
    D = cv.GaussianBlur(Z,(k1,k1),0)-cv.GaussianBlur(Z,(k2,k2),0)
    P = np.ma.masked_array(D, mask=(Z==avg))

    # Set outside pixels to max (so they image white)
    Z[mask] = np.max(Z)
    return Z, P


if (len(sys.argv)<3):
    k1 = 37     # width of center Gaussian DoG
    k2 = 93     # width of surround Gaussian DoG
else:
    k1 = int(sys.argv[1])
    k2 = int(sys.argv[2])


Z1 = getZ('65_7E_id2_L23_3gl.h5')
Z2 = getZ('65_8A_id2_L23_3gl.h5')
Z3 = getZ('66_6B_id2_L23_3gl.h5')

Z1Wht, P1 = getP(Z1,k1=k1,k2=k2)
Z2Wht, P2 = getP(Z2,k1=k1,k2=k2)
Z3Wht, P3 = getP(Z3,k1=k1,k2=k2)


# Plotting params
fs=8
cmp = 'gray'
ax1 = Z1.shape[1]
ax2 = Z1.shape[0]

re = '#FF0000'
gr = '#00FF00'
bl = '#0000FF'

# PLOTTING
F = pl.figure(figsize=(12,6))

f = F.add_subplot(131)
f.imshow(Z1Wht,cmap=cmp)
f.axis('equal')
f.axis('off')

f = F.add_subplot(132)
f.imshow(Z2Wht,cmap=cmp)
f.axis('equal')
f.axis('off')

f = F.add_subplot(133)
f.imshow(Z3Wht,cmap=cmp)
f.axis('equal')
f.axis('off')

F = pl.figure(figsize=(12,6))

f = F.add_subplot(141)
f.imshow(Z1Wht,cmap=cmp)
f.contour(P1,[0],colors=re)
f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

f = F.add_subplot(142)
f.imshow(Z2Wht, cmap=cmp)
f.contour(P2,[0],colors=gr)
f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

f = F.add_subplot(143)
f.imshow(Z3Wht, cmap=cmp)
f.contour(P3,[0],colors=bl)
f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

# OVERLAY
f = F.add_subplot(144)
f.imshow(np.ones(Z1.shape)*np.nan, cmap=cmp)
f.contourf(P1,[0,1],colors=re,alpha=0.5)
f.contourf(P2,[0,1],colors=gr,alpha=0.5)
f.contourf(P3,[0,1],colors=bl,alpha=0.5)
f.plot([ax1/2,ax1/2],[0,ax2],'-',color=(0,0,0))
f.plot([0,ax1],[ax2/2,ax2/2],'-',color=(0,0,0))
f.axis('equal')
f.axis('off')

pl.show()
