#!/usr/bin/env python

file = None
#file = 'output.pdf'

# import libraries
import matplotlib
matplotlib.use('macosx' if file is None else 'PDF')

#import pyfits
import astropy.io.fits as pyfits
import numpy as np
from pylab import *

# load field evolution data
hdu = pyfits.open('output.fits'); q = hdu[0].data

ny,nx = q.shape; print("Input data dimensions is (%i,%i) pixels" % (ny,nx))

# antialias using median filter
#osy = ny/1024; osx = nx/1536
osy = 1; osx = 1

if (osy != 1 or osx != 1):
	from scipy.signal import medfilt
	aakernel = [(osy & (~1)) + 1, (osx & (~1)) + 1]
	q = medfilt(q, aakernel)[::osy,::osx]; ny,nx = q.shape
	print("Antialiased with median kernel %s; output size is (%i,%i) pixels" % (aakernel, ny,nx))

# data extent and meshes
x0 = hdu[0].header['CRVAL1']
y0 = hdu[0].header['CRVAL2']
dx = hdu[0].header['CDELT1']
dy = hdu[0].header['CDELT2']

x = x0 + arange(nx)*dx*osx
y = y0 + arange(ny)*dy*osy

X, Y = meshgrid(x, y)

hdu.close()

# create the figure
figure(figsize=(10,8), frameon=False)
gradient = ["plum", "blue", "palegreen", "red", "white"]
c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", gradient)
imshow(q, extent=(x[0],x[-1],y[0],y[-1]), origin='lower', aspect=(osx*dx)/(osy*dy), cmap=c, interpolation='none')
colorbar(shrink=1.0).set_label('$\\phi$')

# overlay field isolines
overlay = contour(X, Y, q); clabel(overlay, inline=1, fontsize=10)

# overlay field profiles
for t in range(50,ny,100):
	plot(x, y[t]+100.0*q[t,:], 'k-', linewidth=0.5, alpha=0.3)

xlim([x[0],x[-1]]); #xticks(arange(5)-11.0)
ylim([y[0],y[-1]]); #yticks(arange(5)-11.0)

xlabel('$x$')
ylabel('$t$')

if not(file is None): plt.savefig(file, bbox_inches='tight', pad_inches=0.02, transparent=True)
show()
