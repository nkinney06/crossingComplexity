from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import sys
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

def random_point( r=1 ):
    ct = 2*np.random.rand() - 1
    st = np.sqrt( 1 - ct**2 )
    phi = 2* np.pi *  np.random.rand()
    x = r * st * np.cos( phi)
    y = r * st * np.sin( phi)
    z = r * ct
    return np.array( [x, y, z ] )

def near( p, pntList ):
	(cnt,d0) = (0,1000)
	for pj in pntList:
		dist=np.linalg.norm( p - pj[0:2] )
		if dist < d0:
			d0 = dist
			cnt = pj[2]
	return cnt

def strip_first_col(fname, delimiter="\t"):
	with open(fname, 'r') as fin:
		next(fin) # skip header line
		for line in fin:
			try:
				yield line.split(delimiter, 1)[1]
			except IndexError:
				continue

# create the numpy array of points from the file
lines = np.loadtxt(strip_first_col(sys.argv[1]))
pointList = np.zeros(shape=(np.shape(lines)[0],3))
for point in range(len(lines)):
    pointList[point][0] = lines[point][0]
    pointList[point][1] = lines[point][1]
    pointList[point][2] = lines[point][8]

# generate a grid of points and initialize colormap
precision = 100
u = np.linspace( -np.pi, np.pi, precision*2 )
v = np.linspace(      0, np.pi, precision   )
XX, YY = np.meshgrid(u, v, sparse=False, indexing='ij')
WW = np.zeros(shape=(np.shape(XX))) #this will be the colormap

# populate the heatmap with nearest crossing number
for i in range( len( XX ) ):
    for j in range( len( XX[0] ) ):
        x = XX[ i, j ]
        y = YY[ i, j ]
        WW[ i, j ] = near(np.array( [x, y] ), pointList)
	
# generate the plot
fg = plt.figure(figsize=(8,4))
ax = fg.add_subplot( 1, 1, 1 )
levels = MaxNLocator(nbins=64).tick_values(int(np.min(pointList[:,2])), int(np.max(pointList[:,2])))
cmap = plt.get_cmap('seismic')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
im = ax.pcolormesh( XX, YY, WW, cmap=cmap, norm=norm)
plt.ylim(np.pi, 0)
plt.xlim(np.pi, -np.pi)
ax_divider = make_axes_locatable(ax)
cax = ax_divider.append_axes('top', size = '5%', pad = '11%')
fg.colorbar(im,cax=cax,orientation="horizontal")
fg.savefig('flat.png',dpi=300) # save the figure to file
plt.close(fg)


