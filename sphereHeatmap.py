from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import math
import sys

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
		dist=np.linalg.norm( p - pj[0:3] )
		if dist < d0:
			d0 = dist
			cnt = pj[3]
	return cnt

def strip_first_col(fname, delimiter="\t"):
	with open(fname, 'r') as fin:
		next(fin) # skip header line
		for line in fin:
			try:
				yield line.split(delimiter, 1)[1]
			except IndexError:
				continue
	
# read numpy array from file
lines = np.loadtxt(strip_first_col(sys.argv[1]))
pointList = np.zeros(shape=(np.shape(lines)[0],4))
for point in range(len(lines)):
	pointList[point][0] = lines[point][2]
	pointList[point][1] = lines[point][3]
	pointList[point][2] = lines[point][4]
	pointList[point][3] = lines[point][8]
	
#creat pointlist and initialize colormap
precision = 100
u = np.linspace( 0, 2 * np.pi, precision*2)
v = np.linspace( 0, np.pi, precision )
XX = 1 * np.outer( np.cos( u ), np.sin( v ) )
YY = 1 * np.outer( np.sin( u ), np.sin( v ) )
ZZ = 1 * np.outer( np.ones( np.size( u ) ), np.cos( v ) )
WW = np.zeros(shape=(np.shape(XX)))

# populate the colormap and normalize
for i in range( len( XX ) ):
    for j in range( len( XX[0] ) ):
        x = XX[ i, j ]
        y = YY[ i, j ]
        z = ZZ[ i, j ]
        WW[ i, j ] = near(np.array( [x, y, z ] ), pointList)
WW = ( WW - np.amin([int(np.min(pointList[:,3]))]) ) 
WW = WW / np.amax([int(np.max(pointList[:,3]))-int(np.min(pointList[:,3]))])
mycolors = cm.seismic( WW )

# save the figure
fig = plt.figure()
ax = fig.add_subplot( 1, 1, 1, projection='3d')
ax.plot_surface( XX, YY,  ZZ, cstride=1, rstride=1, facecolors=mycolors, shade=False )
#ax.set_aspect(1)
ax.set_xlim([-1.2,1.2])
ax.set_ylim([1.2,-1.2])
ax.set_zlim([-1.2,1.2])
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.axes.zaxis.set_ticks([])
plt.axis('off')
fig.savefig('sphere.png',dpi=300) # save the figure to file
plt.close(fig) # close the figure
