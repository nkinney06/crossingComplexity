from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import cm
import sys

def cuboid_data(o, size=(1,1,1)):
    X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
    X = np.array(X).astype(float)
    for i in range(3):
        X[:,:,i] *= size[i]
    X += np.array(o)
    return X

def plotCubeAt(positions,sizes=None,colors=None, **kwargs):
    if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
    if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)
    g = []
    for p,s,c in zip(positions,sizes,colors):
        g.append( cuboid_data(p, size=s) )
    return Poly3DCollection(np.concatenate(g),  
                            facecolors=np.repeat(colors,6, axis=0), **kwargs)

def near( p, pntList ):
	(cnt,d0) = (0,1000)
	for pj in pntList:
		dist=np.linalg.norm( p - pj[0][0:3] )
		if dist < d0:
			d0 = dist
			cnt = pj[0][3]
	return cnt
							
def strip_first_col(fname, delimiter="\t"):
	with open(fname, 'r') as fin:
		next(fin) # skip header line
		for line in fin:
			try:
				yield line.split(delimiter, 1)[1]
			except IndexError:
				continue
					
# create numpy arrays from inpus files
lines = np.loadtxt(strip_first_col(sys.argv[1]))
pointList = np.zeros(shape=(np.shape(lines)[0],4))
for point in range(len(lines)):
    pointList[point][0] = lines[point][5]
    pointList[point][1] = lines[point][6]
    pointList[point][2] = lines[point][7]
    pointList[point][3] = lines[point][8]
xFace = np.argwhere(pointList[:,0]==1)
yFace = np.argwhere(pointList[:,1]==1)
zFace = np.argwhere(pointList[:,2]==1)
mxFace = np.argwhere(pointList[:,0]==-1)
myFace = np.argwhere(pointList[:,1]==-1)
mzFace = np.argwhere(pointList[:,2]==-1)
				
# map nearby points
precision = 101
u = np.linspace(-1, 1, precision)
v = np.linspace(-1, 1, precision)
W = list()
sizes = list()
positions = list()
cubeSize = 2./(precision-1)

for i in range( len( u ) ):
    for j in range( len( v ) ):
        x = 1.
        y = u[i]
        z = v[j]
        positions.append([x,y,z])
        sizes.append([cubeSize,cubeSize,cubeSize])
        W.append(near(np.array( [x, y, z] ), pointList[xFace,:]))

for i in range( len( u ) ):
    for j in range( len( v ) ):
        x = u[i]
        y = 1.
        z = v[j]
        positions.append([x,y,z])
        sizes.append([cubeSize,cubeSize,cubeSize])
        W.append(near(np.array( [x, y, z] ), pointList[yFace,:]))
		
for i in range( len( u ) ):
    for j in range( len( v ) ):
        x = u[i]
        y = v[j]
        z = 1.
        positions.append([x,y,z])
        sizes.append([cubeSize,cubeSize,cubeSize])
        W.append(near(np.array( [x, y, z] ), pointList[zFace,:]))

# normalize the colors and make the plot
W = ( W - np.amin([int(np.min(pointList[:,3]))]) )
W = W / np.amax([int(np.max(pointList[:,3]))-int(np.min(pointList[:,3]))])
mycolors = cm.seismic( W )
fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.set_aspect('equal')
pc = plotCubeAt(positions, sizes=sizes, colors=mycolors,edgecolor=None)
ax.add_collection3d(pc)
ax.set_xlim([-1.4,1.4])
ax.set_ylim([1.4,-1.4])
ax.set_zlim([-1.4,1.4])
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.axes.zaxis.set_ticks([])
plt.axis('off')
fig.savefig('cube.png',dpi=300) # save the figure to file
plt.close(fig) # close the figure
