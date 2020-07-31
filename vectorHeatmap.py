from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import sys

def strip_first_col(fname, delimiter="\t"):
	with open(fname, 'r') as fin:
		next(fin) # skip header line
		for line in fin:
			try:
				yield line.split(delimiter, 1)[1]
			except IndexError:
				continue
	
lines = np.loadtxt(strip_first_col(sys.argv[1])) # read data from file
pointList = np.zeros(shape=(np.shape(lines)[0],4))
W = list()

for point in range(len(lines)):
	pointList[point][0] = lines[point][5]
	pointList[point][1] = lines[point][6]
	pointList[point][2] = lines[point][7]
	pointList[point][3] = lines[point][8]
	W.append(lines[point][8])

W = ( W - np.amin([int(np.min(pointList[:,3]))]) )
W = W / np.amax([int(np.max(pointList[:,3]))-int(np.min(pointList[:,3]))])
mycolors = cm.seismic( W )
fig = plt.figure()
ax = fig.gca(projection='3d')

for i in range(len(pointList)):
	ax.quiver(0, 0, 0, pointList[i][0], pointList[i][1], pointList[i][2], colors=mycolors[i], length=10., normalize=True, arrow_length_ratio = .1)
#ax.set_aspect(1)
ax.set_xlim([-11.1,11.1])
ax.set_ylim([11.1,-11.1])
ax.set_zlim([-11.1,11.1])
ax.axes.xaxis.set_ticks([])
ax.axes.yaxis.set_ticks([])
ax.axes.zaxis.set_ticks([])
plt.axis('off')
fig.savefig('vector.png',dpi=300) # save the figure to file
plt.close(fig) # close the figure
