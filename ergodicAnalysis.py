import pandas as pd
import pylab as pl
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

########################################################################################
## example of triple product analysis
########################################################################################
directory = sys.argv[1]
datafile = [directory + "/" + x + "/" + y for x in os.listdir(directory) for y in os.listdir(directory + "/" + x) if "curve" in y]

data = dict()
for file in datafile:
	with open(file) as file_in:
		lines = []
		for line in file_in:
			lines.append([float(x) for x in line.strip().split(" ")])
		data[file] = np.array(lines)
		
tripleSums = list()
for j in range(len(datafile)):
	tripleSum = 0
	for i in range(0,data[datafile[j]].shape[0]-3):
		a = np.subtract(data[datafile[j]][i+1],data[datafile[j]][i+0]) # vector A
		b = np.subtract(data[datafile[j]][i+2],data[datafile[j]][i+1]) # vector B
		c = np.subtract(data[datafile[j]][i+3],data[datafile[j]][i+2]) # vector C
		tripleSum += np.dot(a,np.cross(b,c))
	tripleSums.append(tripleSum)

mydata = { 'torsion':tripleSums }	   
df = pd.DataFrame(mydata)
styles = ['b-']
linewidths = [3]

fig, ax = plt.subplots(figsize=(3,3))
for col, style, lw in zip(df.columns, styles, linewidths):
    df[col].plot.kde(style=style, lw=lw, ax=ax)
ax.set_xlabel(u'\u2211[a\u2022(bxc)]',fontsize=18)
ax.set_ylabel("density",fontsize=18)
ax.set_xlim(min(tripleSums),max(tripleSums))
plt.title("Torsion",fontsize=18)
plt.tight_layout()
fig.savefig('torsion.png',dpi=300)

########################################################################################
## bond alignment example
########################################################################################
alignmentsList = list()
for j in range(len(datafile)):
	alignments = [0]*3
	for i in range(0,data[datafile[j]].shape[0]-1):
		a = np.subtract(data[datafile[j]][i+1],data[datafile[j]][i+0]) # vector A
		alignments[np.argmax(abs(a))] += 1
	alignmentsList.append(alignments)

xx = [ float(x[0])/(x[0]+x[1]+x[2]) for x in alignmentsList ]
yy = [ float(x[1])/(x[0]+x[1]+x[2]) for x in alignmentsList ]
zz = [ float(x[2])/(x[0]+x[1]+x[2]) for x in alignmentsList ]

mydata = [xx,yy,zz]
fig7, ax7 = plt.subplots(figsize=(3,3))
ax7.boxplot(mydata,labels=['x','y','z'],showfliers=False)
ax7.set_ylabel("% total bonds",fontsize=18)
ax7.set_xticklabels(['x','y','z'])
ax7.set_ylim(.2,.5)
ax7.set_xlabel("axis",fontsize=18)
plt.title("Bond Alignments",fontsize=18)
plt.tight_layout()
fig7.savefig('alignment.png',dpi=300)

########################################################################################
## end occupanct example
########################################################################################
boxSize = np.max(data[datafile[0]]) + 1
counts = np.zeros(int(boxSize*boxSize*boxSize))
convertLinear = np.array([1,boxSize,boxSize*boxSize])
for j in range(len(datafile)):
	counts[int(np.dot(data[datafile[j]][0],convertLinear))] += 1
	counts[int(np.dot(data[datafile[j]][-1],convertLinear))] += 1
counts = counts/(2*len(datafile))
fig7, ax7 = plt.subplots(figsize=(3,3))
ax7.scatter(range(1,len(counts)+1), counts,s=3);
ax7.set_xticks([1,int(float(boxSize**3)*.25),int(float(boxSize**3)*.5),int(float(boxSize**3)*.75),boxSize**3])
ax7.set_yticks([0,.5/boxSize,1./boxSize,1.5/boxSize,2./boxSize,2.5/boxSize])
ax7.set_yticklabels(['0','0.5/n','1.0/n','1.5/n','2.0/n','2.5/n'])
plt.title("Endpoints",fontsize=18)
ax7.set_xlabel("index",fontsize=18)
ax7.set_ylabel("frequency",fontsize=18)
ax7.set_ylim(0,2.5/boxSize)
plt.tight_layout()
fig7.savefig('endpoints.png',dpi=300)

