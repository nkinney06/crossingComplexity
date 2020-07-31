import pandas as pd
import pylab as pl
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

#########################################################################
## distribution of enumerated crossing for snapshots
#########################################################################

directory = sys.argv[1]
datafile = [directory + "/" + x + "/crossingComplexity.txt" for x in os.listdir(directory)]
dataList = list()
for file in datafile:
	dataList.append(pd.read_csv(file, sep='\t', header=0))
dataA = pd.concat([x for x in dataList], ignore_index=True)

data = { sys.argv[1]:dataA["crossings"].tolist() }
df = pd.DataFrame(data)
styles = ['r-']
linewidths = [3]

fig, ax = plt.subplots(figsize=(6,3))
for col, style, lw in zip(df.columns, styles, linewidths):
    df[col].plot.kde(style=style, lw=lw, ax=ax)

pl.legend(fontsize=15)
ax.set_xlabel("Number of Crossings",fontsize=18)
ax.set_ylabel("Density",fontsize=18)
ax.set_xlim(dataA['crossings'].min(),dataA['crossings'].max())
plt.tight_layout()
fig.savefig('density.png',dpi=300)

data = [dataA["crossings"].tolist()]
fig7, ax7 = plt.subplots(figsize=(3,3))
ax7.boxplot(data,labels=[sys.argv[1]],showfliers=False)
ax7.set_ylabel("Number of Crossings",fontsize=14)
ax7.set_xticklabels([sys.argv[1]],fontsize=15,rotation=0)
plt.tight_layout()
fig7.savefig('boxplot.png',dpi=300)
