import pyknotid.spacecurves as sp
import pyknotid.make as mk
import numpy as np
import sympy as sym
import os
import matplotlib.pyplot as plt
from operator import itemgetter
import math
import sys

directory = sys.argv[1]

smallAngle = .1
systems = os.listdir(directory)
for system in systems:
    curves = [ directory + "/" + system + "/" + x for x in os.listdir(directory + "/" + system) if "curve" in x ]

    chains = []
    for curve in curves:
        lines = []	
        with open(curve) as file_in:
            for line in file_in:
                lines.append([float(x) for x in line.strip().split()])
        points = []
        for point in lines:
            point = [float(point[0]), float(point[1]), float(point[2])]
            rot1 = [math.cos(smallAngle)*point[0]-math.sin(smallAngle)*point[1],math.cos(smallAngle)*point[1]+math.sin(smallAngle)*point[0],point[2]]
            rot2 = [rot1[0],math.cos(smallAngle)*rot1[1]-math.sin(smallAngle)*rot1[2],math.cos(smallAngle)*rot1[2]+math.sin(smallAngle)*rot1[1]]
            points.append([round(rot2[0],2),round(rot2[1],2),round(rot2[2],2)])
        chains.append(np.array(points))
   
    for chainA in range(len(chains)):
        for chainB in range(chainA+1,len(chains)):
            k = sp.link.Link([chains[chainA],chains[chainB]],verbose=False)
            lnk = k.linking_number()
            print(curves[chainA] + "\t" + curves[chainB] + "\t" + str(lnk))
	
