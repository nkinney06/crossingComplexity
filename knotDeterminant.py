import pyknotid.spacecurves as sp
import pyknotid.make as mk
import numpy as np
import sympy as sym
import os
import matplotlib.pyplot as plt
from operator import itemgetter
import math
import sys

outfile = open('determinants.txt', 'a+')

directory = sys.argv[1]
datafiles = [directory + "/" + x + "/" + y for x in os.listdir(directory) for y in os.listdir(directory + "/" + x) if "curve" in y]

smallAngle = .1
for curve in datafiles:

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
	points = np.array(points)
	
	try:
		k = sp.OpenKnot(points)
		fractions = k.alexander_fractions(number_of_samples=100)
		print(curve + " " + str(max(fractions,key=itemgetter(1))[0]), file = outfile, flush = True)
	except:
                pass

outfile.close()

	
	
	
