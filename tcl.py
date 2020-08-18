import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull
import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import re
import subprocess

zoom = .08
color = "Fragment"
size = .1

outfile = open("vmd.tcl", "w")
outfile.write("mol new ./system.vtf\n")
outfile.write("mol representation Licorice " + str(size) + " 42\n")
outfile.write("mol color " + color + "\n")
outfile.write("mol material AOChalky\n")
outfile.write("mol addrep 0\n")
outfile.write("color Display Background 8\n")
outfile.write("rotate x by -90\n")
outfile.write("rotate y by -90\n")
outfile.write("rotate y by -35\n")
outfile.write("rotate x by 35\n")
outfile.write("axes location off\n")
outfile.write("scale to " + str(zoom) + "\n")
outfile.write("render Tachyon structure.dat\n")
outfile.write("quit\n")
outfile.close()

some_command = "vmd -dispdev text -e vmd.tcl"
p = subprocess.Popen(some_command, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()  
p_status = p.wait()

some_command = "tachyon structure.dat -o curve.bmp -format BMP -res 1050 1050"
p = subprocess.Popen(some_command, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()  
p_status = p.wait()

############################################################################################################################
## main programming section
############################################################################################################################

# read the structure.dat file
lines = []
with open("structure.dat") as file_in:
    for line in file_in:
        lines.append(line)

# read the points of one of the space filling curves
colors = [ lines[x+5] for x in range(len(lines)) if "Sphere" in lines[x] ]
colors = list(set([ re.search(r'Color (.*) Tex',x).group(1) for x in colors ]))
sphere = [ lines[x+1] for x in range(len(lines)) if "Sphere" in lines[x] and colors[0] in lines[x+5] ]
sphere = [ re.search(r'Center ([+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)? [+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)? [+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?)',x).group(1) for x in sphere ]
sphere = [ [float(x) for x in y.split()] for y in sphere ]
pts = np.array(sphere)

# compute the convex hull
hull = ConvexHull(pts)

outfile = open("territory.dat", "w")
for line in lines:
        if "End_Scene" not in line:
                outfile.write(line.rstrip() + "\n")

for s in hull.simplices:
        outfile.write('\nTRI\n')
        outfile.write('  V0 {} {} {}\n'.format(pts[s[0],0],pts[s[0],1],pts[s[0],2]))
        outfile.write('  V1 {} {} {}\n'.format(pts[s[1],0],pts[s[1],1],pts[s[1],2]))
        outfile.write('  V2 {} {} {}\n'.format(pts[s[2],0],pts[s[2],1],pts[s[2],2]))
        outfile.write('  TEXTURE\n')
        outfile.write('    AMBIENT 0.1 DIFFUSE 0.85 SPECULAR 0 OPACITY .7\n')
        outfile.write('    COLOR 0.0 0.0 0.5\n')
        outfile.write('    TEXFUNC 0\n')
outfile.write("\nEnd_Scene\n")
outfile.close()

some_command = "tachyon territory.dat -o chull.bmp -format BMP -res 1050 1050"
p = subprocess.Popen(some_command, stdout=subprocess.PIPE, shell=True)
(output, err) = p.communicate()  
p_status = p.wait()

