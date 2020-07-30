from numpy import pi, cos, sin, arccos, arange
import sys
if (len(sys.argv) <= 1):
    print("usage: python testDirections.py <int>")
    exit()
num_pts = int(sys.argv[1])
indices = arange(0, num_pts, dtype=float) + 0.5
phi = arccos(1 - 2*indices/num_pts)
theta = pi * (1 + 5**0.5) * indices
x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);
outfile = open("vectors.txt", "w")
for i in range(num_pts):
    outfile.write("{} {} {}\n".format(x[i],y[i],z[i]))
outfile.close()
