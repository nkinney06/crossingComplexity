import sys
import string
import math

coordinates = list()
for line in sys.stdin:
    coordinates.append(line.strip().split())
chains = set(sorted([int(x[3]) for x in coordinates]))
smallAngle = .1
for chain in chains:
    index = 1
    print("CHN " + string.ascii_lowercase[chain])
    for p in coordinates:
        if int(p[3]) == chain:
            point = [float(p[0]), float(p[1]), float(p[2])]
            rot1 = [math.cos(smallAngle)*point[0]-math.sin(smallAngle)*point[1],math.cos(smallAngle)*point[1]+math.sin(smallAngle)*point[0],point[2]]
            rot2 = [rot1[0],math.cos(smallAngle)*rot1[1]-math.sin(smallAngle)*rot1[2],math.cos(smallAngle)*rot1[2]+math.sin(smallAngle)*rot1[1]]
            rot3 = [round(rot2[0],2),round(rot2[1],2),round(rot2[2],2)]
            print("%d %f %f %f" % (index,rot3[0],rot3[1],rot3[2]))
            index += 1
print("END\n")
