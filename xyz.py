import sys
import string

coordinates = list()
for line in sys.stdin:
    coordinates.append(line.strip().split())

chains = set(sorted([int(x[3]) for x in coordinates]))

for chain in chains:
    index = 1
    print("CHN " + string.ascii_lowercase[chain])
    for point in coordinates:
        if int(point[3]) == chain:
            print("%d %d.000 %d.000 %d.000" % (index,int(point[0]),int(point[1]),int(point[2])))
            index += 1
print("END\n")
