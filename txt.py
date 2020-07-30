import sys
import string

coordinates = list()
for line in sys.stdin:
    coordinates.append(line.strip().split())

for curve in set([int(coordinates[x][3]) for x in range(len(coordinates))]):
    outfile = open("curve_"+str(curve+1)+".txt", "w")
    for point in range(len(coordinates)):
        if ( int(coordinates[point][3])==curve ):
            outfile.write("%d.000 %d.000 %d.000\n" % (int(coordinates[point][0]),int(coordinates[point][1]),int(coordinates[point][2])))
    outfile.close()

