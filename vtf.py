import sys
import string

coordinates = list()
for line in sys.stdin:
    coordinates.append(line.strip().split())

print("\nunitcell 100.0 100.0 100.0")
for atom in range(len(coordinates)):
        print("atom %d radius 1 name %d type %d" % (atom,int(coordinates[atom][3]),int(coordinates[atom][3])))
for atom in range(1,len(coordinates)):
    if (coordinates[atom-1][3]!=coordinates[atom][3]):
        continue
    print("bond %d:%d" % (atom-1,atom))

print("\ntimestep indexed")
for atom in range(len(coordinates)):
    print("%d %d %d %d" % (atom,int(coordinates[atom][0]),int(coordinates[atom][1]),int(coordinates[atom][2])))


