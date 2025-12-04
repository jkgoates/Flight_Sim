import numpy as np
import math as m
from stl import mesh

PRECISION = 1.0e-5

stlfile = 'Aircraft Carrier simplified.stl'
vtkfile = 'aircraft_carrier.vtk'

myMesh = mesh.Mesh.from_file(stlfile)
meshpoints_orig = myMesh.vectors
nfaces = len(myMesh.vectors)

templines = np.zeros((3*nfaces,2,3))
print("Removing duplicate lines for ", nfaces," faces.")
num_unique = 0
for f in range(nfaces):
    print("CHecking face ", f, "of ", nfaces," faces.")
    for i in range(3):
        if (i+1 < 3):
            ip1 = i+1
        else:
            ip1 = 0
        pa = myMesh.vectors[f,i,:]
        pb = myMesh.vectors[f,ip1,:]
        found = False
        for j in range(num_unique):
            p1 = templines[j,0,:]
            p2 = templines[j,1,:]

            dx1 = pa[0] - p1[0]
            dy1 = pa[1] - p1[1]
            dz1 = pa[2] - p1[2]
            
            dx2 = pb[0] - p2[0]
            dy2 = pb[1] - p2[1]
            dz2 = pb[2] - p2[2]

            if(abs(dx1)+abs(dy1)+abs(dz1)+abs(dx2)+abs(dy2)+abs(dz2) < PRECISION):
                found = True
                break
            else:
                dx1 = pa[0] - p2[0]
                dy1 = pa[1] - p2[1]
                dz1 = pa[2] - p2[2]
            
                dx2 = pb[0] - p1[0]
                dy2 = pb[1] - p1[1]
                dz2 = pb[2] - p1[2]
                if(abs(dx1)+abs(dy1)+abs(dz1)+abs(dx2)+abs(dy2)+abs(dz2) < PRECISION):
                    found = True
                    break
        
        if (not found):
            templines[num_unique,0,:] = pa
            templines[num_unique,1,:] = pb
            num_unique += 1

nlines = num_unique
print("Number of unique lines = ", nlines)
lines = np.zeros((nlines,2))
lines = templines[0:nlines,:,:]

# get unique points
points = np.zeros((1,3))
points[0,:] = lines[0,0,:]

lineid = np.zeros((nlines,2), dtype=int)

print("removing duplicate points for ", nlines," lines.")
num_unique = 1
for f in range(nlines):
    print("Checking line ", f, "of ", nlines," lines")
    for i in range(2):
        pa = lines[f,i,:]
        found = False
        for j in range(num_unique):
            p1 = points[j,:]
            dx = pa[0] - p1[0]
            dy = pa[1] - p1[1]
            dz = pa[2] - p1[2]

            if (abs(dx) + abs(dy) + abs(dz) < PRECISION):
                found = True
                lineid[f,i] = j
                break
        if (not found):
            points = np.append(points, [pa], axis=0)
            lineid[f,i] = num_unique
            num_unique +=1

npoints = len(points)

# optional vtk simplifier

with open(vtkfile, 'w') as writer:
    writer.write('# vtk DataFile Version 3.0\n')
    writer.write(f'{vtkfile}\n')
    writer.write('ASCII\n')
    writer.write('DATASET POLYDATA\n')
    writer.write(f'POINTS {npoints} float\n')
    for i in range(npoints):
        writer.write(f'{5.0*points[i,0]:.6e} {5.0*points[i,1]:.6e} {-5.0*points[i,2]+15:.6e}\n')
    writer.write(f'LINES {nlines} {3*nlines}\n')
    for i in range(nlines):
        writer.write(f'2 {lineid[i,0]} {lineid[i,1]}\n')

print('')
print("Finished writing ", vtkfile)
print('')
print("File Statistics:")
print("--> Original number of points    = ", nfaces*3)
print("--> Original number of lines     = ", nfaces*3)
print("--> VTK number of unique points  = ", npoints)
print("--> VTK number of unique lines   = ", nlines)