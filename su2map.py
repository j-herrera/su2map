# Copyright 2017 Javier Herrera Montojo
# j.herrera.montojo@gmail.com

# Imports
import sys
import numpy as np
import warnings
import vtk
from vtk.util import numpy_support

np.seterr(all='warn')

# Test inputs (python su2map.py flow.vtk new_mesh.su2 solution_flow.dat)
if len(sys.argv) != 4:
	print("Program takes three arguments. Input vtk flowfield, new mesh, and output restart file")
	sys.exit()

# Read flow
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(sys.argv[1])
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

data = reader.GetOutput()

p_data = data.GetPointData()

n_arr = p_data.GetNumberOfArrays()

arrays = []
for i in range(n_arr):
    arrays.append(p_data.GetArrayName(i))

# Read mesh
p_su2 = open(sys.argv[2], "r")

ln = p_su2.readline()
ndim =  int(ln.split('=')[1])
ln = p_su2.readline()
nelem = int(ln.split('=')[1])

cells = {}
for i in range(nelem):
	ln = p_su2.readline()
	cells[i] = np.array(np.fromstring(ln, dtype=np.int, sep=' '))
	
ln = p_su2.readline()
npoin = int(ln.split('=')[1])

points = np.zeros((npoin, ndim+1))
for i in range(npoin):
	ln = p_su2.readline()
	points[i,:] = np.array(np.fromstring(ln, dtype=np.float, sep=' '))

p_su2.close()

new_mesh = vtk.vtkUnstructuredGrid()
mpoints = vtk.vtkPoints()
mcells = vtk.vtkCellArray()

mpoints.SetNumberOfPoints(npoin)
for i in range(npoin):
    tmp_pt = points[i,:].copy()
    if (ndim==2):
        tmp_pt = np.array([points[i,0], points[i,1], 0.0])
    else:
        tmp_pt = np.array([points[i,0], points[i,1], points[i,2]])
    mpoints.SetPoint(np.int(points[i,-1]), tmp_pt)

for i in range(nelem):
    idlist = vtk.vtkIdList()
    for j in range(cells[i].size-2):
        idlist.InsertNextId(cells[i][j+1])
    mcells.InsertNextCell(idlist)
    

new_mesh.SetPoints(mpoints)
new_mesh.SetCells(cells[i][0], mcells)

# Probe
probe = vtk.vtkProbeFilter()
probe.SetInput(new_mesh)
probe.SetSource(data)
probe.Update()

probe_data = probe.GetOutputDataObject(0).GetPointData()
valid = np.where(numpy_support.vtk_to_numpy(probe_data.GetAbstractArray(probe.GetValidPointMaskArrayName())) == 0)[0]

# Filter out invalid points (by mean of connected cells)
abstract = {}
cellid = vtk.vtkIdList()
for name in arrays:
    abstract[name] = numpy_support.vtk_to_numpy(probe_data.GetAbstractArray(name))
    for i in valid:
        cellid.Reset()
        new_mesh.GetPointCells(i, cellid)
        nc =  cellid.GetNumberOfIds()
        neighbours = []
        for j in range(nc):
            cell = new_mesh.GetCell(cellid.GetId(j))
            for k in range(cell.GetNumberOfPoints()):
                if not cell.GetPointId(k) in valid:
                    neighbours.append(np.int(cell.GetPointId(k)))
        neighbours = np.unique(neighbours)
        abstract[name][i] = np.sum(abstract[name][neighbours])/len(neighbours)
            
# Export
pos = open(sys.argv[3], "w")
pos.write('"PointID" ')
pos.write('"x" ')
pos.write('"y" ')
if (ndim==3):
    pos.write('"z" ')

for name in arrays:
	pos.write('"'+name+'" ')
pos.write('\n')

for i in range(npoin):
	pos.write(str(i))
	for j in range(ndim):
		pos.write(' ' + str(points[i,j]))
	for name in arrays:
		pos.write('{:e}'.format(abstract[name][i]))
		pos.write(' ')
	pos.write('\n')

pos.write('\n')

pos.close()

	



