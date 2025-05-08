import vtk
import numpy as np

#print(dir(vtk))


#vtk.vtkArrayDataReader('AlMg5_ISF.vtk')
#vtk.vtkArrayReader('AlMg5_ISF.vtk')

#reader = vtk.vtkUnstructuredGridReader()
#reader = vtk.vtkStructuredPointsReader()
reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName('AlMg5_ISF.vtk')
reader.ReadAllVectorsOn()  # read all vector data
reader.ReadAllScalarsOn()  # read all scalar data
reader.Update()
NX,NY,NZ = reader.GetOutput().GetDimensions()
pointData = reader.GetOutput().GetPointData()
scalarArrays = []
# Iterate over the number of arrays in pointData
for i in range(pointData.GetNumberOfArrays()):
    array = pointData.GetArray(i)
    # Check if it's a scalar array 
    if array.GetNumberOfComponents() == 1:
        scalarArrays.append(np.array(array))
noiseScalars = scalarArrays[0]
noiseScalars = np.reshape(noiseScalars, shape=(NX, NY))
print(noiseScalars)


#reader = vtk.vtkGenericDataObjectReader()
#reader.SetFileName(fname)  # declare the vtk filename
#reader.Update()  # update to new file

# read data dimensions
#NX, NY, NZ = reader.GetOutput().GetDimensions()
## creates dynamic point data object
#pointData = reader.GetOutput().GetPointData()
#
## Initialize a list to store scalar arrays
#scalarArrays = []
## Iterate over the number of arrays in pointData
#for i in range(pointData.GetNumberOfArrays()):
#    array = pointData.GetArray(i)
#    # Check if it's a scalar array 
#    if array.GetNumberOfComponents() == 1:
#        scalarArrays.append(np.array(array))
#noiseScalars = scalarArrays[0]
#noiseScalars = np.reshape(noiseScalars, shape=(NX, NY))
#return noiseScalars

