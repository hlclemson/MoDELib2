import vtk
import numpy as np
from scipy.signal import correlate2d

#print(dir(vtk))


#vtk.vtkArrayDataReader('AlMg5_ISF.vtk')
#vtk.vtkArrayReader('AlMg5_ISF.vtk')

#reader = vtk.vtkUnstructuredGridReader()
#reader = vtk.vtkStructuredPointsReader()
reader = vtk.vtkGenericDataObjectReader()
reader.SetFileName('sf_noise_AlMg5.vtk')
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
print(f"np.max(noiseScalars) = {np.max(noiseScalars)}")
print(f"np.min(noiseScalars) = {np.min(noiseScalars)}")
print(noiseScalars)
print(noiseScalars.shape)

correlation = correlate2d(noiseScalars, noiseScalars, mode='same')
# Plot the tiled heatmap with coolwarm colormap and colorbar
plt.figure(figsize=(10, 10), dpi=200)
im = plt.imshow(correlation, cmap='coolwarm')
plt.colorbar(im, label='Noise Intensity')  # Add colorbar
plt.xlabel('X-axis')
plt.ylabel('Z-axis')
#plt.show()
plt.savefig("noise.png")


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

