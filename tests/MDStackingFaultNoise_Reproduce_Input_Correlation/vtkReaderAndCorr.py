import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate2d
import vtk

def readVTKnoise(fname: str):
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fname)  # declare the vtk filename
    reader.ReadAllVectorsOn()  # read all vector data
    reader.ReadAllScalarsOn()  # read all scalar data
    reader.Update()  # update to new file
    # read data dimensions
    dims = [0, 0, 0]
    reader.GetOutput().GetDimensions(dims)
    NX, NY, NZ = dims
    # creates dynamic point data object
    pointData = reader.GetOutput().GetPointData()

    # Initialize a list to store scalar arrays
    scalarArrays = []
    # Iterate over the number of arrays in pointData
    for i in range(pointData.GetNumberOfArrays()):
        array = pointData.GetArray(i)
        # Check if it's a scalar array 
        if array.GetNumberOfComponents() == 1:
            scalarArrays.append(np.array(array))
    noiseScalars = scalarArrays[0]
    noiseScalars = np.reshape(noiseScalars, shape=(NX, NY))
    return noiseScalars


def main():
    dataDir = "../build/OUTPUT/"

    file_xz = f"{dataDir}/ModelNoiseOutput_xz.vtk"
    file_yz = f"{dataDir}/ModelNoiseOutput_yz.vtk"
    noiseScalars_xz = readVTKnoise(file_xz)
    noiseScalars_yz = readVTKnoise(file_yz)

    # Plot the tiled heatmap with coolwarm colormap and colorbar
    plt.figure(figsize=(10, 10), dpi=200)
    im = plt.imshow(noiseScalars_xz, cmap='coolwarm')
    plt.colorbar(im, label='Noise Intensity')  # Add colorbar
    plt.xlabel('X-axis')
    plt.ylabel('Z-axis')
    plt.savefig("sampledNoise_xz.png")
    plt.close()

    # Plot the tiled heatmap with coolwarm colormap and colorbar
    plt.figure(figsize=(10, 10), dpi=200)
    im = plt.imshow(noiseScalars_yz, cmap='coolwarm')
    plt.colorbar(im, label='Noise Intensity')  # Add colorbar
    plt.xlabel('X-axis')
    plt.ylabel('Z-axis')
    plt.savefig("sampledNoise_yz.png")
    plt.close()

    # Plot the tiled heatmap with coolwarm colormap and colorbar
    plt.figure(figsize=(10, 10), dpi=200)
    correlation_xz = correlate2d(noiseScalars_xz, noiseScalars_xz, mode='same')
    im = plt.imshow(correlation_xz, cmap='coolwarm')
    plt.colorbar(im, label='Correlation Intensity')  # Add colorbar
    plt.xlabel('X-axis')
    plt.ylabel('Z-axis')
    plt.savefig("correlationFromNoise_xz.png")
    plt.close()

    # Plot the tiled heatmap with coolwarm colormap and colorbar
    plt.figure(figsize=(10, 10), dpi=200)
    correlation_yz = correlate2d(noiseScalars_yz, noiseScalars_yz, mode='same')
    im = plt.imshow(correlation_yz, cmap='coolwarm')
    plt.colorbar(im, label='Correlation Intensity')  # Add colorbar
    plt.xlabel('X-axis')
    plt.ylabel('Z-axis')
    plt.savefig("correlationFromNoise_yz.png")
    plt.close()


if __name__ == "__main__":
    main()

