import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate2d
from scipy.stats import norm
import vtk


def readVTKnoise(fname: str):
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fname)  # declare the vtk filename
    reader.ReadAllVectorsOn()  # read all vector data
    reader.ReadAllScalarsOn()  # read all scalar data
    reader.Update()  # update to new file
    # read data dimensions
    NX, NY, NZ = reader.GetOutput().GetDimensions()
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

    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 24}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)
    plt.rcParams['figure.constrained_layout.use'] = True

    # plot data
    #plt.figure(figsize=(10, 10), dpi=200)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=200)

    mean_xz = np.mean(noiseScalars_xz)
    std_xz = np.std(noiseScalars_xz)
    mean_yz = np.mean(noiseScalars_yz)
    std_yz = np.std(noiseScalars_yz)

    print(f"xz noise mean = {mean_xz}, xz noise standard deviation = {std_xz}")
    print(f"yz noise mean = {mean_yz}, yz noise standard deviation = {std_yz}")

    # Plot histogram with density=False
    #pDensity, bin_edges, _ = ax.hist(noiseScalars_xz, bins=40, density=True, alpha=0.5, label='')
    #bin_widths = np.diff(bin_edges)
    #probability = pDensity*bin_widths
    #print(bin_edges)
    #print(probability)
    #ax.bar(bin_edges[:-1], probability)
    #exit()

# Ca#lculate bin widths
    #bin_widths = np.diff(bin_edges)

# Co#nvert counts to probabilities
    #probabilities = counts / np.sum(counts)

# Pl#ot the probabilities
    #ax.clear()  # Clear the previous histogram
    ##ax.bar(bin_edges[:-1], probabilities, width=bin_widths, alpha=0.5, label='')
    #ax.bar(bin_edges[:-1], probabilities, alpha=0.5, label='')

# Ad#d Gaussian fit (unchanged)
    #xmin, xmax = np.min(noiseScalars_xz), np.max(noiseScalars_xz)
    #gaussRange = np.linspace(xmin, xmax, 100)
    #p = norm.pdf(gaussRange, meanSFE, stdSFE)
    #ax.plot(gaussRange, p, '-k', linewidth=2, alpha=0.4)

# Ad#d annotation (unchanged)
    #ax.annotate(f'$\\sigma$ = {stdSFE:.2f} $Pa^2$', xy=(0, 0.01))

    #plt.show()


    #plt.hist(noiseScalars_xz, bins='auto', density=True, alpha=0.5, label=f'')
    #ax.hist(noiseScalars_xz, bins=20, density=True, alpha=0.5, label=f'')
    ##ax.annotate(f'$\\sigma$ = {stdSFE:.2f} $Pa^2$', xy=(meanSFE, 0), xytext=(meanSFE, 0.1), arrowprops=dict(facecolor='black', shrink=0.05))
    #ax.annotate(f'$\\sigma$ = {stdSFE:.2f} $Pa^2$', xy=(0, 0.01))

    ## add gauss fit
    ##xmin, xmax = plt.xlim()
    #xmin, xmax = np.min(noiseScalars_xz), np.max(noiseScalars_xz)
    #gaussRange = np.linspace(xmin, xmax, 100)
    #p = norm.pdf(gaussRange, meanSFE, stdSFE)
    #ax.plot(gaussRange, p, '-k', linewidth=2, alpha=0.4)
    #plt.show()


if __name__ == "__main__":
    main()
