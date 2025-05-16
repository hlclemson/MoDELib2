import sys, os, json
import numpy as np
import scipy as scp

sys.path.insert(1,'./src') # add library directory path
import loadRawAtomData as ra
import distributeAtomEnergy as da
import userInputs as ui 
import plot2D as p2
#import plotDist as pd
#import calcCorrelation2D_ti as cc
import countAtom as ca
import calcCorrelation2D as cc
import sample_fr as sf
import vtkDataWrite as vtkdata
import plot1D as p1

import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
from matplotlib import colors
import matplotlib as mpl

def calculateBoxVectors(dumpFile):
    xAxis = dumpFile[0]
    yAxis = dumpFile[1]
    zAxis = dumpFile[2]
    xlowBound, xhighBound, xy = xAxis[0], xAxis[1], xAxis[2]
    ylowBound, yhighBound, xz = yAxis[0], yAxis[1], yAxis[2]
    zlowBound, zhighBound, yz = zAxis[0], zAxis[1], zAxis[2]
    xVectorMagnitude = (xhighBound-xlowBound) - xy
    yVectorMagnitude = np.sqrt((yhighBound-ylowBound)**2 + xy**2)
    return xVectorMagnitude, yVectorMagnitude;

def calcDistFirstNeighbor(atomData: np.ndarray) -> np.float16:
    ## array indexes for atom data
    xAxisIndex = 2
    yAxisIndex = 3
    zAxisIndex = 4
    firstAtomPos = atomData[0][xAxisIndex:zAxisIndex+1]
    numOfAtomsToCheck = 50
    distBetweenAtoms = np.zeros(numOfAtomsToCheck)
    for count, atom in enumerate(atomData):
        atomPos = atom[xAxisIndex:zAxisIndex+1]
        distBetweenAtoms[count] = np.linalg.norm(atomPos-firstAtomPos)
        if count == numOfAtomsToCheck-1:
            break
    # remove zero atom distance in the data
    distBetweenAtoms = distBetweenAtoms[distBetweenAtoms!=0]
    firstNearNeighborDist = np.round(np.min(distBetweenAtoms), 3)
    return firstNearNeighborDist;

def main(verbose, alloy, corCalcMode, startSeed, endSeed):
    # verbose mode function
    vprint = print if verbose else lambda *a, **k: None

    # read SFE MD data directory paths from the config file
    with open('config.json') as f:
        data = json.load(f)
        for key, value in data.items():
            if key=='mainDataDirectory':
                mainDir = value
            elif key=='subDataDirectory':
                dataDir = value
            elif key=='enableZeroPad':
                enableZeroPad = bool(int(value))
            elif key=='gridSizePaddedWithZero':
                gridSizePaddedWithZero = value
            else:
                continue;

    # number of samples
    numOfSamples = endSeed-startSeed+1 if endSeed-startSeed!=0 else 1
    distDataSFE = np.array([])
    seedRange = np.arange(startSeed, endSeed+1);
    for counter, currentSeed in enumerate(seedRange):
        dataPath, nonMinDataPath, outputPath, latParameter = ui.getDataPath(mainDir, dataDir, alloy, currentSeed)
        #vprint(f'dataPath = {dataPath}\n nonMinDataPath = {nonMinDataPath}\n outputPath = {outputPath}\n latParameter = {latParameter}')
        print(f'dataPath = {dataPath}\n nonMinDataPath = {nonMinDataPath}\n outputPath = {outputPath}\n latParameter = {latParameter}')

        # call non-minimized structure for atomID extraction
        nonMinFileList = ra.scanDirs(nonMinDataPath)
        nonMinDumpFileName = ra.loadNonMinStruct(nonMinFileList, alloy)
        atomNonMin = ra.loadAtomDat(nonMinDataPath, nonMinDumpFileName)
        # load box data
        boxNonMin = ra.loadBoxDat(nonMinDataPath, nonMinDumpFileName)
        firstNearNeighborDist = calcDistFirstNeighbor(atomNonMin)
        #vprint(f'boxNonMin = {boxNonMin}')
        #vprint(f'nonMinFileList = {nonMinFileList}')
        #vprint(f'nonMinDumpFileName = {nonMinDumpFileName}')
        #vprint(f'atomNonMin = {atomNonMin}')
        #vprint(f'firstNearNeighborDist = {firstNearNeighborDist}')

        print(f'boxNonMin = {boxNonMin}')
        print(f'nonMinFileList = {nonMinFileList}')
        print(f'nonMinDumpFileName = {nonMinDumpFileName}')
        print(f'atomNonMin = {atomNonMin}')
        print(f'firstNearNeighborDist = {firstNearNeighborDist}')
        # call filelist inside data directory
        fileList = ra.scanDirs(dataPath)
        # call min energy structure fileName(0 : perfect crystal, ISF : intrinsic stacking fault)
        dumpFileName_E0 = ra.findMinimumEStructName(fileList, alloy, Epos='0')
        dumpFileName_EISF = ra.findMinimumEStructName(fileList, alloy, Epos='ISF')
        vprint(f'{dumpFileName_E0} {dumpFileName_EISF}')

        # load raw atom data
        atomE0 = ra.loadAtomDat(dataPath, dumpFileName_E0)
        atomEISF = ra.loadAtomDat(dataPath, dumpFileName_EISF)

        # degree from 110 to 011
        theta = 60*(np.pi/180) # from 110 x axis
        # Grid Parameters
        xGridWidth = firstNearNeighborDist-.001
        yGridWidth = firstNearNeighborDist-.001
        # basis 
        basis1 = np.array([1, 0]) * xGridWidth
        basis2 = np.array([np.cos(theta), np.sin(theta)]) * yGridWidth
        vprint(f'basis1 = {basis1}, basis2 = {basis2}')

        # calculate the num of grids
        xVectorMagnitude, yVectorMagnitude = calculateBoxVectors(boxNonMin)
        vprint(f'xVector = {xVectorMagnitude} yVector = {yVectorMagnitude}')
        # one grid number is subtracted since A stacking of atoms do not exist 
        # at the end of x and y boundary (handled by periodic BC)
        xGridNum = int(np.floor(xVectorMagnitude/xGridWidth))
        yGridNum = int(np.floor(yVectorMagnitude/yGridWidth))
        #xGridNum = int(np.ceil(xVectorMagnitude/xGridWidth))
        #yGridNum = int(np.ceil(yVectorMagnitude/yGridWidth))
        #print(f'xGridNum = {xGridNum}, yGridNum = {yGridNum}')

        atomPerGrid = da.countAtomsPerGrid(atomNonMin, atomE0, atomEISF, xGridNum, yGridNum, xGridWidth, yGridWidth, basis1, basis2, theta, debug=0)
        isNumOfAtomsPerGridEqual = np.all(atomPerGrid==atomPerGrid[0])
        #print(f'atomPerGrid = {atomPerGrid}')
        if not isNumOfAtomsPerGridEqual:
            exit("Number of atoms per grid is not equal")

        # distribute and calculate GSFE
        Egrid = da.calcGSFE(atomNonMin, atomE0, atomEISF, xGridNum, yGridNum, xGridWidth, yGridWidth, basis1, basis2, theta, debug=0)

        # calculate input function \gamma - \gamma_mean
        F = Egrid - np.mean(Egrid)

        # print standard deviation of F

        # pad the input data with zeros (if enabled)
        if enableZeroPad:
            xNewGridSize, yNewGridSize = gridSizePaddedWithZero[0], gridSizePaddedWithZero[1]
            xGridSize, yGridSize = F.shape[1], F.shape[0]
            numOfZerosToPad_x = int((xNewGridSize-xGridSize)/2)
            numOfZerosToPad_y = int((yNewGridSize-yGridSize)/2)
            F = np.pad(F,((numOfZerosToPad_y, numOfZerosToPad_x)), constant_values=0)

        # accumulate all the SFE data in a single table for distribution study
        if counter==0:
            distDataSFE = F
        else:
            distDataSFE = np.vstack((distDataSFE, F))

        if corCalcMode=='FFT':
            Fhat = np.fft.fft2(F)
            Chat_i = Fhat*np.conjugate(Fhat)
            # for ensemble average
            if counter==0:
                # record the grid size of the first seed
                firstSeedXgridSize = Chat_i.shape[0]
                firstSeedYgridSize = Chat_i.shape[1]
                # initialize the ensemble average grid
                ChatAVG = np.zeros((Chat_i.shape[0], Chat_i.shape[1]), dtype=np.complex_)
                # add calculated correlation
                ChatAVG += Chat_i
            else:
                # skip the seed if the grid size does not match
                if Chat_i.shape[0] != firstSeedXgridSize or Chat_i.shape[1] != firstSeedYgridSize:
                    numOfSamples -= 1
                    vprint(f'seed number {currentSeed} is skipped')
                    continue
                ChatAVG += Chat_i
        elif corCalcMode=='Real':
            # calculate spatial correlation in real space
            Cx_i = cc.calAutoCor2D(F)
            vprint(f'{Cx_i}')
            # for ensemble average
            if counter==0:
                # record the grid size of the first seed
                firstSeedXgridSize = Cx_i.shape[0]
                firstSeedYgridSize = Cx_i.shape[1]
                # initialize the ensemble average grid
                Cx_AVG = np.zeros((Cx_i.shape[0], Cx_i.shape[1]), dtype=np.float128)
                # add calculated correlation
                Cx_AVG += Cx_i
                # for errorBar 1D plot (standard error), extract 1 row/column
                dxCorr = Cx_i[0,:]
                dyCorr = Cx_i[:,0]
            else:
                # skip the seed if the grid size does not match
                if Cx_i.shape[0] != firstSeedXgridSize or Cx_i.shape[1] != firstSeedYgridSize:
                    numOfSamples -= 1
                    print(f'seed number {currentSeed} is skipped')
                    continue
                Cx_AVG += Cx_i
                # for errorBar 1D plot (standard error), extract 1 row/column, stack them vertically
                dxCorr = np.vstack((dxCorr, Cx_i[0,:]))
                dyCorr = np.vstack((dyCorr, Cx_i[:,0]))
        #elif corCalcMode=='Scipy':
        else:
            Cx_i = scp.signal.correlate2d(F, F, mode='same')
            #Cx_i = scp.signal.correlate2d(F, F, mode='same')/np.sqrt(np.sum(F**2)*np.sum(F**2)) # normalized correlation
            vprint(f'{Cx_i}')
            # for ensemble average
            if counter==0:
                # record the grid size of the first seed
                firstSeedXgridSize = Cx_i.shape[0]
                firstSeedYgridSize = Cx_i.shape[1]
                # initialize the ensemble average grid
                Cx_AVG = np.zeros((Cx_i.shape[0], Cx_i.shape[1]), dtype=np.float128)
                # add calculated correlation
                Cx_AVG += Cx_i
                # for errorBar 1D plot (standard error), extract 1 row/column
                dxCorr = Cx_i[0,:]
                dyCorr = Cx_i[:,0]
            else:
                # skip the seed if the grid size does not match
                if Cx_i.shape[0] != firstSeedXgridSize or Cx_i.shape[1] != firstSeedYgridSize:
                    numOfSamples -= 1
                    print(f'seed number {currentSeed} is skipped')
                    continue
                Cx_AVG += Cx_i
                # for errorBar 1D plot (standard error), extract 1 row/column, stack them vertically
                dxCorr = np.vstack((dxCorr, Cx_i[0,:]))
                dyCorr = np.vstack((dyCorr, Cx_i[:,0]))

    # create data folder
    dataFolder = './data'
    if not os.path.exists(dataFolder):
        os.system(f'mkdir {dataFolder}')

    # Save stacking fault energy data for distribution study
    np.savetxt(f'./data/{alloy}_F_R{numOfSamples}.txt', distDataSFE)

    # Construct non-orthogonal dot grid based on the non-minimized atom structure
    xGridNumWithPad = F.shape[1]
    yGridNumWithPad = F.shape[0]
    dotGrid = np.zeros((xGridNumWithPad*yGridNumWithPad, 2)) # each row contains x and y position of the dot
    unitBasis1Vector = np.array([1, 0])*firstNearNeighborDist;
    unitBasis2Vector = np.array([np.cos(np.deg2rad(60)), np.sin(np.deg2rad(60))])*firstNearNeighborDist;
    for yIdx in range(yGridNumWithPad):
        for xIdx in range(xGridNumWithPad):
            dotGrid[yIdx*xGridNumWithPad + xIdx] = unitBasis1Vector*xIdx + unitBasis2Vector*yIdx

    # average the correlation data
    if corCalcMode=='FFT':
        ChatAVG /= numOfSamples
        # inverse fourier transform of averaged data in fourier space
        Fsq = np.fft.ifft2(ChatAVG)
        #Fsq = Fsq.real
        Fsq = Fsq.real /np.sqrt(np.sum(F**2)*np.sum(F**2)) # normalized correlation
    #elif corCalcMode=='Real' or corCalcMode=='Scipy':
    else:
        Cx_AVG /= numOfSamples

    if corCalcMode=='FFT':
        # save data in VTK format
        #vtkdata.saveInVTKstructGrid(atomNonMin, xGridNum, yGridNum, latParameter, xVectorMagnitude, yVectorMagnitude, Fsq, f'{alloy}_CxFFT_R{numOfSamples}.vtk')
        vtkdata.saveInVTKstructGrid2(dotGrid, xGridNumWithPad, yGridNumWithPad, latParameter, xVectorMagnitude, yVectorMagnitude, Fsq, f'{alloy}_CxFFT_R{numOfSamples}.vtk')
        # save spatial correlation data 
        np.savetxt(f'./data/{alloy}_CxFFT_{numOfSamples}.txt', Fsq)
        # plot spatial correlation
        dx = np.arange(Fsq.shape[1])
        dy = np.arange(Fsq.shape[0])
        title = f'$C(x)$ from FFT, {numOfSamples} realizations'
        plotlabels = {'title': title, 'xlabel': r'$dx$', 'ylabel': r'$dy$'}
        p2.mapPlot(dx, dy, Fsq, plotlabels, f'{alloy}_correlationFFT_{numOfSamples}R.png', dataDir)
        plotlabels = {'title': f'{alloy}', 'xlabel': r'$dx$', 'ylabel': r'$dy$'}
        #p2.drawDotCorrelationPlot(atomNonMin, basis1, basis2, Fsq, plotlabels, f'{alloy}_dot_correlationFFT_{numOfSamples}R.png', dataDir)
        p2.drawDotCorrelationPlot(dotGrid, basis1, basis2, Fsq, plotlabels, f'{alloy}_dot_correlationFFT_{numOfSamples}R.png', dataDir)

    elif corCalcMode=='Real' or corCalcMode=='Scipy':
        # save data in VTK format
        #vtkdata.saveInVTKstructGrid(atomNonMin, xGridNum, yGridNum, latParameter, xVectorMagnitude, yVectorMagnitude, Cx_AVG, f'{alloy}_Cx_R{numOfSamples}.vtk')
        if corCalcMode=='Real':
            vtkdata.saveInVTKstructGrid2(dotGrid, xGridNumWithPad, yGridNumWithPad, latParameter, xVectorMagnitude, yVectorMagnitude, Cx_AVG, f'{alloy}_Cx_R{numOfSamples}.vtk')
        else:
            vtkdata.saveInVTKstructGrid2(dotGrid, xGridNumWithPad, yGridNumWithPad, latParameter, xVectorMagnitude, yVectorMagnitude, Cx_AVG, f'{alloy}_Cx_R{numOfSamples}_S.vtk')
        # save spatial correlation data 
        if corCalcMode=='Real':
            np.savetxt(f'./data/{alloy}_Cx_R{numOfSamples}.txt', Cx_AVG)
        else:
            np.savetxt(f'./data/{alloy}_Cx_R{numOfSamples}_S.txt', Cx_AVG)
        if endSeed != 1:
            # calculate standard error
            dxCorrSTE = np.std(dxCorr, axis=0)/np.sqrt(numOfSamples)
            dyCorrSTE = np.std(dyCorr, axis=0)/np.sqrt(numOfSamples)
            # save standard error data
            if corCalcMode=='Real':
                np.savetxt(f'./data/{alloy}_STE_x_{numOfSamples}.txt', dxCorrSTE)
                np.savetxt(f'./data/{alloy}_STE_y_{numOfSamples}.txt', dyCorrSTE)
            else:
                np.savetxt(f'./data/{alloy}_STE_x_{numOfSamples}_S.txt', dxCorrSTE)
                np.savetxt(f'./data/{alloy}_STE_y_{numOfSamples}_S.txt', dyCorrSTE)
        # plot spatial correlation
        dx = np.arange(Cx_AVG.shape[1])
        dy = np.arange(Cx_AVG.shape[0])
        title = f'$ C(x) $ from Real Space, {numOfSamples} realizations'
        plotlabels = {'title': title, 'xlabel': r'$dx$', 'ylabel': r'$dy$'}
        if corCalcMode=='Real':
            p2.mapPlot(dx, dy, Cx_AVG, plotlabels, f'{alloy}_correlation_{numOfSamples}R.png', dataDir)
        else:
            p2.mapPlot(dx, dy, Cx_AVG, plotlabels, f'{alloy}_correlation_{numOfSamples}R_S.png', dataDir)
        plotlabels = {'title': f'{alloy}', 'xlabel': r'$dy$', 'ylabel': r'$dx$'}
        #p2.drawDotCorrelationPlot(atomNonMin, basis1, basis2, Cx_AVG, plotlabels, f'{alloy}_dot_correlation_{numOfSamples}R.png', dataDir)
        if corCalcMode=='Real':
            p2.drawDotCorrelationPlot(dotGrid, basis1, basis2, Cx_AVG, plotlabels, f'{alloy}_dot_correlation_{numOfSamples}R.png', dataDir)
        else:
            p2.drawDotCorrelationPlot(dotGrid, basis1, basis2, Cx_AVG, plotlabels, f'{alloy}_dot_correlation_{numOfSamples}R_S.png', dataDir)

def printHelpMessage(availableAlloys: list, availCorrelationCalcModes: list):
    print(f'code.py [options] [alloyName] [correlation calc mode] [startSeed] [endSeed]')
    print(f'Options: ')
    print(f'    -v: verbose mode, print out debug information')
    print(f'    -h: print help message')
    print(f'available alloys : {str(availableAlloys).strip("[]")}')
    print(f'available correlation calculation mode: {str(availCorrelationCalcModes).strip("[]")}')
    print(f'       - Real  : Calculate spatial correlation in real space')
    print(f'       - FFT   : Calculate spatial correlation with FFT (Wiener-Khinchin Theorem)')
    print(f'       - Scipy : Calculate spatial correlation in real space and shift the zero frequency element to the center')
    print(f'Example: ')
    print(f'       - code.py AlMg5 FFT 1 10 ')
    print(f'       - code.py AlMg10 Real 1 10 ')
    print(f'       - code.py AlMg10 Scipy 1 10 ')

if __name__=='__main__':
    availAlloys = ['AlMg5', 'AlMg10', 'AlMg15', 'AlMg25', 'AlMg50', 'NiCoCr0', 'NiCoCr350', 'NiCoCr650', 'NiCoCr950','NiCoCr1250']
    availCorrelationCalcModes = ['Real', 'Scipy', 'FFT']
    # default verbose value
    verbose = False
    alloy = None
    corCalcMode = None
    startSeed = None
    endSeed = None
    for arg in sys.argv[1:]:
        if arg=='-v':
            verbose = True
        elif arg=='-h':
            printHelpMessage(availAlloys, availCorrelationCalcModes)
            sys.exit(1)
        elif arg in availAlloys:
            alloy = arg
        elif arg in availCorrelationCalcModes:
            corCalcMode = arg
        else:
            if startSeed==None:
                startSeed = int(arg)
            else:
                endSeed = int(arg)

    if alloy==None or corCalcMode==None or startSeed==None or endSeed==None:
        printHelpMessage(availAlloys, availCorrelationCalcModes)
        sys.exit(1)

    # run the main program
    main(verbose, alloy, corCalcMode, startSeed, endSeed)


