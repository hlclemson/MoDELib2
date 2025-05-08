#!/bin/python3

import matplotlib.pyplot as plt
from matplotlib import cm
#import seaborn as sns
import numpy as np
from scipy.signal import correlate2d
import os
import sys
sys.path.append("../../python")
from modlibUtils import *

#def draw_histogram(data, data2, bin):
#    # set plot properties
#    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
#    mathfont = {'fontset': 'stix'}
#    plt.rc('font', **font)
#    plt.rc('mathtext', **mathfont)
#
#    # Create a figure with two subplots (1 row, 2 columns)
#    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots
#
#    # Set common properties for both plots
#    plotLabels = {'title': '', 'xlabel': 'SFE $(mJ/m^2)$', 'ylabel': 'Probability Density'}
#    #colors = sns.color_palette("husl")
#
#    # Plotting on the first subplot
#    axs[0].grid(True)
#    axs[0].set_title(f"{plotLabels['title']}")
#    axs[0].set_xlabel(f"{plotLabels['xlabel']}")
#    axs[0].set_ylabel(f"{plotLabels['ylabel']}")
#    #axs[0].hist(data, facecolor=colors[0], edgecolor='k', density=True, bins=bin)
#    axs[0].hist(data, density=True, bins=bin)
#    axs[0].set_title('Original')
#
#    # Plotting on the second subplot
#    axs[1].grid(True)
#    axs[1].set_title(f"{plotLabels['title']}")
#    axs[1].set_xlabel(f"{plotLabels['xlabel']}")
#    axs[1].set_ylabel(f"{plotLabels['ylabel']}")
#    #axs[1].hist(data2, facecolor=colors[1], edgecolor='k', density=True, bins=bin)  # Changed color for distinction
#    axs[1].hist(data2, density=True, bins=bin)  # Changed color for distinction
#    axs[1].set_title('Sampled (n=1000)')
#
#    # Save figure 
#    figName = "histogram.png"
#    figFolder = './figures'
#    if not os.path.exists(figFolder):
#        os.makedirs(figFolder)  # Changed to os.makedirs for better practice
#    plt.savefig(f'{figFolder}/{figName}')
#    plt.close()
#
#def drawHeatMap(data, data2):
#    # set plot properties
#    font = {'family': 'serif', 'weight': 'normal', 'size': 5}
#    mathfont = {'fontset': 'stix'}
#    plt.rc('font', **font)
#    plt.rc('mathtext', **mathfont)
#
#    # Create a figure with two subplots (1 row, 2 columns)
#    fig, axs = plt.subplots(1, 2, figsize=(16,6), dpi=200)  # Adjusted for two subplots
#    x = np.arange(data.shape[1])
#    y = np.arange(data.shape[0])
#    xTickPos = np.arange(len(x))
#    yTickPos = np.arange(len(y))
#    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
#    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
#    axs[0].set_xticks(xTickPos)
#    axs[0].set_xticklabels(labels=xTicks, rotation=45, ha="right")
#    axs[0].set_yticks(yTickPos)
#    axs[0].set_yticklabels(labels=yTicks)
#    axs[0].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
#    axs[0].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
#    axs[0].set_title('Original')  # Set the title for each axes
#    X = x
#    Y = y
#    X, Y = np.meshgrid(X, Y)
#    Z = data
#    c = axs[0].pcolormesh(X, Y, Z, cmap=cm.viridis)
#    fig.colorbar(c, ax=axs[0])
#
#    x = np.arange(data2.shape[1])
#    y = np.arange(data2.shape[0])
#    xTickPos = np.arange(len(x))
#    yTickPos = np.arange(len(y))
#    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
#    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
#    axs[1].set_xticks(xTickPos)
#    axs[1].set_xticklabels(labels=xTicks, rotation=45, ha="right")
#    axs[1].set_yticks(yTickPos)
#    axs[1].set_yticklabels(labels=yTicks)
#    axs[1].xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
#    axs[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
#    axs[1].set_title('Sampled (n=1000)')  # Set the title for each axes
#    X = x
#    Y = y
#    X, Y = np.meshgrid(X, Y)
#    Z = data2
#    c = axs[1].pcolormesh(X, Y, Z, cmap=cm.viridis)
#    fig.colorbar(c, ax=axs[1])
#
#    # Save figure 
#    figName = "heatmap.pdf"
#    figFolder = './figures'
#    if not os.path.exists(figFolder):
#        os.makedirs(figFolder)  # Changed to os.makedirs for better practice
#    plt.savefig(f'{figFolder}/{figName}')
#    plt.close()
#
#def shiftCenter(data: np.ndarray) -> np.ndarray:
#    # swap the first half of the rows with the second half
#    dummyArray = data.copy()
#    dummyArray[:int(data.shape[0]/2), :], dummyArray[int(data.shape[0]/2):, :] = data[int(data.shape[0]/2):, :], data[:int(data.shape[0]/2), :]
#    # swap the first half of the columns with the second half
#    dummyArray2 = dummyArray.copy()
#    dummyArray2[:, int(data.shape[1]/2):], dummyArray2[:, :int(data.shape[1]/2)] = dummyArray[:, :int(data.shape[1]/2)], dummyArray[:, int(data.shape[1]/2):]
#    return dummyArray2;

def main():
    matFile = str('AlMg5.txt')
    b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
    rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')
    cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
    convertTimeUnit = b_SI / cs  # [sec]
    convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))
    convertUnitDDDtoSI = b_SI*mu0_SI

    # Load data
    #data = np.loadtxt('./noiseDistributionR10.txt')*convertUnitDDDtoSI
    data = np.loadtxt('./noiseDistribution.txt')*convertUnitDDDtoSI

    mean = np.mean(data)
    std = np.std(data)
    print(f"mean = {mean} J/m^2")
    print(f"std = {std} J/m^2")

    # draw histogram
    #bin = 'auto'
    #draw_histogram(data, data2, bin)

    # draw heatmap
    #data = np.reshape(data, (27, 30))
    #data2 = np.reshape(data2, (100, 100))
    #data = np.reshape(data, (100, 100))
    #data2 = np.reshape(data2, (100, 100))

    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 5}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    # Plotting on the first subplot
    fig, axs = plt.subplots(1, 1, figsize=(8,6), dpi=200)  # Adjusted for two subplots
    axs.grid(True)
    axs.set_title(f"")
    axs.set_xlabel(f"")
    axs.set_ylabel(f"")
    #axs[0].hist(data, facecolor=colors[0], edgecolor='k', density=True, bins=bin)
    #axs.hist(data, density=True, bins=bin)
    axs.hist(data, density=True, bins='auto')
    fig.savefig(f'testNoiseDistribution.png')
    plt.close()


if __name__ == "__main__":
    main()
