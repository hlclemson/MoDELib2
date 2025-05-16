import matplotlib.pyplot as plt
import sys
import numpy as np
from matplotlib import cm
sys.path.append("../../python")
from modlibUtils import *


def main():
# set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 5}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    #matFile = str(file_templates['material_file']).split('/')[-1]
    matFile = str('AlMg5.txt')
    b_SI=getValueInFile(f'inputFiles/{matFile}', 'b_SI')
    mu0_SI = getValueInFile(f'inputFiles/{matFile}', 'mu0_SI')
    rho_SI = getValueInFile(f'inputFiles/{matFile}', 'rho_SI')
    cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
    convertTimeUnit = b_SI / cs  # [sec]
    convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))
    convertUnitDDDtoSI = (b_SI*mu0_SI)**2

    #data = np.loadtxt('./ensembleCorrelationR100.txt')
    data = np.loadtxt('./test.txt')
    data = np.reshape(data, (60,60))*convertUnitDDDtoSI
    #data = np.reshape(data, (27,30))*convertUnitDDDtoSI
    #data = np.reshape(data, (60,80))*convertUnitDDDtoSI

    #data = np.reshape(data, (27,30))
    fig, axs = plt.subplots(1, 1, figsize=(8,6), dpi=200)  # Adjusted for two subplots
    x = np.arange(data.shape[1])
    y = np.arange(data.shape[0])
    xTickPos = np.arange(len(x))
    yTickPos = np.arange(len(y))
    xTicks = np.linspace(min(x), max(x), num=len(xTickPos))
    yTicks = np.linspace(min(y), max(y), num=len(yTickPos))
    axs.set_xticks(xTickPos)
    axs.set_xticklabels(labels=xTicks, rotation=45, ha="right")
    axs.set_yticks(yTickPos)
    axs.set_yticklabels(labels=yTicks)
    axs.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: '{:.0f}'.format(val)))
    axs.set_title('Original')  # Set the title for each axes
    X = x
    Y = y
    X, Y = np.meshgrid(X, Y)
    Z = data
    c = axs.pcolormesh(X, Y, Z, cmap=cm.coolwarm)
    fig.colorbar(c, ax=axs)
# Save figure 
    figName = "test.png"

    plt.savefig(f'{figName}')
    plt.close()

if __name__ == "__main__":
    main()
