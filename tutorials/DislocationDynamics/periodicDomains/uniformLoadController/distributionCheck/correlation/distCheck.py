import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    # set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)
    fig, ax = plt.figure(figsize=(8,6), dpi=200), plt.subplot(1,1,1)
    ax.grid(True)
    plotLabels = {'title': f'', 'xlabel': 'SFE $(mJ/m^2)$', 'ylabel': 'Probability Density'}
    title = f"{plotLabels['title']}"
    xlabel = f"{plotLabels['xlabel']}"
    ylabel = f"{plotLabels['ylabel']}"
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # set plot properties
    colors = sns.color_palette("husl")
    font = {'family': 'serif', 'weight': 'normal', 'size': 12}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

    data = np.loadtxt('correlation_recovered.txt')
    data2 = np.loadtxt('original_Correlation.txt')
    bin = 'auto'
    #bin = 10
    ax.hist(data, facecolor=colors[0], edgecolor='k', density=True, bins=bin)
    ax.hist(data2, facecolor=colors[0], edgecolor='k', density=True, bins=bin)

    # save figure 
    figName = f"hist_discretized.png"
    figFolder = './figures'
    if not os.path.exists(figFolder):
        os.system(f'mkdir {figFolder}')
    plt.savefig(f'{figFolder}/{figName}')
    plt.close()

main()
