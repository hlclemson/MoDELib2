import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def main():
    data = np.loadtxt('./paddedCorrelation.txt')

# draw histogram
#bin = 'auto'
#draw_histogram(data, data2, bin)

# draw heatmap
#data = np.reshape(data, (27, 30))
#data2 = np.reshape(data2, (100, 100))
#data = np.reshape(data, (100, 100))
#data2 = np.reshape(data2, (100, 100))

    data = np.reshape(data, (60,60))
    #data = np.reshape(data, (27,30))
# set plot properties
    font = {'family': 'serif', 'weight': 'normal', 'size': 5}
    mathfont = {'fontset': 'stix'}
    plt.rc('font', **font)
    plt.rc('mathtext', **mathfont)

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
    figName = "paddedCorrelation.png"
    plt.savefig(f'{figName}')
    plt.close()

if __name__ == "__main__":
    main()
