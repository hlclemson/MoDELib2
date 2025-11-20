import re
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams


# Configure global plot settings (applies to all figures)
rcParams.update({
    'figure.dpi': 200,
    'figure.autolayout': True,  # Prevent label clipping
    'axes.grid': True,
    'grid.alpha': 0.5,
    'text.usetex': True,   # Enable LaTeX
    'font.size': 30,     # Default font size for text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    'font.weight': 'bold',          # Base font weight
    'text.latex.preamble': r'\usepackage{amsmath} \boldmath'  # Additional packages
})

#def plot_k_data_avg(fname: str) -> None:
#    df = pd.read_csv(fname)
#
#    dchars = df['dislocationCharacter'].unique()
#    noise_types = df['noiseType'].unique()
#
#    for dchar in dchars:
#        fig, ax = plt.subplots(figsize=(11, 8))
#        for noise in noise_types:
#            data = df[(df['dislocationCharacter'] == dchar) & (df['noiseType'] == noise)]
#
#            # Extract and convert alloy values to integers
#            x_values = data['alloy'].str.extract(r'AlMg(\d+)').astype(int).squeeze()
#            zeta = -(data['kExponent'] + 1) / 2
#
#            # Group by unique x_values and average y_values
#            grouped = data.groupby(x_values)
#            avg_zeta = grouped.apply(lambda x: -(x['kExponent'] + 1).mean() / 2)
#
#            # Plot the averaged values
#            ax.plot(avg_zeta.index, avg_zeta.values, 'o-', label=noise)
#
#        ax.set_xlabel('Alloy Composition (Mg%)')
#        ax.set_ylabel('ζ')
#        ax.set_title(f'Dislocation Character: {dchar}')
#        ax.legend()
#        plt.show()

def plot_k_data_avg(fname: str) -> None:
    df = pd.read_csv(fname)

    def return_noise_str(k: str) -> str:
        #dic = {'SF-noise': 'SF', 'SS-C-noise': 'SS-C', 'SS-SF-noise':'SS-SF', 'SS-W-noise': 'SS-U', 'W-noise': 'U'}
        dic = {'SF-noise': '\\textbf{SF}', 'SS-C-noise': '\\textbf{SS-C}', 'SS-SF-noise':'\\textbf{SS-SF}', 'SS-W-noise': '\\textbf{SS-U}', 'W-noise': '\\textbf{U}'}
        return dic[k]

    dchars = df['dislocationCharacter'].unique()
    noise_types = df['noiseType'].unique()
    for dchar in dchars:
        fig, ax = plt.subplots(figsize=(11, 8))
        for noise in noise_types:
            data = df[(df['dislocationCharacter'] == dchar) & (df['noiseType'] == noise)]

            # Set x-axis ticks to only show the plotted values
            x_values = data['alloy'].str.extract(r'AlMg(\d+)').astype(int).squeeze()
            zeta = -(data['kExponent']+1) /2
            y_values = zeta

            # Group by unique x_values and average y_values
            grouped = data.groupby(x_values)
            avg_zeta = grouped.apply(lambda x: -(x['kExponent'] + 1).mean() / 2)

            # Plot the averaged values
            #ax.plot(avg_zeta.index, avg_zeta.values, 'o-', linewidth=3, markersize=13, label=noise)
            ax.plot(avg_zeta.index, avg_zeta.values, 'o-', linewidth=3, markersize=13, label=return_noise_str(noise))

        #avg_temp_data = np.mean(temp_data, axis=1)
        #print(avg_temp_data)

        #ax.scatter(x_values, y_values, s=50, label=noise)

        ax.set_xticks(x_values.squeeze().unique())  # Ensures no extra ticks
        ax.set(title=f'',
            xlabel='\\textbf{Mg at-\\%}',
            ylabel='\\textbf{Roughness}')

        # Place legend outside the plot (right side)
        ax.legend(
            frameon=False,
            bbox_to_anchor=(1.05, 1),  # Coordinates (x, y) for legend position
            loc='upper left',           # Anchor point for the legend
            borderaxespad=0.,            # Padding between legend and axes
            prop={'weight': 'bold'}
        )
        plt.tight_layout()
        save_path = './exponent_plot'
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(f"{save_path}/{fname}_{dchar}_roughness_avg.png", bbox_inches='tight', transparent=True)
        plt.close()

def plot_k_data_sep(fname: str) -> None:
    #df = pd.read_csv('k_power_data_crss.csv')
    df = pd.read_csv(fname)

    dchars = df['dislocationCharacter'].unique()
    noise_types = df['noiseType'].unique()
    for dchar in dchars:
        fig, ax = plt.subplots(figsize=(11, 8))
        for noise in noise_types:
            data = df[(df['dislocationCharacter'] == dchar) & (df['noiseType'] == noise)]

            # Set x-axis ticks to only show the plotted values
            x_values = data['alloy'].str.extract(r'AlMg(\d+)').astype(int)
            zeta = -(data['kExponent']+1) /2
            y_values = zeta

            ax.scatter(x_values, y_values, s=50, label=noise)


        ax.set_xticks(x_values.squeeze().unique())  # Ensures no extra ticks
        ax.set(title=f'{dchar}',
            xlabel='Mg Concentration',
            ylabel='Roughness')

        # Place legend outside the plot (right side)
        ax.legend(
            frameon=False,
            bbox_to_anchor=(1.05, 1),  # Coordinates (x, y) for legend position
            loc='upper left',           # Anchor point for the legend
            borderaxespad=0.            # Padding between legend and axes
        )
        save_path = './exponent_plot'
        plt.tight_layout()
        os.makedirs(save_path, exist_ok=True)
        fig.savefig(f"{save_path}/{fname}_{dchar}_roughness.png", bbox_inches='tight', transparent=True)
        plt.close()

def main():
    #plot_k_data_sep('k_power_data_0MPa.csv')
    plot_k_data_avg('k_power_data_0MPa.csv')

    #plot_k_data_sep('k_power_data_200MPa.csv')
    plot_k_data_avg('k_power_data_200MPa.csv')

if __name__ == "__main__":
    main()
