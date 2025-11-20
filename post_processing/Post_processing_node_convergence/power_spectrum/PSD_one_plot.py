import os
import re
import csv
import sys
import json
import string
import sqlite3
import itertools
import subprocess
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import linregress
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D

from python.modlibUtils import readEVLtxtPlus
from python.modlibUtils import readEVLtxtLoopNode

#plt.rcParams['text.usetex'] = True
sys.path.append("./python")
from modlibUtils import *
from pathlib import Path
#sys.path.append("../../build/tools/pyMoDELib/")
#import pyMoDELib
from multiprocessing import Pool, Queue, Process, allow_connection_pickling


# Configure global plot settings (applies to all figures)
rcParams.update({
    'figure.dpi': 200,
    'figure.autolayout': True,  # Prevent label clipping
    'axes.grid': True,
    'grid.alpha': 0.6,
    #'text.usetex': True,   # Enable LaTeX
    'text.usetex': False,   # Enable LaTeX
    'font.size': 33,     # Default font size for text
    'mathtext.fontset': 'stix',  # Use STIX font for math text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    #'font.weight': 'bold',          # Base font weight
    #'text.latex.preamble': r'\usepackage{amsmath} \boldmath',  # Additional packages
    #'text.latex.preamble': r'\usepackage{amsmath}',  # Additional packages
})


def plot_psd_in_one(
    appliedStr: str,
    mainDataDir: str,
    disloc_type: str,
    noise_type: str,
    disloc_character: str,
    alloy: str,
    disloc_len: int,
    dislocGlideLineTangentIdx: dict,
    xAxisDataLabel: str,
    yAxisDataLabel: str,
    minimizationSteps: int,
    sql_fname: str,
    pid: int,
) -> None:

    dataTypePath = Path(disloc_type)/noise_type/disloc_character/alloy/f"{alloy}-IMG1-{disloc_character[0].upper()}-{disloc_len}b"

    genDataDir = os.path.join(
        mainDataDir,
        dataTypePath,
        "MoDELib2",
        "tutorials",
        "DislocationDynamics",
        "periodicDomains",
        "generatedData",
    )

    # check if the directory exist or not
    if not os.path.isdir(genDataDir):
        print(f"Directory not found: {genDataDir}")
        return

    dataDirs = []
    # Vectorized version of scanning generated data directory using list comprehension
    dataDirs = [
        entry
        for entry in os.listdir(genDataDir)
        if os.path.isdir(os.path.join(genDataDir, entry))
    ]

    """
    Regex pattern to match the number between 'sID' and 'exID' in the first dir entry
    this is to count how many slip system IDs are specified. 
    Partial dislocation has more than one slip system ID.
    """
    pattern = r"sID(\d+)exID"
    match = re.search(pattern, dataDirs[0])
    if match:
        slipID = match.group(1)
        slipIDnum = [int(id) for id in slipID]

    """
    Find the unique seed numbers
    """
    seeds = np.unique([int(x.split("S")[-1]) for x in dataDirs])

    # detect the number of nodes


def plot_all_PSD(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr: str) -> None:
    for dchar in dislocChars:
        for alloy in alloys:
            used_pairs = set()
            def return_noise_str(k: str) -> str:
                #dic = {'SF-noise': '\\textbf{SF}', 'SS-C-noise': '\\textbf{SS-C}', 'SS-SF-noise':'\\textbf{SS-SF}', 'SS-W-noise': '\\textbf{SS-U}', 'W-noise': '\\textbf{U}'}
                dic = {'SF-noise': 'SF', 'SS-C-noise': 'SS-C', 'SS-SF-noise':'SS-SF', 'SS-W-noise': 'SS-U', 'W-noise': 'U'}

                pair = (k, dic[k])
                if pair in used_pairs:
                    return None
                else:
                    used_pairs.add(pair)
                    return dic[k]

            #fig, ax = plt.subplots(figsize=(14, 8))
            fig, ax = plt.subplots(figsize=(13, 8))

            is_guideline = False

            for dislocT, noiseT in dislocTypes.items():
                for noise in noiseT["noiseTypes"]:
                    for dlen in dislocLens:
                        datDir = Path(mainDataFolder)/dislocT/noise/dchar/alloy/f"{alloy}-IMG1-{dchar[0].upper()}-{dlen}b"
                        wave_vectors = np.load(datDir/f'wave_vectors_avg_{appliedStr}.npy')
                        power_spectra = np.load(datDir/f'power_spectra_avg_{appliedStr}.npy')

                        # Filter out zeros/negative values before log transform
                        #k = wave_vectors[len(wave_vectors)//2:]
                        #sxx = power_spectra[len(wave_vectors)//2:]
                        k = wave_vectors
                        sxx = power_spectra
                        mask = (k > 0)
                        k_filtered = k[mask]
                        sxx_filtered = sxx[mask]

                        # Log transform the data
                        log_k = np.log10(k_filtered)
                        log_sxx = np.log10(sxx_filtered)

                        # Perform linear regression on the log-transformed data
                        slope, intercept, r_value, p_value, std_err = linregress(log_k, log_sxx)
                        #print(intercept)
                        #intercept-=2

                        # Create 1/k^4 guideline using the fitted intercept
                        k_line = np.linspace(min(k_filtered), max(k_filtered), 100)
                        sxx_line = 10**intercept * k_line**(-4)

                        if not is_guideline:
                            # Plot the guideline
                            plt.plot(k_line, sxx_line, 'k--', label=r'$1/k^4$', alpha=0.4)
                            is_guideline = True

                        # Generate trendline points (in log space)
                        #trendline = k_exponent * log_k + intercept

                        # Plot the trendline
                        #ax.plot(k_filtered, 10**trendline, "r--", alpha=0.6)  # Convert the trendline back to linear scale

                        # Calculate the equation of the trendline
                        #equation = f'P(k) = {slope:.6f}k + {intercept:.6f}'
                        #equation = f'$P(k) \\propto k^{{{k_exponent:.2f}}}$'
                        #ax.text(0.75, 0.95, equation, transform=ax.transAxes, 
                        #    fontsize=26, va='top')
                        #ax.text(0.45, 0.95, equation, transform=ax.transAxes, va='top')

                        # Add the equation as an annotation on the plot
                        #ax.annotate(equation, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, ha='left', va='top')
                        #ax.scatter(wave_vectors, power_spectra, s=25, c='k', alpha=1,
                        #        edgecolors='none', label=f'')
                        ax.scatter(wave_vectors, power_spectra, s=60, alpha=0.8,
                                #edgecolors='none', label=f'')
                                edgecolors='none', label=return_noise_str(noise))

                        # Log-scale settings
                        ax.set(xscale='log', yscale='log',
                            title=f'',
                            #xlabel='\\textbf{Wave Vector (k)}',
                            #ylabel='\\textbf{Power Spectral Density}')
                            xlabel='Wave Vector (k)',
                            ylabel='Power Spectral Density')

                        # Add scientific notation for tick labels
                        #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
                        ax.legend(
                        frameon=False,
                        bbox_to_anchor=(1.05, 1),  # Coordinates (x, y) for legend position
                        loc='upper left',           # Anchor point for the legend
                        borderaxespad=0.,            # Padding between legend and axes
                        #prop={'weight': 'bold'},
                        markerscale=2.2
                        )
                        ax.set_xlim(5e-3, 5e-1)
                        ax.set_ylim(5e-5, 5e+1)

            fig.tight_layout()
            plt.savefig(Path(mainDataFolder)/dislocT/f"all_psd_plot_{dchar}_{alloy}_{appliedStr}_g.png", bbox_inches='tight', transparent=True)
            plt.close()


def plot_all_PSD_filtered(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr: str, filter: list) -> None:
    for dchar in dislocChars:
        for alloy in alloys:
            used_pairs = set()
            def return_noise_str(k: str) -> str:
                #dic = {'SF-noise': '\\textbf{SF}', 'SS-C-noise': '\\textbf{SS-C}', 'SS-SF-noise':'\\textbf{SS-SF}', 'SS-W-noise': '\\textbf{SS-U}', 'W-noise': '\\textbf{U}'}
                dic = {'SF-noise': 'SF', 'SS-C-noise': 'SS-C', 'SS-SF-noise':'SS-SF', 'SS-W-noise': 'SS-U', 'W-noise': 'U'}

                pair = (k, dic[k])
                if pair in used_pairs:
                    return None
                else:
                    used_pairs.add(pair)
                    return dic[k]

            #fig, ax = plt.subplots(figsize=(14, 8))
            fig, ax = plt.subplots(figsize=(13, 8))

            is_guideline = False

            for dislocT, noiseT in dislocTypes.items():
                for noise in noiseT["noiseTypes"]:
                    if noise in filter:
                        continue
                    for dlen in dislocLens:
                        datDir = Path(mainDataFolder)/dislocT/noise/dchar/alloy/f"{alloy}-IMG1-{dchar[0].upper()}-{dlen}b"
                        wave_vectors = np.load(datDir/f'wave_vectors_avg_{appliedStr}.npy')
                        power_spectra = np.load(datDir/f'power_spectra_avg_{appliedStr}.npy')

                        # Filter out zeros/negative values before log transform
                        #k = wave_vectors[len(wave_vectors)//2:]
                        #sxx = power_spectra[len(wave_vectors)//2:]
                        k = wave_vectors
                        sxx = power_spectra
                        mask = (k > 0)
                        k_filtered = k[mask]
                        sxx_filtered = sxx[mask]

                        # Log transform the data
                        log_k = np.log10(k_filtered)
                        log_sxx = np.log10(sxx_filtered)

                        # Perform linear regression on the log-transformed data
                        slope, intercept, r_value, p_value, std_err = linregress(log_k, log_sxx)
                        #print(intercept)
                        #intercept-=2

                        # Create 1/k^4 guideline using the fitted intercept
                        k_line = np.linspace(min(k_filtered), max(k_filtered), 100)
                        sxx_line = 10**intercept * k_line**(-4)

                        if not is_guideline:
                            # Plot the guideline
                            plt.plot(k_line, sxx_line, 'k--', label=r'$1/k^4$', alpha=0.4)
                            is_guideline = True

                        # Generate trendline points (in log space)
                        #trendline = k_exponent * log_k + intercept

                        # Plot the trendline
                        #ax.plot(k_filtered, 10**trendline, "r--", alpha=0.6)  # Convert the trendline back to linear scale

                        # Calculate the equation of the trendline
                        #equation = f'P(k) = {slope:.6f}k + {intercept:.6f}'
                        #equation = f'$P(k) \\propto k^{{{k_exponent:.2f}}}$'
                        #ax.text(0.75, 0.95, equation, transform=ax.transAxes, 
                        #    fontsize=26, va='top')
                        #ax.text(0.45, 0.95, equation, transform=ax.transAxes, va='top')

                        # Add the equation as an annotation on the plot
                        #ax.annotate(equation, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, ha='left', va='top')
                        #ax.scatter(wave_vectors, power_spectra, s=25, c='k', alpha=1,
                        #        edgecolors='none', label=f'')
                        ax.scatter(wave_vectors, power_spectra, s=60, alpha=0.8,
                                #edgecolors='none', label=f'')
                                edgecolors='none', label=return_noise_str(noise))

                        # Log-scale settings
                        ax.set(xscale='log', yscale='log',
                            title=f'',
                            #xlabel='\\textbf{Wave Vector (k)}',
                            #ylabel='\\textbf{Power Spectral Density}')
                            xlabel='Wave Vector (k)',
                            ylabel='Power Spectral Density')

                        # Add scientific notation for tick labels
                        #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
                        ax.legend(
                        frameon=False,
                        bbox_to_anchor=(1.05, 1),  # Coordinates (x, y) for legend position
                        loc='upper left',           # Anchor point for the legend
                        borderaxespad=0.,            # Padding between legend and axes
                        #prop={'weight': 'bold'},
                        markerscale=2.2
                        )
                        #ax.set_xlim(2e-4, 1)
                        #ax.set_ylim(1e-1, 20)
                        ax.set_xlim(5e-3, 5e-1)
                        ax.set_ylim(5e-5, 6)

            fig.tight_layout()
            plt.savefig(Path(mainDataFolder)/dislocT/f"sf_w_psd_plot_{dchar}_{alloy}_{appliedStr}_g.png", bbox_inches='tight', transparent=True)
            plt.close()


def main() -> None:
    mainDataFolder = Path("./Powerspectrum/")

    with open("config.json", "r") as f:
        config = json.load(f)
        xAxisDataLabel = config["xAxisData"]
        yAxisDataLabel = config["yAxisData"]
        minimizationSteps = config["minimizationSteps"]
        dislocTypes = config["dislocTypes"]
        dislocChars = config["dislocCharacters"]
        dislocGlideLineTangentIdx = config["dislocGlideLineTangentIdx"]
        alloys = config["alloys"]
        dislocLens = config["dislocLengths"]

    #figDataDir = f"{mainDataFolder}/Figure/"
    figDataDir = f"./Figure/"

    # Using itertools.product to generate all combinations
    dataPathComponents = itertools.product(
        dislocTypes.items(), dislocChars, alloys, dislocLens
    )

    appliedStr1 = "0MPa"
    appliedStr2 = "200MPa"

    plot_all_PSD(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr1)
    plot_all_PSD(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr2)
    plot_all_PSD_filtered(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr1, ["SS-C-noise", "SS-SF-noise","SS-W-noise"])
    plot_all_PSD_filtered(mainDataFolder, dislocTypes, dislocChars, alloys, dislocLens, appliedStr2, ["SS-C-noise", "SS-SF-noise","SS-W-noise"])


if __name__ == "__main__":
    main()
