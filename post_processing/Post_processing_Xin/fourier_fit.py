import re
import os
import sys
import string

from matplotlib.artist import get
import scipy as sp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#sys.path.append("../../build/tools/pyMoDELib/")
#import pyMoDELib

from pathlib import Path
from typing import Union
from typing import DefaultDict
from matplotlib import rcParams
from collections import defaultdict

#sys.path.append("../../python")
#from modlibUtils import *


# Configure global plot settings (applies to all figures)
rcParams.update({
    'figure.dpi': 100,
    'figure.autolayout': True,  # Prevent label clipping
    'axes.grid': True,
    'grid.alpha': 0.3,
    'text.usetex': True,   # Enable LaTeX
    'font.size': 20,     # Default font size for text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    'text.latex.preamble': r'\usepackage{amsmath}'  # Additional packages
})


class EVL:
    nodes=np.empty([0,0])


class AUX:
    meshNodes=np.empty([0,0])
    gaussPoints=np.empty([0,0])
    periodicPatches=np.empty([0,0])


def get_idx(idx_name: str) -> Union[int, slice]:
    """ returns constant index """
    index_mapping = {
    'LOOP_ID_COL': 0,
    'LOOP_TYPE_COL': 11,
    'NETWORK_NODE_COL': 5,
    'SESSILE_TYPE': 1,
    'CORE_E_IDX': 31,
    'ELASTIC_E_IDX': 32,
    'POS_COLS': slice(2, 5) 
    }

    if idx_name not in index_mapping:
        raise KeyError(f"invalid index name: '{idx_name}' Options are: {list(index_mapping.keys())}");

    return index_mapping.get(idx_name, None)


def readAUXtxt(filename):
    auxFile = open(filename+'.txt', "r")
    numNodes=int(auxFile.readline().rstrip())
    numGPs=int(auxFile.readline().rstrip())
    numPGPP=int(auxFile.readline().rstrip())
    aux=AUX();
    aux.gaussPoints=np.empty([numGPs, 41])
    for k in np.arange(numNodes):
        np.meshNodes=np.fromstring(auxFile.readline().rstrip(), sep=' ')
    for k in np.arange(numGPs):
        aux.gaussPoints[k,:]=np.fromstring(auxFile.readline().rstrip(), sep=' ')
    return aux


def readEVLtxtLoopNode(filename: str) -> EVL:
    evlFile = open(filename+'.txt', "r")
    numNetNodes=int(evlFile.readline().rstrip())
    numLoops=int(evlFile.readline().rstrip())
    numLinks=int(evlFile.readline().rstrip())
    numLoopNodes=int(evlFile.readline().rstrip())
    numSpInc=int(evlFile.readline().rstrip())
    numPolyInc=int(evlFile.readline().rstrip())
    numPolyIncNodes=int(evlFile.readline().rstrip())
    numPolyIncEdges=int(evlFile.readline().rstrip())
    numEDrow=int(evlFile.readline().rstrip())
    numCDrow=int(evlFile.readline().rstrip())
    evl=EVL();

    dislocLoopDataSpan = 17
    evl.dislocLoop=np.empty([numLoops, dislocLoopDataSpan], dtype=np.float64)
    # Read file line by line until we fill the array
    for n in range(numLoops):
        while True:  # Keep reading until we get a valid line
            line = evlFile.readline()
            # End of file check
            if not line:
                raise ValueError(f"Reached EOF before reading {numLoops} loops")
            data = np.fromstring(line.rstrip(), sep=' ')
            if len(data) == dislocLoopDataSpan:
                evl.dislocLoop[n] = data
                # Move to next node after successful read 
                break

    # init connectivity data
    loopNodeDataSpan = 11
    evl.loopNodes=np.empty([numLoopNodes, loopNodeDataSpan], dtype=np.float64)
    # Read file line by line until we fill the array
    for m in range(numLoopNodes):
        while True:  # Keep reading until we get a valid line
            line = evlFile.readline()
            # End of file check
            if not line:
                raise ValueError(f"Reached EOF before reading {numLoopNodes} nodes")
            data = np.fromstring(line.rstrip(), sep=' ')
            if len(data) == loopNodeDataSpan:
                evl.loopNodes[m] = data
                # Move to next node after successful read 
                break

    return evl


def getFarray(F,Flabels,label):
    k=0;
    for line in Flabels:
        if line==label:
            if F.ndim==1:
                return F[k]
            else:
                return F[:,k]
        k=k+1
    return np.zeros(shape=(0,0))


def readFfile(folder):
    F=np.loadtxt(Path(folder)/'F_0.txt');
    with open(Path(folder)/'F_labels.txt', 'r') as f:
        lines = f.readlines()
        for idx in range(len(lines)):
            lines[idx] = lines[idx].rstrip()
    return F,lines


def extract_evl_number(fname: str) -> int:
    """Extract numeric portion from 'evl_XXXXX.txt' filenames"""
    match = re.search(r'evl_(\d+)', fname)
    return int(match.group(1)) if match else 0


def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
    #with open(f"{matDir}/{alloy}.txt", "r") as mFile:
    with open(Path(matDir)/f'{alloy}.txt', "r") as mFile:
        for line in mFile:
            # strip down comments from the data
            line = line.split(";")[0]
            if line.startswith(f"{var}"):
                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
                value = float(line.split("=")[1].strip())
                break
    return value


def fourier_series(t, amplitudes, phases, freqs, n_terms, a0=0):
    reconstructed = np.zeros_like(t, dtype=float)
    #reconstructed = a0 * np.ones_like(t)
    for i in range(n_terms):
        omega = 2 * np.pi * freqs[i]
        #reconstructed += amplitudes[i] * np.cos(omega * t + phases[i])
        reconstructed += amplitudes[i] * np.sin(omega * t + phases[i])
    reconstructed -= a0
    return reconstructed


def find_glissile_nodes(dataDir, runID):
    evl = readEVLtxtLoopNode(str(Path(dataDir)/'evl'/f'evl_{runID}'))
    #evl = readEVLtxtLoopNode(str(Path(dataDir)/'evl'/f'{evl_file}'))
    loop_data = evl.dislocLoop
    loop_nodes = evl.loopNodes

    # Get sessile loop IDs (vectorized)
    sessile_mask = loop_data[:, get_idx('LOOP_TYPE_COL')] == get_idx('SESSILE_TYPE')
    sessile_loop_ids = loop_data[sessile_mask, get_idx('LOOP_ID_COL')].astype(int)

    # it seems there is a bug in MoDELib2 not marking sessile loops accurately
    # if sessile loop is not found
    if not len(sessile_loop_ids):
        unique_loop_ids, node_count_per_loop = np.unique(loop_nodes[:, get_idx('LOOP_ID_COL')], return_counts=True)
        all_nodes_per_loop = { int(k): v for k, v in zip(unique_loop_ids, node_count_per_loop) }
        # hardcoded threshold connecting node number, they are usually 3+3+3 of them, so 10
        # if you use less than 10 nodes on your glissile dislocations, it breaks the code
        sessile_loop_ids = [ key for key, count in all_nodes_per_loop.items() if count < 10 ]

    # Get all unique loop IDs (vectorized)
    all_loop_ids = np.unique(loop_data[:, get_idx('LOOP_ID_COL')]).astype(int)

    # Filter non-sessile loop IDs (vectorized)
    non_sessile_loop_ids = all_loop_ids[~np.isin(all_loop_ids, sessile_loop_ids)]

    # Get sessile network nodes (vectorized)
    sessile_network_mask = np.isin(loop_nodes[:, get_idx('LOOP_ID_COL')], sessile_loop_ids)
    sessile_network_nodes = np.unique(loop_nodes[sessile_network_mask, get_idx('NETWORK_NODE_COL')]).astype(int)

    # Precompute masks for fast filtering
    non_sessile_node_mask = ~np.isin(loop_nodes[:, get_idx('NETWORK_NODE_COL')], sessile_network_nodes)
    valid_loop_mask = np.isin(loop_nodes[:, get_idx('LOOP_ID_COL')], non_sessile_loop_ids)
    combined_mask = non_sessile_node_mask & valid_loop_mask

    # Group nodes by loop ID
    filtered_nodes = loop_nodes[combined_mask]
    return filtered_nodes



def plot_fourier_approx(loopNode_pos_per_loopID: dict, save_path: str, runID: str) -> None:
    n_cols_plot = 2
    n_rows_plot = int(np.ceil(len(loopNode_pos_per_loopID.keys()))/n_cols_plot)

    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(14, 10), dpi=200)

    # Flatten axes for easy iteration
    axes = axes.ravel()

    # Plot each loop's data
    for ax, (loop_id, node_pos) in zip(axes, loopNode_pos_per_loopID.items()):
        signal = node_pos[:, 0]  # Signal (1st column)
        time = node_pos[:, 1]  # Time (2nd column)

        N = len(signal)  # Number of samples
        sampling_period = time[1] - time[0]  # Time step (assuming uniform sampling)
        sampling_freq = 1 / sampling_period  # Sampling frequency

        # Compute FFT
        signal_fft = np.fft.fft(signal)
        frequencies = np.fft.fftfreq(N, d=sampling_period)  # Frequency bins
        # Take only positive frequencies (Nyquist theorem)
        positive_freq = frequencies[:N//2]
        amplitudes = 2/N * np.abs(signal_fft[:N//2])  # Normalized amplitudes
        phases = np.angle(signal_fft[:N//2]) + np.pi/2   # Phases (in radians)

        n_terms = len(signal)//2  # Adjust based on how many harmonics to include
        a0 = np.real(signal_fft[0]) / N
        signal_reconstructed = fourier_series(time, amplitudes, phases, positive_freq, n_terms, a0)

        ax.set(title=f'loop {loop_id}')
        ax.plot(time, signal, 'b-', label='input')
        ax.plot(time, signal_reconstructed, 'r--', label=f'Fourier Fit ({n_terms} terms)')

    #plt.legend()
    #plt.show()
    os.makedirs(save_path, exist_ok=True)
    plt.savefig(Path(save_path)/f"{runID}.png", bbox_inches='tight')
    plt.close()

def main():

    alloy = "AlMg5"
    #alloy_mat_file = f"{alloy}1350K"
    alloy_mat_file = f"{alloy}"
    dataDir = f"../xin_SF_noise_{alloy}/generatedData/seed1/"
    output_dir = Path(dataDir)/"fourier_approximation"

    f,fLabels=readFfile(str(Path(dataDir)/'F'))
    runIDs=getFarray(f,fLabels,'runID').astype(int)

    # Constants for column indices
    for runID in runIDs:
        filtered_nodes = find_glissile_nodes(dataDir, runID)

        unique_loop_ids, indices = np.unique(filtered_nodes[:, get_idx('LOOP_ID_COL')], return_inverse=True)

        loopNode_pos_per_loopID = {
            int(loop_id): filtered_nodes[indices == i, get_idx('POS_COLS')]  # Extract network nodes on each loop ID
            for i, loop_id in enumerate(unique_loop_ids)
        }

        plot_fourier_approx(loopNode_pos_per_loopID, output_dir, runID)

    return 0


if __name__ == "__main__":
    main()
