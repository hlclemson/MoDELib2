import re
import os
import sys
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import linregress
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
    #'LINE_TANGENT_AXIS': 1,
    #'GLIDE_DIR_AXIS': 0,
    'LINE_TANGENT_AXIS': 0,
    'GLIDE_DIR_AXIS': 1,
    'PLANE_NORM_DIR_AXIS': 2,
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

def fit_power_law(k: np.ndarray, sxx: np.ndarray) -> float:
    # Filter out zeros/negative values before log transform
    mask = (k > 0) & (sxx > 0)
    k_filtered = k[mask]
    sxx_filtered = sxx[mask]

    # Log transform the data
    log_k = np.log10(k_filtered)
    log_sxx = np.log10(sxx_filtered)

    # Perform linear regression on the log-transformed data
    slope, intercept, r_value, p_value, std_err = linregress(log_k, log_sxx)

    return slope, intercept


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


def plot_distribution(distance_data: np.ndarray, save_path: str, runID: str) -> None:
    n_cols_plot = 2
    n_rows_plot = int(np.ceil(len(distance_data))/n_cols_plot)

    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(14, 10), dpi=200)  # 2x2 grid (adjust to 4x4 if needed)

    # Flatten axes for easy iteration
    axes = axes.ravel()

    # Plot each loop's data
    for ax, partial_sep_dist in zip(axes, distance_data):
        mean_val = np.mean(partial_sep_dist)
        std_val = np.std(partial_sep_dist)
        #ax.plot(runIDs, en, 'o-')
        ax.hist(partial_sep_dist, bins='auto', density=True, alpha=0.3)
        # Add mean line (dotted red line)
        ax.axvline(mean_val, color='red', linestyle=':', linewidth=2)

        # Add annotations (upper right corner)
        stats_text = f'$\\mu$ = {mean_val:.2f} \\AA \n $\\sigma$ = {std_val:.2f} \\AA'
        ax.annotate(stats_text, xy=(0.95, 0.95), xycoords='axes fraction',
                    ha='right', va='top', bbox=None)

        #ax.set_title(f'Loop ID {loop_id}', fontsize=20)
        ax.set_xlabel('separation [\\AA]')
        ax.set_ylabel('Probabilty Density')
        #ax.grid(True, linestyle='--', alpha=0.6)
        #ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))  # Scientific notation

    os.makedirs(save_path, exist_ok=True)
    plt.savefig(Path(save_path)/f"{runID}.png", bbox_inches='tight')
    plt.close()


def plot_power_spectrum_one(k_exponent: float, intercept: float, k: np.ndarray, sxx: np.ndarray, disloc_length: float, pname: str, save_path=None) -> None:
    """Vectorized plotting function for a single loop."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Filter out zeros/negative values before log transform
    mask = (k > 0)
    k_filtered = k[mask]
    #sxx_filtered = sxx[mask]

    ## Log transform the data
    log_k = np.log10(k_filtered)
    #log_sxx = np.log10(sxx_filtered)

    ## Perform linear regression on the log-transformed data
    #slope, intercept, r_value, p_value, std_err = linregress(log_k, log_sxx)

    # Generate trendline points (in log space)
    trendline = k_exponent * log_k + intercept

    # Plot the trendline
    ax.plot(k_filtered, 10**trendline, "b--", alpha=0.6)  # Convert the trendline back to linear scale

    # limit
    ax.axvline(x=2*np.pi/disloc_length, color='red', linestyle='--', linewidth=1)

    # Calculate the equation of the trendline
    #equation = f'P(k) = {slope:.6f}k + {intercept:.6f}'
    equation = f'$P(k) \\propto k^{{{k_exponent:.2f}}}$'
    #ax.text(0.75, 0.95, equation, transform=ax.transAxes, 
    #    fontsize=26, va='top')
    ax.text(0.45, 0.95, equation, transform=ax.transAxes, va='top')

    # Add the equation as an annotation on the plot
    #ax.annotate(equation, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, ha='left', va='top')
    ax.scatter(k, sxx, s=25, c='k', alpha=1,
               edgecolors='none', label=f'')
   
    # Log-scale settings
    ax.set(xscale='log', yscale='log',
           title=f'',
           xlabel='k/b',
           ylabel='Power Spectral Density')

    # Add scientific notation for tick labels
    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    fig.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        #plt.savefig(f"{save_path}/average_psd.png", bbox_inches='tight', transparent=True)
        plt.savefig(f"{save_path}/{pname}", bbox_inches='tight', transparent=True)
    plt.close()


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


def main():

    alloy = "AlMg5"
    #alloy_mat_file = f"{alloy}1350K"
    alloy_mat_file = f"{alloy}"
    #dataDir = f"../xin_SF_noise_{alloy}/generatedData/seed1/"
    #dataDir = f"../../PSD_2025_Aug10_Mg5_E_LT1.5_0MPa_300/generatedData/"
    dataType = "full_french_AlMg50_E_0MPa"
    dataDir = Path("../")/f"{dataType}"/"generatedData"
    #output_dir = Path(dataDir)/"separation_distance"
    fig_dir = Path("figures")/f"{dataType}"
    os.makedirs(fig_dir, exist_ok=True)

    seed_dirs = []
    with os.scandir(dataDir) as dirs:
        for i in dirs:
            if i.is_dir() and "seed" in i.name:
                seed_dirs.append(i.name)
    seed_dirs = sorted(seed_dirs, key=lambda x: int(x.replace('seed', '')))

    for seed_dir in seed_dirs:
        node_dirs = []
        with os.scandir(dataDir/seed_dir) as dirs:
            for i in dirs:
                if i.is_dir() and "node" in i.name:
                    node_dirs.append(i.name)
        node_dirs = sorted(node_dirs, key=lambda x: int(x.replace('node', '')))
        for node_dir in node_dirs:
            str_dirs = []
            with os.scandir(dataDir/seed_dir/node_dir) as dirs:
                for i in dirs:
                    if i.is_dir() and "MPa" in i.name:
                        str_dirs.append(i.name)
            str_dirs = sorted(str_dirs, key=lambda x: int(x.replace('MPa', '')))
            tested_strs = np.array([ int(x.replace('MPa', '')) for x in str_dirs ])
            for str_dir in str_dirs:
                data_dir_temp = dataDir/seed_dir/node_dir/str_dir

                f,fLabels=readFfile(str(data_dir_temp/'F'))
                runIDs=getFarray(f,fLabels,'runID').astype(int)

                wave_vectors_avg = defaultdict(list)
                power_spectra_avg = defaultdict(list)
                for runID in runIDs:
                    ## Group nodes by loop ID
                    filtered_nodes = find_glissile_nodes(data_dir_temp, runID)

                    unique_loop_ids, indices = np.unique(filtered_nodes[:, get_idx('LOOP_ID_COL')], return_inverse=True)
                    loopNodesPerLoopID = {
                        int(loop_id): filtered_nodes[indices == i, get_idx('POS_COLS')]  # Extract network nodes on each loop ID
                        for i, loop_id in enumerate(unique_loop_ids)
                    }

                    # Vectorized center of mass calculation
                    loop_ids = np.array(list(loopNodesPerLoopID.keys()))
                    #node_positions = np.array(list(loopNodesPerLoopID.values()))  # Shape: (n_loops, n_nodes, 3)
                    node_positions = list(loopNodesPerLoopID.values())  # Shape: (n_loops, n_nodes, 3)
                    centers_of_mass = np.array([np.mean(pos, axis=0) for pos in node_positions])
                    center_mass_per_loopID = dict(zip(loop_ids, centers_of_mass))

                    # get positional index for node parsing
                    # 0, 1, 2 are x y z indices, they are just named as line tangent and glide direction
                    lineTangentIdx = get_idx("LINE_TANGENT_AXIS")
                    glideDirIdx = get_idx("GLIDE_DIR_AXIS")

                    # Stack all node positions and centers for batch processing
                    all_nodes = np.concatenate(list(loopNodesPerLoopID.values()))  # Shape: (total_nodes, 3)
                    all_centers = np.repeat(centers_of_mass, [len(pos) for pos in node_positions], axis=0)
                    # Precompute all deviations in the glide direction
                    deviations = np.abs(all_nodes[:, glideDirIdx] - all_centers[:, glideDirIdx])

                    # Split deviations back per loop ID using segment lengths
                    split_indices = np.cumsum([len(pos) for pos in node_positions])[:-1]
                    deviation_per_loop = np.split(deviations, split_indices)

                    # Vectorized FFT calculations
                    for loop_id, dev in zip(loop_ids, deviation_per_loop):
                        n = len(dev) # data number
                        if n == 0: # skip the empty data
                            continue

                        # Get dislocation length for this loop
                        loop_nodes = loopNodesPerLoopID[loop_id]
                        disloc_length = np.max(loop_nodes[:, lineTangentIdx]) - np.min(loop_nodes[:, lineTangentIdx])

                        # Compute FFT
                        k_deviation = np.fft.fft(dev)
                        dt = disloc_length / n
                        k = np.fft.fftfreq(n, dt)
                        sxx = np.abs(k_deviation)**2 / (n * dt)

                        # append the data
                        wave_vectors_avg[loop_id].append(k)
                        power_spectra_avg[loop_id].append(sxx)

                # average the data (each loop)
                for loop_id in loop_ids:
                    wave_vectors_avg[loop_id] = np.mean(np.array(wave_vectors_avg[loop_id]), axis=0)
                    power_spectra_avg[loop_id] = np.mean(np.array(power_spectra_avg[loop_id]), axis=0)
                    print(wave_vectors_avg)

                # average over all the loops
                #wave_vectors_avg_loop_avg = [ np.array(x) for x in wave_vectors_avg.values() ]
                wave_vectors_avg_loop_avg = np.mean(np.array([v for v in wave_vectors_avg.values()]), axis=0)
                power_spectra_avg_loop_avg = np.mean(np.array([v for v in power_spectra_avg.values()]), axis=0)

                # fit power law curve
                k_exponent, intercept = fit_power_law(wave_vectors_avg_loop_avg, power_spectra_avg_loop_avg)

                print(f"disloc_length : {disloc_length}")
                plot_power_spectrum_one(k_exponent, intercept, wave_vectors_avg_loop_avg, power_spectra_avg_loop_avg, disloc_length,
                                f"PSD_{seed_dir}_{node_dir}_{str_dir}.png", save_path=fig_dir)

                    ## find the partial pair (on the same glide plane)
                    #plane_height_per_loopid = {}
                    #for loopID, node_pos in loopNodesPerLoopID.items():
                    #    plane_height_per_loopid[loopID] = np.mean(node_pos[:, get_idx('PLANE_NORM_DIR_AXIS')], dtype=int)

                    #loopid_pair_on_plane = defaultdict(list)
                    #for loopID, height in plane_height_per_loopid.items():
                    #    loopid_pair_on_plane[height].append(loopID)

                    ## find the node pair (on the same line tangent direction coordinate)
                    #partial_dists_all_planes = []
                    #for height, pair_loop_ids in loopid_pair_on_plane.items():
                    #    #partial_pair_pos = []
                    #    loop1, loop2 = pair_loop_ids
                    #    node_pos_loop1 = loopNodesPerLoopID[loop1]
                    #    node_pos_loop2 = loopNodesPerLoopID[loop2]

                    #    partial_dists = np.abs(node_pos_loop1[:, get_idx('GLIDE_DIR_AXIS')] - node_pos_loop2[:, get_idx('GLIDE_DIR_AXIS')])
                    #    partial_dists_all_planes.append(partial_dists)
                    #partial_dists_all_planes = np.array(partial_dists_all_planes)

                    #matDir = str(Path(dataDir)/'inputFiles') # extract material info
                    #b_SI = readValFromMaterialFile(matDir, alloy_mat_file, "b_SI")
                    #to_Ang = 1e-10

                    # convert distance unit from DDD to Angstrom
                    #partial_dists_all_planes *= b_SI/to_Ang
                    #plot_distribution(partial_dists_all_planes, output_dir, runID)

    return 0


if __name__ == "__main__":
    main()
