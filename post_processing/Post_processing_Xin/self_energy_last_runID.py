import re
import os
import sys
import string
from matplotlib.artist import get
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#sys.path.append("../../../build/tools/pyMoDELib/")
#import pyMoDELib

#from typing import Optional
from typing import Union
from typing import DefaultDict
from pathlib import Path
from matplotlib import rcParams
from collections import defaultdict

#sys.path.append("../../../python")
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
    'LINE_TANGENT_AXIS': 1,
    'GLIDE_DIR_AXIS': 0,
    'PLANE_NORM_DIR_AXIS': 2,
    'LOOP_ID_COL': 0,
    'LOOP_TYPE_COL': 11,
    'NETWORK_NODE_COL': 5,
    'SESSILE_TYPE': 1,
    'CORE_E_IDX': 31,
    'ELASTIC_E_IDX': 32,
    'POS_COLS': slice(2, 5),
    'AUX_SOURCE_SINK': slice(0, 2)
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


def plot_data(runIDs: np.ndarray, energy_dict: dict, title: str, fig_fname: str) -> None:
    n_cols_plot = 2
    n_rows_plot = int(np.ceil(len(energy_dict.keys()))/n_cols_plot)
    # Create figure and subplots
    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(12, 10), dpi=200)  # 2x2 grid (adjust to 4x4 if needed)
    fig.suptitle(title, fontsize=20, y=1.00)

    # Flatten axes for easy iteration (if using 2x2)
    axes = axes.ravel()

    # Plot each loop's data
    for ax, (loop_id, energies) in zip(axes, energy_dict.items()):
        ax.plot(runIDs, energies, 'o-')
        ax.set_title(f'Loop ID {loop_id}', fontsize=20)
        ax.set_xlabel('run ID', fontsize=16)
        ax.set_ylabel('Energy (J)', fontsize=16)
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))  # Scientific notation

    # Adjust layout
    plt.tight_layout()
    fig.savefig(fig_fname)
    plt.close()


def plot_dislocation_loop(loop_pos_per_loopID: dict, save_path: str, runID: str) -> None:
    n_cols_plot = 2
    n_rows_plot = int(np.ceil(len(loop_pos_per_loopID.keys()))/n_cols_plot)

    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(12, 10), dpi=200)  # 2x2 grid (adjust to 4x4 if needed)

    # Flatten axes for easy iteration
    axes = axes.ravel()

    # Plot each loop's data
    for ax, (loop_id, node_pos) in zip(axes, loop_pos_per_loopID.items()):
        # Vectorized scatter plot with optimized markers
        xNodePos = node_pos[:, 0]
        yNodePos = node_pos[:, 1]
        ax.scatter(xNodePos, yNodePos, s=30, c='k', alpha=0.7, 
                edgecolors='none', label=f'Loop {loop_id}')
        ax.set(title=f'loop {loop_id}')

    # Add scientific notation for tick labels
    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    os.makedirs(save_path, exist_ok=True)
    plt.savefig(Path(save_path)/f"{runID}.png", bbox_inches='tight')
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

    alloy = "AlMg15"
    #alloy_mat_file = f"{alloy}1350K"
    alloy_mat_file = f"{alloy}"

    #dataDir = f"../xin_SF_noise_{alloy}/generatedData/seed{seed}/"
    dataDir = f"../xin_SF_noise_{alloy}/generatedData/"
    output_dir = f"../xin_SF_noise_{alloy}/dislocation_self_energy"

    seed_dirs = []
    with os.scandir(dataDir) as dirs:
        for i in dirs:
            if i.is_dir() and "seed" in i.name:
                seed_dirs.append(i.name)
    seed_dirs = sorted(seed_dirs, key=lambda x: int(x.replace('seed', '')))

    elastic_data_all = []
    core_data_all = []
    for seed_dir in seed_dirs:
        data_dir_temp = Path(dataDir)/seed_dir
        f,fLabels=readFfile(str(data_dir_temp/'F'))
        runIDs=getFarray(f,fLabels,'runID').astype(int)
        runID = runIDs[-1]

        t_elastic_energies = defaultdict(list)
        t_core_energies = defaultdict(list)

        filtered_nodes = find_glissile_nodes(data_dir_temp, runID)
        unique_loop_ids, indices = np.unique(filtered_nodes[:, get_idx('LOOP_ID_COL')], return_inverse=True)

        loop_nodes_per_loopID = {
            int(loop_id): filtered_nodes[indices == i, get_idx('NETWORK_NODE_COL')]  # Extract network nodes on each loop ID
            for i, loop_id in enumerate(unique_loop_ids)
        }

        loop_pos_per_loopID = {
            int(loop_id): filtered_nodes[indices == i, get_idx('POS_COLS')]  # Extract network nodes on each loop ID
            for i, loop_id in enumerate(unique_loop_ids)
        }

        aux =readAUXtxt(str(data_dir_temp/'evl'/f'ddAux_{runID}'))
        gauss_data = np.array(aux.gaussPoints)

        # Convert to NumPy arrays
        source_sink_pairs = gauss_data[:, get_idx('AUX_SOURCE_SINK')]
        # init dictionary
        gauss_elastic_per_loopID = {loopID: [] for loopID in loop_nodes_per_loopID}
        gauss_core_per_loopID = {loopID: [] for loopID in loop_nodes_per_loopID}
        for loopID, netNodes in loop_nodes_per_loopID.items():
            # create lookup hash map out of network nodes
            node_set = set(netNodes)
            # create mask boolean list for loop ID
            mask = np.isin(source_sink_pairs, list(node_set)).all(axis=1)
            gauss_elastic_per_loopID[loopID] = gauss_data[mask, get_idx('ELASTIC_E_IDX')]
            gauss_core_per_loopID[loopID] = gauss_data[mask, get_idx('CORE_E_IDX')]

        matDir = str(data_dir_temp/'inputFiles') # extract material info
        b_SI = readValFromMaterialFile(matDir, alloy_mat_file, "b_SI")
        mu0_SI = readValFromMaterialFile(matDir, alloy_mat_file, "mu0_SI")
        isf_SI = readValFromMaterialFile(matDir, alloy_mat_file, "ISF_SI")

        # unit conversion
        mub3_to_Nm = mu0_SI * b_SI**3
        loop_ids = np.array(list(loop_nodes_per_loopID.keys()))

        elastic_sums = np.array([np.sum(gauss_elastic_per_loopID[loop_id]) for loop_id in loop_ids]) * mub3_to_Nm
        core_sums = np.array([np.sum(gauss_core_per_loopID[loop_id]) for loop_id in loop_ids]) * mub3_to_Nm

        # plot dislocation lines on each output
        #plot_dislocation_loop(loop_pos_per_loopID, Path(output_dir)/'disloc_line', runID)
        for loop_id, elastic_e, core_e in zip(loop_ids, elastic_sums, core_sums):
            t_elastic_energies[loop_id].append(elastic_e)
            t_core_energies[loop_id].append(core_e)

        # find the partial pair (on the same glide plane)
        plane_height_per_loopid = {}
        for loopID, node_pos in loop_pos_per_loopID.items():
            plane_height_per_loopid[loopID] = np.mean(node_pos[:, get_idx('PLANE_NORM_DIR_AXIS')], dtype=int)

        # find the loop ids on each glide plane dipole
        loopid_pair_on_plane = defaultdict(list)
        for loopID, height in plane_height_per_loopid.items():
            loopid_pair_on_plane[height].append(loopID)

        # find the node pair (on the same line tangent direction coordinate)
        partial_dists_all_planes = []
        for height, pair_loop_ids in loopid_pair_on_plane.items():
            #partial_pair_pos = []
            loop1, loop2 = pair_loop_ids
            node_pos_loop1 = loop_pos_per_loopID[loop1]
            node_pos_loop2 = loop_pos_per_loopID[loop2]

            partial_dists = np.abs(node_pos_loop1[:, get_idx('GLIDE_DIR_AXIS')] - node_pos_loop2[:, get_idx('GLIDE_DIR_AXIS')])
            partial_dists_all_planes.append(partial_dists)
        partial_dists_all_planes = np.array(partial_dists_all_planes)
        # convert distance unit from DDD to m
        partial_dists_all_planes *= b_SI
        # calculate the mean separation distance of each plane
        mean_sep_all_planes = np.mean(partial_dists_all_planes, axis=1)
        # calculate the standard deviation of each plane
        std_sep_all_planes = np.std(partial_dists_all_planes, axis=1)
        # calculate the average of the mean separations of each plane
        avg_mean_sep = np.mean(mean_sep_all_planes)
        # calculate the avg of the std of each plane
        avg_std_sep = np.mean(std_sep_all_planes)

        # average the energy over loops
        avg_over_loops_elastic_e = np.mean(np.array(list(t_elastic_energies.values())))
        avg_over_loops_core_e = np.mean(np.array(list(t_core_energies.values())))

        elastic_data = [alloy, seed_dir.replace('seed', ''), b_SI, mu0_SI, isf_SI, avg_mean_sep, avg_std_sep, avg_over_loops_elastic_e]
        core_data = [alloy, seed_dir.replace('seed', ''), b_SI, mu0_SI, isf_SI, avg_mean_sep, avg_std_sep, avg_over_loops_core_e]
        elastic_data_all.append(elastic_data)
        core_data_all.append(core_data)

    #elastic_data_all = np.array(elastic_data_all)
    #print(elastic_data_all.shape)

    # Convert dictionaries to DataFrames
    header_core = ["alloy", "seed", "b [m]", "mu [Pa]", "ISF [J/m2]", "mean separation [m]", "std separation [m]", "core energy [J]"]
    header_elastic = ["alloy", "seed", "b [m]", "mu [Pa]", "ISF [J/m2]", "mean separation [m]", "std separation [m]", "elastic energy [J]"]
    df_core = pd.DataFrame(core_data_all, columns=header_core)
    df_elastic = pd.DataFrame(elastic_data_all, columns=header_elastic)

    # save data
    os.makedirs(output_dir, exist_ok=True)
    df_elastic.to_csv(Path(output_dir)/f'converged_elastic_energies_{alloy}.csv', index=False)
    df_core.to_csv(Path(output_dir)/f'converged_core_energies_{alloy}.csv', index=False)

if __name__ == "__main__":
    main()
