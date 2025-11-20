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
import config

#from typing import Optional
from typing import Union
from typing import DefaultDict
from pathlib import Path
from matplotlib import rcParams
from collections import defaultdict

from scipy.sparse import data

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
    #nodes2=np.empty([0,0])


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
    'NODE_ID_COL_IDX_IN_NODE': 0,
    'VELOCITY_COLS_IDX_IN_NODE': slice(4, 6),
    'POS_COLS': slice(2, 5),
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


def readEVLtxt(filename):
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
    evl.nodes=np.empty([numNetNodes, 7])
    for k in np.arange(numNetNodes):
        data=np.fromstring(evlFile.readline().rstrip(), sep=' ')
        #evl.nodes[k,:]=data[1:4]
        evl.nodes[k,:]=data[0:7]
    return evl


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

    #alloy = "FeCr8"
    ##alloy_mat_file = f"{alloy}1350K"
    ##alloy_mat_file = f"{alloy}"
    #alloy_mat_file = f"FeCrAl_Fe"
    ##dataDir = f"../mobility_Fe_FeCr8_edge/generatedData/seed{seed}/"
    #dataDir = f"../mobility_Fe_FeCr8_edge/generatedData/"
    #output_dir = Path(dataDir)/".."/"mobility"

    alloy = config.alloy
    alloy_mat_file = config.alloy_mat_file
    dataDir = config.dataDir
    output_dir = config.output_dir

    seed_data_dir = []
    with os.scandir(dataDir) as seed_dir:
        for i in seed_dir:
            if i.is_dir() and "seed" in i.name:
                seed_data_dir.append(i.name)
    seed_data_dir = sorted(seed_data_dir, key=lambda x: int(x.replace('seed', '')))

    record_data_output = []
    for seed in seed_data_dir:

        stress_data_dir = []
        with os.scandir(Path(dataDir)/seed) as dirs:
            for i in dirs:
                if i.is_dir() and "MPa" in i.name:
                    stress_data_dir.append(i.name)
        stress_data_dir = sorted(stress_data_dir, key=lambda x: int(x.replace('MPa', '')))

        for stress_dir in stress_data_dir:
            ddd_data_dir = Path(dataDir)/seed/stress_dir

            # extract material info
            b_SI = readValFromMaterialFile(ddd_data_dir/"inputFiles", alloy_mat_file, "b_SI")
            mu0_SI = readValFromMaterialFile(ddd_data_dir/"inputFiles", alloy_mat_file, "mu0_SI")
            rho_SI = readValFromMaterialFile(ddd_data_dir/"inputFiles", alloy_mat_file, "rho_SI")
            cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
            convertTimeUnit = b_SI / cs  # [sec]

            f,fLabels=readFfile(str(ddd_data_dir/'F'))
            runIDs=getFarray(f,fLabels,'runID').astype(int)

            t_elastic_energies = defaultdict(list)
            t_core_energies = defaultdict(list)

            record_data_mean_temp = []
            #runIDs=runIDs[20:]
            half = len(runIDs) // 2  # Integer division (Python 3)
            runIDs = runIDs[half:]
            for runID in runIDs:
            #runID = runIDs[-1]

                filtered_nodes = find_glissile_nodes(ddd_data_dir, runID)
                unique_loop_ids, indices = np.unique(filtered_nodes[:, get_idx('LOOP_ID_COL')], return_inverse=True)

                loop_nodes_per_loopID = {
                    int(loop_id): filtered_nodes[indices == i, get_idx('NETWORK_NODE_COL')]  # Extract network nodes on each loop ID
                    for i, loop_id in enumerate(unique_loop_ids)
                }

                evl_nodes = readEVLtxt(str(ddd_data_dir/'evl'/f'evl_{runID}'))
                node_ids = evl_nodes.nodes[:, get_idx('NODE_ID_COL_IDX_IN_NODE')]
                node_velocity_per_loopID = {loopID: [] for loopID in loop_nodes_per_loopID}
                for loopID, network_nodes in loop_nodes_per_loopID.items():
                    # create lookup hash map out of network nodes
                    node_set = set(network_nodes)
                    # create mask boolean list for loop ID
                    mask = np.isin(node_ids, list(node_set))
                    node_velocity_per_loopID[loopID] = np.mean(evl_nodes.nodes[mask, get_idx('VELOCITY_COLS_IDX_IN_NODE')])*b_SI/convertTimeUnit

                #record_data_temp.extend(node_velocity_per_loopID.values())
                if len(list(node_velocity_per_loopID.values())) == 2:
                    record_data_mean_temp.append(list(node_velocity_per_loopID.values()))
                #print(list(node_velocity_per_loopID.values()))
                #print(np.array(record_data_mean_temp).shape)

                # Skip inhomogeneous arrays
                #if all(len(sublist) == len(record_data_mean_temp[0]) for sublist in record_data_mean_temp):
                #    record_data_mean_temp = np.mean(np.array(record_data_mean_temp), axis=0)


            record_data_mean_temp =np.array(record_data_mean_temp)
            record_data_mean_temp = np.mean(np.array(record_data_mean_temp), axis=0)

            record_data_temp = [seed, int(stress_dir.replace('MPa', ''))]
            record_data_temp.extend(list(record_data_mean_temp))
            record_data_output.append(record_data_temp)

    # Define column headers
    headers = ['Seed', 'Stress', 'mean_velocity_plane1', 'mean_velocity_plane2']
    
    # Create DataFrame
    df = pd.DataFrame(record_data_output, columns=headers)
    
    # Save to CSV
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(output_dir/'velocity_data.csv', index=False)
    print(df)
    
    # Group by Stress and calculate mean for both velocity columns
    df_avg = df.groupby('Stress')[['mean_velocity_plane1', 'mean_velocity_plane2']].mean().reset_index()

    # Add a column for the combined average if needed
    #df_avg['mean_velocity_combined'] = df_avg[['mean_velocity_plane1', 'mean_velocity_plane2']].mean(axis=1)
    print(df_avg)
    df_avg.to_csv(output_dir/'velocity_data_avg.csv', index=False)
    exit()


    # Group by 'Seed' and average all numeric columns
    # Step 2: Average across all seeds
    df_seed_avg = df.groupby('Seed')
    print(df_seed_avg)

    final_avg_plane1 = df['mean_velocity_plane1'].mean()
    final_avg_plane2 = df['mean_velocity_plane2'].mean()
    print(final_avg_plane1)
    print(final_avg_plane2)
    exit()

    # Save averaged data
    os.makedirs(output_dir, exist_ok=True)
    df_avg.to_csv(Path(output_dir)/'velocity_data_averaged.csv', index=False)


if __name__ == "__main__":
    main()
