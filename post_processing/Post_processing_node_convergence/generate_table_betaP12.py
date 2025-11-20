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
    'figure.dpi': 200,
    'figure.autolayout': True,  # Prevent label clipping
    'axes.grid': True,
    'grid.alpha': 0.6,
    #'text.usetex': True,   # Enable LaTeX
    'text.usetex': False,   # Enable LaTeX
    'font.size': 20,     # Default font size for text
    'mathtext.fontset': 'stix',  # Use STIX font for math text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    #'font.weight': 'bold',          # Base font weight
    #'text.latex.preamble': r'\usepackage{amsmath} \boldmath',  # Additional packages
    #'text.latex.preamble': r'\usepackage{amsmath}',  # Additional packages
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


#def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
#    #with open(f"{matDir}/{alloy}.txt", "r") as mFile:
#    with open(Path(matDir)/f'{alloy}.txt', "r") as mFile:
#        for line in mFile:
#            # strip down comments from the data
#            line = line.split(";")[0]
#            if line.startswith(f"{var}"):
#                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
#                value = float(line.split("=")[1].strip())
#                break
#    return value


#def plot_data(runIDs: np.ndarray, energy_dict: dict, title: str, fig_fname: str) -> None:
#    n_cols_plot = 2
#    n_rows_plot = int(np.ceil(len(energy_dict.keys()))/n_cols_plot)
#    # Create figure and subplots
#    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(12, 10), dpi=200)  # 2x2 grid (adjust to 4x4 if needed)
#    fig.suptitle(title, fontsize=20, y=1.00)
#
#    # Flatten axes for easy iteration (if using 2x2)
#    axes = axes.ravel()
#
#    # Plot each loop's data
#    for ax, (loop_id, energies) in zip(axes, energy_dict.items()):
#        ax.plot(runIDs, energies, 'o-')
#        ax.set_title(f'Loop ID {loop_id}', fontsize=20)
#        ax.set_xlabel('run ID', fontsize=16)
#        ax.set_ylabel('Energy (J)', fontsize=16)
#        ax.grid(True, linestyle='--', alpha=0.6)
#        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))  # Scientific notation
#
#    # Adjust layout
#    plt.tight_layout()
#    fig.savefig(fig_fname)
#    plt.close()


#def plot_dislocation_loop(loop_pos_per_loopID: dict, save_path: str, runID: str) -> None:
#    n_cols_plot = 2
#    n_rows_plot = int(np.ceil(len(loop_pos_per_loopID.keys()))/n_cols_plot)
#
#    fig, axes = plt.subplots(n_rows_plot, n_cols_plot, figsize=(12, 10), dpi=200)  # 2x2 grid (adjust to 4x4 if needed)
#
#    # Flatten axes for easy iteration
#    axes = axes.ravel()
#
#    # Plot each loop's data
#    for ax, (loop_id, node_pos) in zip(axes, loop_pos_per_loopID.items()):
#        # Vectorized scatter plot with optimized markers
#        xNodePos = node_pos[:, 0]
#        yNodePos = node_pos[:, 1]
#        ax.scatter(xNodePos, yNodePos, s=30, c='k', alpha=0.7, 
#                edgecolors='none', label=f'Loop {loop_id}')
#        ax.set(title=f'loop {loop_id}')
#
#    # Add scientific notation for tick labels
#    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
#
#    os.makedirs(save_path, exist_ok=True)
#    plt.savefig(Path(save_path)/f"{runID}.png", bbox_inches='tight')
#    plt.close()

#def find_glissile_nodes(dataDir, runID):
#    evl = readEVLtxtLoopNode(str(Path(dataDir)/'evl'/f'evl_{runID}'))
#    #evl = readEVLtxtLoopNode(str(Path(dataDir)/'evl'/f'{evl_file}'))
#    loop_data = evl.dislocLoop
#    loop_nodes = evl.loopNodes
#
#    # Get sessile loop IDs (vectorized)
#    sessile_mask = loop_data[:, get_idx('LOOP_TYPE_COL')] == get_idx('SESSILE_TYPE')
#    sessile_loop_ids = loop_data[sessile_mask, get_idx('LOOP_ID_COL')].astype(int)
#
#    # it seems there is a bug in MoDELib2 not marking sessile loops accurately
#    # if sessile loop is not found
#    if not len(sessile_loop_ids):
#        unique_loop_ids, node_count_per_loop = np.unique(loop_nodes[:, get_idx('LOOP_ID_COL')], return_counts=True)
#        all_nodes_per_loop = { int(k): v for k, v in zip(unique_loop_ids, node_count_per_loop) }
#        # hardcoded threshold connecting node number, they are usually 3+3+3 of them, so 10
#        # if you use less than 10 nodes on your glissile dislocations, it breaks the code
#        sessile_loop_ids = [ key for key, count in all_nodes_per_loop.items() if count < 10 ]
#
#    # Get all unique loop IDs (vectorized)
#    all_loop_ids = np.unique(loop_data[:, get_idx('LOOP_ID_COL')]).astype(int)
#
#    # Filter non-sessile loop IDs (vectorized)
#    non_sessile_loop_ids = all_loop_ids[~np.isin(all_loop_ids, sessile_loop_ids)]
#
#    # Get sessile network nodes (vectorized)
#    sessile_network_mask = np.isin(loop_nodes[:, get_idx('LOOP_ID_COL')], sessile_loop_ids)
#    sessile_network_nodes = np.unique(loop_nodes[sessile_network_mask, get_idx('NETWORK_NODE_COL')]).astype(int)
#
#    # Precompute masks for fast filtering
#    non_sessile_node_mask = ~np.isin(loop_nodes[:, get_idx('NETWORK_NODE_COL')], sessile_network_nodes)
#    valid_loop_mask = np.isin(loop_nodes[:, get_idx('LOOP_ID_COL')], non_sessile_loop_ids)
#    combined_mask = non_sessile_node_mask & valid_loop_mask
#
#    # Group nodes by loop ID
#    filtered_nodes = loop_nodes[combined_mask]
#    return filtered_nodes


def main():

    glidePlasticStrain = 1e-10 # abstract value found through trial and error
    alloy = "AlMg5"
    disloc_length = 109 # in burgers vector

    #alloy_mat_file = f"{alloy}1350K"
    #alloy_mat_file = f"{alloy}"

    ddd_dir = Path(sys.argv[1])

    #dataDir = Path(f"./generatedData")
    dataDir = ddd_dir/"generatedData"
    fig_dir = ddd_dir/"figures"
    os.makedirs(fig_dir, exist_ok=True)

    seed_dirs = []
    with os.scandir(dataDir) as dirs:
        for i in dirs:
            if i.is_dir() and "seed" in i.name:
                seed_dirs.append(i.name)
    seed_dirs = sorted(seed_dirs, key=lambda x: int(x.replace('seed', '')))
    
    print(seed_dirs)

    #fig, ax = plt.subplots(figsize=(8, 6))

    total_data = []
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
                betaP_12=getFarray(f,fLabels,'betaP_12').astype(float)
                betaP_13=getFarray(f,fLabels,'betaP_13').astype(float)

                # filter the duplicate data (minimization creates duplicates on runID)
                # Find duplicates and their indices
                unique_values, counts = np.unique(runIDs, return_counts=True)
                duplicates = unique_values[counts > 1]
                # remove the first duplicate data entry
                dup_idx = np.where(runIDs==duplicates)[0][0]
                runIDs = np.delete(runIDs, dup_idx)
                betaP_12 = np.delete(betaP_12, dup_idx)
                #betaP_13 = np.delete(betaP_13, dup_idx)

                # calculate the absolute plastic distortion
                abs_betaP = np.abs(betaP_12 - betaP_12[0])
                #abs_betaP = np.abs(betaP_13 - betaP_13[0])
                #print("Duplicate values:", duplicates)

                dydx = np.gradient(abs_betaP, runIDs)
                dydx_mean = np.mean(dydx[-3:])

                #ax.plot(runIDs, abs_betaP, label=str_dir)
                #ax.plot(runIDs, dydx, label=str_dir)

                # find the stresses above CRSS
                if dydx_mean > glidePlasticStrain:
                    crss_upper_bound = int(str_dir.replace('MPa', ''))
                    break
            non_glide_stresses = tested_strs[tested_strs < crss_upper_bound]
            if len(non_glide_stresses)<1:
                crss = crss_upper_bound
            else:
                #crss_lower_bound = tested_strs[tested_strs < crss_upper_bound].max()
                crss_lower_bound = non_glide_stresses.max()
                search_resolution = crss_upper_bound - crss_lower_bound
                # skip the unfinished realization (larger than 10MPa)
                if search_resolution > 10:
                    #plt.close()
                    continue
                crss = crss_upper_bound - (search_resolution/2)

            print(f"node: {node_dir}, seed {seed_dir}")
            print(f"crss_lower_bound: {crss_lower_bound}, crss_upper_bound: {crss_upper_bound}")
            print(f"crss: {crss}")
            #ax.plot(node_num, crss, label=str_dir)

            node_num = node_dir.replace('node', '')
            seed_num = seed_dir.replace('seed', '')
            total_data.append([seed_num, node_num, crss])

            #ax.legend(frameon=False, fontsize=13)
            #fig.savefig(fig_dir/f"seed{seed_num}_node{node_num}.png")
            #plt.close()

    # save raw data
    header = ['Seed', 'Node', 'CRSS']
    df = pd.DataFrame(total_data, columns=header)
    # Assign data types for each column
    df['Seed'] = df['Seed'].astype(int)
    df['Node'] = df['Node'].astype(int)
    df['CRSS'] = df['CRSS'].astype(float)

    # save as csv file
    df.to_csv(
        ddd_dir/'node_conv_CRSS_table.csv',
        index=False,  # Exclude index column
        float_format='%.3f'
    )

    # plot node number vs CRSS
    unique_seed = df['Seed'].unique()
    fig, ax = plt.subplots(figsize=(8, 6))
    for seed in unique_seed:
        subset = df[df['Seed'] == seed]
        ax.scatter(subset['Node'].astype(float)/disloc_length, subset['CRSS'], label=f"seed {seed}")
    ax.legend()
    ax.set_xlabel('Nodal Density [#/b]')
    ax.set_ylabel('CRSS [MPa]')
    fig.savefig(fig_dir/f"node_vs_CRSS_per_seed.png")
    plt.close()

    # ensemble average node number vs CRSS 
    fig, ax = plt.subplots(figsize=(8, 6))

    # Calculate average CRSS over different seed numbers
    #average_crss = df.groupby(['Seed', 'Node'])['CRSS'].mean().reset_index()
    #average_crss = df.groupby(['Node','Seed'])['CRSS'].mean().reset_index()
    average_crss = df.groupby(['Node'])['CRSS'].mean().reset_index()
    # sort based on the node number and prevent new column generation
    average_crss.sort_values(by='Node').reset_index(drop=True)
    ax.plot(average_crss['Node'].astype(float)/disloc_length, average_crss['CRSS'], '-o')
    ax.set_xlabel('Nodal Density [#/b]')
    ax.set_ylabel('CRSS [MPa]')
    fig.savefig(fig_dir/f"node_vs_CRSS_ensemble.png")


if __name__ == "__main__":
    main()
