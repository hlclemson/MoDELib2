import re
import os
import sys
import string
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


# Configure global plot settings (applies to all figures)
rcParams.update({
    'figure.dpi': 200,
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


def main():

    alloy = "NiCoCr"
    alloy_mat_file = f"{alloy}1350K"
    #alloy_mat_file = f"{alloy}"
    #dataDir = f"../xin_SF_noise_{alloy}_Initial_separation_study/generatedData/seed1/"
    #dataDir = Path(f"../xin_SF_noise_{alloy}_Initial_separation_study/generatedData/seed1/")
    dataDir = Path(f"../xin_SF_noise_{alloy}_Initial_separation_study/generatedData/")
    output_dir = dataDir/"separation_distance"

    ensem_seed_dirs = []
    with os.scandir(dataDir) as seed_dirs:
        for dir in seed_dirs:
            if dir.is_dir() and "seed" in dir.name:
                ensem_seed_dirs.append(dir.name)

    # sort the directory list based on the seed number
    ensem_seed_dirs = sorted(ensem_seed_dirs, key=lambda x: int(x.replace('seed', '')))

    data_table_recorded = []
    for seed_dir in ensem_seed_dirs:

        init_sep_dirs = []
        with os.scandir(dataDir/seed_dir) as sep_dirs:
            for dir in sep_dirs:
                if dir.is_dir() and "b" in dir.name:
                    init_sep_dirs.append(dir.name)

        # Using regex to find all digits between 'sep' and 'b' and sort based on burgers vector
        init_sep_dirs = sorted(init_sep_dirs, key=lambda x: int(re.search(r'sep(\d+)b', x).group(1)))
        print(init_sep_dirs)

        for init_sep_dir in init_sep_dirs:
            # temporary data dir
            temp_data_dir = dataDir/seed_dir/init_sep_dir

            # load f file and runIDs
            f,fLabels=readFfile(str(temp_data_dir/'F'))
            runIDs=getFarray(f,fLabels,'runID').astype(int)

            # get the final runID
            runID = runIDs[-1]
            #for runID in runIDs:

            ## Group nodes by loop ID
            filtered_nodes = find_glissile_nodes(temp_data_dir, runID)

            unique_loop_ids, indices = np.unique(filtered_nodes[:, get_idx('LOOP_ID_COL')], return_inverse=True)
            loopNodesPerLoopID = {
                int(loop_id): filtered_nodes[indices == i, get_idx('POS_COLS')]  # Extract network nodes on each loop ID
                for i, loop_id in enumerate(unique_loop_ids)
            }

            # find the partial pair (on the same glide plane)
            plane_height_per_loopid = {}
            for loopID, node_pos in loopNodesPerLoopID.items():
                plane_height_per_loopid[loopID] = np.mean(node_pos[:, get_idx('PLANE_NORM_DIR_AXIS')], dtype=int)

            loopid_pair_on_plane = defaultdict(list)
            for loopID, height in plane_height_per_loopid.items():
                loopid_pair_on_plane[height].append(loopID)

            # find the node pair (on the same line tangent direction coordinate)
            partial_dists_all_planes = []
            for height, pair_loop_ids in loopid_pair_on_plane.items():
                #partial_pair_pos = []
                loop1, loop2 = pair_loop_ids
                node_pos_loop1 = loopNodesPerLoopID[loop1]
                node_pos_loop2 = loopNodesPerLoopID[loop2]

                partial_dists = np.abs(node_pos_loop1[:, get_idx('GLIDE_DIR_AXIS')] - node_pos_loop2[:, get_idx('GLIDE_DIR_AXIS')])
                partial_dists_all_planes.append(partial_dists)
            partial_dists_all_planes = np.array(partial_dists_all_planes)

            matDir = str(temp_data_dir/'inputFiles') # extract material info
            b_SI = readValFromMaterialFile(matDir, alloy_mat_file, "b_SI")
            to_Ang = 1e-10

            # convert distance unit from DDD to Angstrom
            partial_dists_all_planes *= b_SI/to_Ang
            #plot_distribution(partial_dists_all_planes, output_dir, init_sep_dir)

            # combine both glide plane data
            partial_dists_all_planes = partial_dists_all_planes.flatten()

            # calculate statistics
            mean = np.mean(partial_dists_all_planes) 
            std = np.std(partial_dists_all_planes)

            # initial separation distance
            init_sep = float(re.search(r'sep(\d+)b', init_sep_dir).group(1))*b_SI/to_Ang

            # seed number
            seed_num = int(seed_dir.replace('seed', ''))

            data_record = [alloy, seed_num, init_sep, mean, std]
            data_table_recorded.append(data_record)

    # dataframe
    # alloy separation_dist glideplane1 mean gp1 std  gp2 mean gp2 std  
    os.makedirs(output_dir, exist_ok=True)

    # Create DataFrame with column names
    header = ['Alloy', 'Noise_seed', 'Initial_sep [AA]', 'Sep_dist_mean [AA]', 'Sep_dist_std [AA]']
    df = pd.DataFrame(data_table_recorded, columns=header)

    # save as csv file
    df.to_csv(
        Path(output_dir)/'initial_separation_data.csv',
        index=False,  # Exclude index column
        float_format='%.3f'  # Scientific notation with 4 decimal places
    )


    print(df)

    # Plot with error bars
    plt.errorbar(
        x=df['Initial_sep [AA]'],
        y=df['Sep_dist_mean [AA]'],
        yerr=df['Sep_dist_std [AA]'],
        fmt='o--',  # circle markers with dashed line
        capsize=5,  # error bar cap width
        elinewidth=2,  # error bar line width
        markersize=4,
        label='NiCoCr'
    )

    # Formatting
    #plt.xscale('log')  # Log scale for x-axis (common for strain rates)
    #plt.xlabel('Initial Separation [Å]', fontsize=12)
    #plt.ylabel('Mean Separation Distance [Å]', fontsize=12)
    ##plt.title('', fontsize=14)
    #plt.grid(True, which="both", ls="--", alpha=0.3)
    #plt.legend(fontsize=12)
    plt.xlabel('Initial Separation [Å]')
    plt.ylabel('Mean Separation Distance [Å]')
    #plt.title('', fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.3)
    plt.legend()

    # Show plot
    plt.tight_layout()
    plt.savefig(output_dir/'init_vs_mean.png')


if __name__ == "__main__":
    main()
