#from logging import logProcesses
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
matplotlib.use('Agg')
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
    'grid.alpha': 0.3,
    'text.usetex': True,   # Enable LaTeX
    'font.size': 30,     # Default font size for text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    'text.latex.preamble': r'\usepackage{amsmath}'  # Additional packages
})

def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
    with open(f"{matDir}/{alloy}.txt", "r") as mFile:
        for line in mFile:
            # strip down comments from the data
            line = line.split(";")[0]
            if line.startswith(f"{var}"):
                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
                value = float(line.split("=")[1].strip())
                break
    return value


def parseNumNodes(dir: str, slipIDnum: list) -> str:
    # nodeEntryString = dir.split('N')[1].split('DP')[0]
    nodeEntryPattern = r"N(\d+)DP"
    match = re.search(nodeEntryPattern, dir)
    if match:
        nodeStrTemp = match.group(1)
        nodeStr = [x for x in nodeStrTemp]
        # partial dislocation has more than 1 slipsys ID
        if len(slipIDnum) > 1:
            nodeNum = nodeStr[0 : int(len(nodeStr) / 2)]
            nodeNum = "".join(nodeNum)
        # full dislocation has only 1 slipsys ID
        else:
            nodeNum = "".join(nodeStr)
    return nodeNum


def extract_evl_number(fname: str) -> int:
    """Extract numeric portion from 'evl_XXXXX.txt' filenames"""
    match = re.search(r'evl_(\d+)', fname)
    return int(match.group(1)) if match else 0


def plot_dislocation_loop(evl_runID: str, loop_id: int, node_positions: np.ndarray, save_path=None) -> None:
    """Vectorized plotting function for a single loop."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Vectorized scatter plot with optimized markers
    xNodePos = node_positions[:, 0]
    yNodePos = node_positions[:, 1]
    ax.scatter(xNodePos, yNodePos, s=30, c='k', alpha=0.7, 
               edgecolors='none', label=f'Loop {loop_id}')

    ax.set(title=f'Dislocation Loop Nodes (Loop {loop_id})',
           xlabel='x',
           ylabel='y')

    # Add scientific notation for tick labels
    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(f"{save_path}/{evl_runID}_loop_{loop_id}.png", bbox_inches='tight', transparent=True)
    plt.close()

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


def plot_power_spectrum_one(k_exponent: float, intercept: float, k: np.ndarray, sxx: np.ndarray, pname: str, save_path=None) -> None:
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
    ax.plot(k_filtered, 10**trendline, "r--", alpha=0.6)  # Convert the trendline back to linear scale

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
           xlabel='Wave Vector (k)',
           ylabel='Power Spectral Density')

    # Add scientific notation for tick labels
    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    fig.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        #plt.savefig(f"{save_path}/average_psd.png", bbox_inches='tight', transparent=True)
        plt.savefig(f"{save_path}/{pname}", bbox_inches='tight', transparent=True)
    plt.close()

def plot_power_spectrum(loop_id: int, k: np.ndarray, sxx: np.ndarray, save_path=None) -> None:
    """Vectorized plotting function for a single loop."""
    fig, ax = plt.subplots(figsize=(10, 8))

    # Filter out zeros/negative values before log transform
    mask = (k > 0) & (sxx > 0)
    k_filtered = k[mask]
    sxx_filtered = sxx[mask]

    # Log transform the data
    log_k = np.log10(k_filtered)
    log_sxx = np.log10(sxx_filtered)

    # Perform linear regression on the log-transformed data
    slope, intercept, r_value, p_value, std_err = linregress(log_k, log_sxx)

    # Generate trendline points (in log space)
    trendline = slope * log_k + intercept

    # Plot the trendline
    ax.plot(k_filtered, 10**trendline, "r--", alpha=0.6)  # Convert the trendline back to linear scale

    # Calculate the equation of the trendline
    #equation = f'P(k) = {slope:.6f}k + {intercept:.6f}'
    equation = f'$P(k) \\propto k^{{{slope:.2f}}}$'
    ax.text(0.55, 0.95, equation, transform=ax.transAxes, 
        fontsize=26, va='top')

    # Add the equation as an annotation on the plot
    #ax.annotate(equation, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12, ha='left', va='top')
    ax.scatter(k, sxx, s=13, c='k', alpha=0.9, 
               edgecolors='none', label=f'Loop {loop_id}')
   
    # Log-scale settings
    ax.set(xscale='log', yscale='log',
           title=f'Power Spectrum (Loop {loop_id})',
           xlabel='Wave Vector (k)',
           ylabel='Power Spectral Density')

    # Add scientific notation for tick labels
    #ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    fig.tight_layout()

    if save_path:
        os.makedirs(save_path, exist_ok=True)
        plt.savefig(f"{save_path}/loop_{loop_id}.png", bbox_inches='tight', transparent=True)
    plt.close()

# Combine all temporary databases into one final database
def merge_databases(dName: str, process_ids: list):
    final_db = f"{dName}.db"

    with sqlite3.connect(final_db) as main_conn:
        main_cursor = main_conn.cursor()
        # nuke exisitng data
        #main_cursor.execute("DROP TABLE IF EXISTS CRSSdata")
        # Create final table
        #main_cursor.execute("""
        #    CREATE TABLE IF NOT EXISTS PSD (
        #        dislocationType TEXT NOT NULL,
        #        noiseType TEXT NOT NULL,
        #        noiseSeed INTEGER,
        #        dislocationCharacter TEXT NOT NULL,
        #        alloy TEXT NOT NULL,
        #        dislocationLength INTEGER,
        #        stress INTEGER,
        #        kExponent REAL
        #    )
        #""")
        main_cursor.execute("""
            CREATE TABLE IF NOT EXISTS PSD (
                dislocationType TEXT NOT NULL,
                noiseType TEXT NOT NULL,
                dislocationCharacter TEXT NOT NULL,
                alloy TEXT NOT NULL,
                dislocationLength INTEGER,
                stress INTEGER,
                kExponent REAL
            )
        """)

        # Attach and merge each temporary database
        for pid in process_ids:
            temp_db = f"{dName}_{pid}.db"
            # Skip if file doesn't exist
            if not os.path.exists(temp_db):
                print(f"Skipping missing file: {temp_db}")
                continue
            try:
                # Attach the database
                main_cursor.execute(f"ATTACH DATABASE '{temp_db}' AS temp_db")

                # Check if CRSSdata table exists
                main_cursor.execute("""
                    SELECT count(*) 
                    FROM temp_db.sqlite_master 
                    WHERE type='table' AND name='PSD'
                """)
                table_exists = main_cursor.fetchone()[0] > 0

                if table_exists:
                    # Count rows before merge
                    main_cursor.execute("SELECT COUNT(*) FROM temp_db.PSD")
                    row_count = main_cursor.fetchone()[0]

                    # Perform the merge
                    main_cursor.execute("""
                        INSERT INTO main.PSD
                        SELECT * FROM temp_db.PSD
                    """)
                    main_conn.commit()
                    print(f"Merged {row_count} rows from {temp_db}")
                else:
                    print(f"Skipping {temp_db} - no PSD table found")

            except sqlite3.Error as e:
                print(f"Error processing {temp_db}: {str(e)}")
            finally:
                # Always detach and clean up
                main_cursor.execute("DETACH DATABASE temp_db")
                if os.path.exists(temp_db):
                    os.remove(temp_db)

def save_SQL_db_crss(fname: str, data_set: list, pid: int) -> None:
    #with sqlite3.connect(f"k_power_data_crss_{pid}.db") as conn:
    with sqlite3.connect(f"{fname}_{pid}.db") as conn:
        cursor = conn.cursor()
        # Create table if not exists (each process needs its own schema)
        #cursor.execute("""
        #    CREATE TABLE IF NOT EXISTS PSD (
        #        dislocationType TEXT NOT NULL,
        #        noiseType TEXT NOT NULL,
        #        noiseSeed INTEGER,
        #        dislocationCharacter TEXT NOT NULL,
        #        alloy TEXT NOT NULL,
        #        dislocationLength INTEGER,
        #        stress INTEGER,
        #        kExponent REAL
        #    )
        #""")
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS PSD (
                dislocationType TEXT NOT NULL,
                noiseType TEXT NOT NULL,
                dislocationCharacter TEXT NOT NULL,
                alloy TEXT NOT NULL,
                dislocationLength INTEGER,
                stress INTEGER,
                kExponent REAL
            )
        """)
        # Insert data
        cursor.executemany(
            "INSERT INTO PSD VALUES (?, ?, ?, ?, ?, ?, ?)",
            data_set
        )
        conn.commit()


def export_from_sql_to_csv(fname: str):
    # Connect to the SQLite database
    conn = sqlite3.connect(f'{fname}.db')
    cursor = conn.cursor()
    
    # Query the table you want to export
    cursor.execute("SELECT * FROM PSD")
    rows = cursor.fetchall()
    
    # Write to CSV
    with open(f'{fname}.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # Write header (column names)
        csvwriter.writerow([i[0] for i in cursor.description])
        # Write data rows
        csvwriter.writerows(rows)
    
    conn.close()


def save_SQL_db_below_crss(data_set: list, pid: int):
    with sqlite3.connect(f"k_power_data_below_crss_{pid}.db") as conn:
        cursor = conn.cursor()
        # Create table if not exists (each process needs its own schema)
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS PSD (
                dislocationType TEXT NOT NULL,
                noiseType TEXT NOT NULL,
                noiseSeed INTEGER,
                dislocationCharacter TEXT NOT NULL,
                alloy TEXT NOT NULL,
                dislocationLength INTEGER,
                stress INTEGER,
                kExponent REAL
            )
        """)
        # Insert data
        cursor.executemany(
            "INSERT INTO PSD VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            data_set
        )
        conn.commit()

def calcShape(
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
    nodeNumbers = []
    with os.scandir(genDataDir) as it:
        for entry in it:
            if entry.is_dir():
                node = parseNumNodes(entry.name, slipIDnum)
                if node not in nodeNumbers:
                    nodeNumbers.append(node)
                else:
                    continue

    binnedPerNodes = defaultdict(list)
    for nod in nodeNumbers:
        temp = []
        for dir in dataDirs:
            # check if the directory has the matching node
            dirNodeNumPattern = r"N(\d+)DP"
            match = re.search(dirNodeNumPattern, dir)
            if match:
                dirNodeNum = match.group(1)
            # partial or full
            if len(slipIDnum) > 1:
                if f"{nod}{nod}" == dirNodeNum:
                    binnedPerNodes[nod].append(dir)
            else:
                if f"{nod}" == dirNodeNum:
                    binnedPerNodes[nod].append(dir)
                # temp.append(dir)
        # binnedPerNodes[nod].append(temp)

    binnedPerNodeAndSeedDirs = {}
    for nod in nodeNumbers:
        tempDirsPerSeed = {}
        dirBin = binnedPerNodes[nod]
        for seed in seeds:
            temp = []
            for dir in dirBin:
                seed = int(seed)
                if int(dir.split("S")[-1]) == seed:
                    # tempDirs[seed].append(dir)
                    temp.append(dir)
            tempDirsPerSeed[seed] = temp
        binnedPerNodeAndSeedDirs[nod] = tempDirsPerSeed

    for nod in nodeNumbers:
        nodeDirs = binnedPerNodeAndSeedDirs[nod]
        wave_vectors_avg_over_seed = []
        power_spectra_avg_over_seed = []
        for seed in seeds:
            dirsPerSeed = nodeDirs[seed]
            # skip the empty dictionary
            if not dirsPerSeed:
                continue

            # Sort the directories based on the stress values
            sortedDirsOnSTR = sorted(
                dirsPerSeed, key=lambda x: int(x.split("T")[0].split("r")[1])
            )

            #prevStressVal = None
            #dydx_stress = {}
            for dIdx, dir in enumerate(sortedDirsOnSTR):
                alloyMatch = re.search(r"AlMg(\d+)", dir)
                alloy = alloyMatch.group(0)
                matDir = f"{mainDataDir}/{dataTypePath}/MoDELib2/tutorials/DislocationDynamics/MaterialsLibrary/"
                # extract material info
                b_SI = readValFromMaterialFile(matDir, alloy, "b_SI")
                mu0_SI = readValFromMaterialFile(matDir, alloy, "mu0_SI")
                rho_SI = readValFromMaterialFile(matDir, alloy, "rho_SI")
                cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
                convertTimeUnit = b_SI / cs  # [sec]
                convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

                copyFlag = False
                counter = 0
                fTensor = np.zeros((3, 3))

                try:
                    with open(
                        f"{genDataDir}/{dir}/inputFiles/polycrystal.txt", "r"
                    ) as poly:
                        for line in poly:
                            line = line.strip("\n")
                            if "F=" in line:
                                copyFlag = True
                            if copyFlag:
                                if counter == 0:
                                    tensorLine = np.array(
                                        line.split("=")[1].split(" "), dtype=np.float64
                                    )
                                elif counter == 1:
                                    tensorLine = np.array(
                                        line.strip().split(" "), dtype=np.float64
                                    )
                                else:
                                    tensorLine = np.array(
                                        line.split(";")[0].strip().split(" "),
                                        dtype=np.float64,
                                    )
                                fTensor[counter] += tensorLine
                                counter += 1
                            if counter == 3:
                                copyFlag = False
                                break
                except FileNotFoundError:
                    print(f"File not found: {genDataDir}/{dir}/inputFiles/polycrystal.txt")
                    continue
                except IOError as e:
                    print(f"Error reading file: {e}")
                    continue

                V = np.linalg.det(fTensor)
                fDiag = np.diag(fTensor)
                A = fDiag[0] * fDiag[1]

                # read labels
                with open(f"{genDataDir}/{dir}/F/F_labels.txt", "r") as label:
                    fLabels = label.read()
                # remove empty element and store it as a list
                fLabels = np.array([x for x in fLabels.split("\n") if x], dtype=str)

                # last 33 element indexes
                lastCols = np.arange(-33, 0)
                # open F file
                fData = defaultdict(list)
                with open(f"{genDataDir}/{dir}/F/F_0.txt", "r") as f:
                    for line in f:
                        line = [float(x) for x in line.split(" ") if x and x != "\n"]
                        # parse the first 14 elements
                        for i in range(len(fLabels[:14])):
                            fData[str(fLabels[i])].append(line[i])
                        # parse the last 33 elements
                        for j in lastCols:
                            fData[str(fLabels[j])].append(line[j])

                xAxisData = np.array(fData[xAxisDataLabel])
                yAxisData = np.array(fData[yAxisDataLabel])

                minRunidIdx = np.squeeze(
                    np.where(xAxisData == np.float64(minimizationSteps))
                )
                stressVal = int(re.search(r"Str(\d+)", dir).group(1))

                # devide the plastic distortion by two (dipole)
                # subtract the plastic distortion caused by minimization step
                yAxisData = yAxisData - yAxisData[minRunidIdx]

                # absolute value
                xAxisData = np.abs(xAxisData)
                yAxisData = np.abs(yAxisData)

                print(f"processing {dir}..")

                # if there is a data error, skip the data
                if len(xAxisData) != len(yAxisData):
                    continue

                # Get all evl file names
                local_evl_dir = Path(genDataDir) / dir / 'evl'
                # drop the .txt
                all_evl_fnames = [f.name.split('.')[0] for f in local_evl_dir.glob('evl_*.txt') if f.is_file()]
                # Get last 50% of sorted files
                sorted_evl_runIDs = sorted(all_evl_fnames, key=extract_evl_number)[len(all_evl_fnames) // 2:]

                wave_vectors_avg = defaultdict(list)
                power_spectra_avg = defaultdict(list)
                for evl_runID in sorted_evl_runIDs:
                    evl_file = str(Path(local_evl_dir) / evl_runID)
                    evl=readEVLtxtLoopNode(evl_file)
                    loop_data = evl.dislocLoop
                    loop_nodes = evl.loopNodes

                    # Constants for column indices
                    LOOP_ID_COL = 0
                    LOOP_TYPE_COL = 11
                    NETWORK_NODE_COL = 5
                    SESSILE_TYPE = 1

                    # Get sessile loop IDs (vectorized)
                    sessile_mask = loop_data[:, LOOP_TYPE_COL] == SESSILE_TYPE
                    sessile_loop_ids = loop_data[sessile_mask, LOOP_ID_COL].astype(int)
                    # if sessile loop is not found
                    if not len(sessile_loop_ids):
                        unique_loop_ids, node_count_per_loop = np.unique(loop_nodes[:, LOOP_ID_COL], return_counts=True)
                        all_nodes_per_loop = { int(k): v for k, v in zip(unique_loop_ids, node_count_per_loop) }
                        # hardcoded threshold connecting node number, they are usually 3+3+3 of them, so 10
                        sessile_loop_ids = [ key for key, count in all_nodes_per_loop.items() if count < 10 ]

                    # Get all unique loop IDs (vectorized)
                    all_loop_ids = np.unique(loop_data[:, LOOP_ID_COL]).astype(int)

                    # Filter non-sessile loop IDs (vectorized)
                    non_sessile_loop_ids = all_loop_ids[~np.isin(all_loop_ids, sessile_loop_ids)]

                    # Get sessile network nodes (vectorized)
                    sessile_network_mask = np.isin(loop_nodes[:, LOOP_ID_COL], sessile_loop_ids)
                    sessile_network_nodes = np.unique(loop_nodes[sessile_network_mask, NETWORK_NODE_COL]).astype(int)

                    # Precompute masks for fast filtering
                    non_sessile_node_mask = ~np.isin(loop_nodes[:, NETWORK_NODE_COL], sessile_network_nodes)
                    valid_loop_mask = np.isin(loop_nodes[:, LOOP_ID_COL], non_sessile_loop_ids)
                    combined_mask = non_sessile_node_mask & valid_loop_mask

                    # Group nodes by loop ID (vectorized)
                    filtered_nodes = loop_nodes[combined_mask]
                    unique_loop_ids, indices = np.unique(filtered_nodes[:, LOOP_ID_COL], return_inverse=True)
                    loopNodesPerLoopID = {
                        int(loop_id): filtered_nodes[indices == i, 2:5]  # Extract positions (columns 2-4)
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
                    lineTangentIdx = dislocGlideLineTangentIdx[disloc_character]["lineTangentIdx"]
                    glideDirIdx = dislocGlideLineTangentIdx[disloc_character]["glideDirectionIdx"]

                    # Stack all node positions and centers for batch processing
                    all_nodes = np.concatenate(list(loopNodesPerLoopID.values()))  # Shape: (total_nodes, 3)
                    all_centers = np.repeat(centers_of_mass, [len(pos) for pos in node_positions], axis=0)
                    # Precompute all deviations in the glide direction
                    deviations = np.abs(all_nodes[:, glideDirIdx] - all_centers[:, glideDirIdx])

                    # Split deviations back per loop ID using segment lengths
                    split_indices = np.cumsum([len(pos) for pos in node_positions])[:-1]
                    deviation_per_loop = np.split(deviations, split_indices)

                    # plot raw node position at a given evl file (for debugging)
                    #for loop_id in loop_ids:
                    #    if loop_id in loopNodesPerLoopID:  # Check existence first
                    #        plot_dislocation_loop(evl_runID, loop_id, loopNodesPerLoopID[loop_id],
                    #                        save_path=Path('./DislocationLoopShape')/dataTypePath/dir)

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

                # average over all the loops
                #wave_vectors_avg_loop_avg = [ np.array(x) for x in wave_vectors_avg.values() ]
                wave_vectors_avg_loop_avg = np.mean(np.array([v for v in wave_vectors_avg.values()]), axis=0)
                power_spectra_avg_loop_avg = np.mean(np.array([v for v in power_spectra_avg.values()]), axis=0)

                # fit power law curve
                k_exponent, intercept = fit_power_law(wave_vectors_avg_loop_avg, power_spectra_avg_loop_avg)

                # parse current stress value
                match = re.search(r'Str(\d+)T', dir)
                stress = int(match.group(1))

                # create SQL data
                #k_slope_data = [
                #    (
                #        disloc_type,
                #        noise_type,
                #        int(seed),
                #        disloc_character,
                #        alloy,
                #        int(disloc_length),
                #        int(stress),
                #        k_exponent
                #    )
                #]

                #save_SQL_db_crss(sql_fname, k_slope_data, pid)

                plot_power_spectrum_one(k_exponent, intercept, wave_vectors_avg_loop_avg, power_spectra_avg_loop_avg,
                                f"PSD_{appliedStr}.png", save_path=Path('./Powerspectrum')/dataTypePath/dir)
            wave_vectors_avg_over_seed.append(wave_vectors_avg_loop_avg)
            power_spectra_avg_over_seed.append(power_spectra_avg_loop_avg)

        # average over seeds
        wave_vectors_avg_over_seed = np.mean(np.array(wave_vectors_avg_over_seed), axis=0)
        power_spectra_avg_over_seed= np.mean(np.array(power_spectra_avg_over_seed), axis=0)
        np.save(Path('./Powerspectrum')/dataTypePath/f'wave_vectors_avg_{appliedStr}.npy', wave_vectors_avg_over_seed)
        np.save(Path('./Powerspectrum')/dataTypePath/f'power_spectra_avg_{appliedStr}.npy', power_spectra_avg_over_seed)

        # fit power law curve
        k_exponent, intercept = fit_power_law(wave_vectors_avg_over_seed, power_spectra_avg_over_seed)
        plot_power_spectrum_one(k_exponent, intercept, wave_vectors_avg_over_seed, power_spectra_avg_over_seed,
                        f"PSD_{appliedStr}.png", save_path=Path('./Powerspectrum')/dataTypePath)

        k_slope_data = [
            (
                disloc_type,
                noise_type,
                disloc_character,
                alloy,
                int(disloc_length),
                int(stress),
                k_exponent
            )
        ]
        save_SQL_db_crss(sql_fname, k_slope_data, pid)

def main() -> None:

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
        mainDataFolder = config["dataFolder"]

    #figDataDir = f"{mainDataFolder}/Figure/"
    figDataDir = f"./Figure/"
    dd2OvitoBin = "./MoDELib2/tools/DD2OvitoVtk/build/DD2OvitoVtk"

    # Using itertools.product to generate all combinations
    dataPathComponents = itertools.product(
        dislocTypes.items(), dislocChars, alloys, dislocLens
    )

    appliedStr = "200MPa"
    sql_fname = f"k_power_data_{appliedStr}"

    # Prepare argument tuples for starmap
    args_list = [
        (*args, pid)  # Append PID to each argument tuple
        for pid, args in enumerate([
            (
                appliedStr,
                mainDataFolder,
                dtype,
                noise,
                dchar,
                alloy,
                dlen,
                dislocGlideLineTangentIdx,
                xAxisDataLabel,
                yAxisDataLabel,
                minimizationSteps,
                sql_fname,
            )
            for (dtype, noise_types), dchar, alloy, dlen in dataPathComponents
            for noise in noise_types["noiseTypes"]
        ])
    ]

    # serial run (useful for testing)
    #for args in args_list:
    #    calcShape(*args)
    #exit()

    # multiprocessing
    with Pool() as pool:
        pool.starmap(calcShape, args_list)

    # Get unique process IDs based on the argument lists
    pids = [arg[-1] for arg in args_list]

    # Merge all databases at the end
    merge_databases(sql_fname, pids)
    export_from_sql_to_csv(sql_fname)


if __name__ == "__main__":
    main()
