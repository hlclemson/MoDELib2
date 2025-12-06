# standard lib
from os import path
import re
import sys
import json
import shutil
# 3rd party lib
import numpy as np
from pathlib import Path
from collections import defaultdict
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 200,
        #"figure.autolayout": True,  # Prevent label clipping
        #"axes.grid": True,
        #"grid.alpha": 0.6,
        "text.usetex": False,
        "font.size": 10,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
    }
)

# ----- MoDELib / utils paths -----
sys.path.append("../../python")
#sys.path.append("../../python/lib")
from visUtils import *
from modlibUtils import *

sys.path.append("../../build/tools/pyMoDELib")
import pyMoDELib

def readValFromMaterialFile(matDir: str, alloy: str, var: str) -> float:
    with open(matDir / alloy, "r") as mFile:
        for line in mFile:
            # strip down comments from the data
            line = line.split(";")[0]
            if line.startswith(f"{var}"):
                # Split the line by '=' and take the second part, then remove leading/trailing whitespace
                value = float(line.split("=")[1].strip())
                break
    return value


def main():
    box_axis_idx = {"x": 0, "y": 1, "z": 2}

    with open("config.json", "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    glide_axis = box_axis_idx[config["glide_axis"]]
    line_tangent_axis = box_axis_idx[config["line_tangent_axis"]]

    data_paths = [p for pat in ("*MPa.txt", "MPa.txt") for p in simulationDir.rglob(pat)]
    # sort by seed number
    data_paths = sorted(
        data_paths, key=lambda p: int(re.search(r"seed(\d+)", str(p)).group(1))
    )

    # sort each paths based on seed
    paths_per_seed = defaultdict(list)
    for i in data_paths:
        seed_num = re.search(r"seed(\d+)", str(i.parent)).group(1)
        paths_per_seed[int(seed_num)].append(i)

    # sort stress value
    #total_data = defaultdict(list)
    total_data_mean = {}
    total_data_std = {}
    # s1: {"stress": 1MPa:, mean":"std":}
    for seed, path in paths_per_seed.items():
        paths_stress_sorted = sorted(path, key=lambda p: int(re.search(r"(\d+)MPa.txt", str(p)).group(1)))
        tmp_dict_mean = {}
        tmp_dict_std = {}
        for data_path in paths_stress_sorted:
            # extract stress value
            stress = int(re.search(r"(\d+)MPa", str(data_path.name)).group(1))
            # open data
            separation_data = np.loadtxt(data_path)
            tmp_dict_mean[stress] = np.mean(separation_data)
            tmp_dict_std[stress] = np.std(separation_data)
        total_data_mean[seed] = tmp_dict_mean
        total_data_std[seed] = tmp_dict_std

    # Number of colors you want
    n_colors = 10
    # Use a qualitative colormap like 'tab10' for distinct colors
    colors = cm.tab10(np.linspace(0, 1, n_colors))

    for i, seed in enumerate(total_data_mean.keys()):
        data = total_data_mean[seed] 
        stresses = np.array(list(data.keys()))
        sep_means = np.array(list(data.values()))
        plt.plot(stresses, sep_means, 'o', color=colors[i], label=f"seed{seed}")
        #exit()
        #for stress in data.keys():
        #    mean = data[stress][0]
        #    std = data[stress][1]
        #    #plt.plot(stress,, label=seed)
        #    #plt.errorbar(stress, mean, yerr=std)
        #    plt.plot(stress, mean, 'o', color=colors[i], label=f"seed{seed}")

    plt.axvline(32, linestyle="--", color='k', label="CRSS")
    plt.title("Mean separation distance (AlMg5, combined edge)")
    plt.grid()
    plt.xlabel("Stress [MPa]")
    plt.xlim(-10, 150)
    plt.ylabel("Distance [b]")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
