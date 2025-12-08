# standard lib
from os import path
from posixpath import exists
import re
import sys
import json
import shutil
# 3rd party lib
import numpy as np
from pathlib import Path
from scipy import stats
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


def ci95(data):
    """return mean and 95 % margin for a 1-D array"""
    n = len(data)
    m = np.mean(data)
    se = stats.sem(data)
    h = se * stats.t.ppf((1 + 0.95) / 2, n - 1)
    return m, h

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

    alloy = "AlMg5"
    sid_pair = "sID0-sID1" 
    title_type = f"{alloy}, combined edge"

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
    # for per seed
    total_data_mean = {}
    total_data_std = {}
    # for overall average
    total_data_overall = defaultdict(list)
    for seed, path in paths_per_seed.items():
        paths_stress_sorted = sorted(path, key=lambda p: int(re.search(r"(\d+)MPa.txt", str(p)).group(1)))
        # for per seed
        tmp_dict_mean = {}
        tmp_dict_std = {}
        # for overall average
        #tmp_dict_overall = {}
        for data_path in paths_stress_sorted:
            # extract stress value
            stress = int(re.search(r"(\d+)MPa", str(data_path.name)).group(1))
            # open data
            separation_data = np.loadtxt(data_path)
            tmp_dict_mean[stress] = np.mean(separation_data)
            tmp_dict_std[stress] = np.std(separation_data)
            #tmp_dict_overall[stress] = separation_data
            total_data_overall[stress].append(separation_data)

        total_data_mean[seed] = tmp_dict_mean
        total_data_std[seed] = tmp_dict_std

    # Number of colors you want
    n_colors = 10
    # Use a qualitative colormap like 'tab10' for distinct colors
    colors = cm.tab10(np.linspace(0, 1, n_colors))

    for i, seed in enumerate(total_data_mean.keys()):
        data = total_data_mean[seed]
        data_std = total_data_std[seed]
        stresses = np.array(list(data.keys()))
        sep_means = np.array(list(data.values()))
        sep_std = np.array(list(data_std.values()))
        plt.errorbar(stresses, sep_means, yerr=sep_std, fmt='o', color=colors[i], capsize=5, label=f"seed{seed}")

    plt.axvline(32, linestyle="--", color='k', label="CRSS")
    plt.title(f"Mean Separation Distance {title_type}")
    plt.grid()
    plt.xlabel("Stress [MPa]")
    plt.xlim(-10, 150)
    plt.ylabel("Distance [b]")
    plt.legend()
    #plt.show()
    fig_dir = Path("figures")
    os.makedirs(fig_dir, exist_ok=True)
    plt.savefig(fig_dir/f"{alloy}_{sid_pair}_separation_per_seed.png", transparent=True)
    plt.close("all")

    # Number of colors you want
    #n_colors = 10
    # Use a qualitative colormap like 'tab10' for distinct colors
    #colors = cm.tab10(np.linspace(0, 1, n_colors))
    realization_num = 5
    for stress, sep_data in total_data_overall.items():
        sep_data = np.array(sep_data)
        if len(sep_data)<realization_num:
            continue
        # flatten array
        sep_data = sep_data.reshape(-1)
        print(f"{stress}MPa, {sep_data.shape}")
        mean, error95 = zip(ci95(sep_data))
        #exit()
        plt.errorbar(stress, mean, yerr=error95, fmt='o', color='k', capsize=5, label=f"{stress}MPa")
    plt.axvline(32, linestyle="--", color='r', label="CRSS")
    plt.title(f"Mean Separation Distance {title_type}")
    plt.grid()
    plt.xlabel("Stress [MPa]")
    plt.xlim(-10, 150)
    plt.ylabel("Distance [b]")
    plt.legend()
    plt.savefig(fig_dir/f"{alloy}_{sid_pair}_separation_dist.png", transparent=True)


if __name__ == "__main__":
    main()
