# standard lib
import os
from os import error, path
from posixpath import exists
import re
import sys
import json
import shutil
# 3rd party lib
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from collections import defaultdict

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

    exit1_01_lead0 = {
        "dirs":[
            "exitID1_partial_orient_615AA_crss_0-1_almg5_total",
            "exitID1_partial_orient_615AA_crss_0-1_almg10_total",
            "exitID1_partial_orient_615AA_crss_0-1_almg15_total"
        ],
        "crss_dirs":[
            "exitID1_partial_orient_615AA_crss_0-1_almg5_total",
            "exitID1_partial_orient_615AA_crss_0-1_almg10_total",
            "exitID1_partial_orient_615AA_crss_0-1_almg15_total"
        ],
        "leading_partial": 0
    }
    exit1_23_lead3 = {
        "dirs":[
            "exitID1_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg15_total",
        ], 
        "crss_dirs":[
            "exitID1_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg15_total"
        ],
        "leading_partial": 3
    }
    exit1_23_lead2 = {
        "dirs":[
            "exitID1_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID1_partial_orient_615AA_crss_2-3_almg15_total",
        ], 
        "crss_dirs":[
            "exitID1_partial_orient_615AA_crss_2-3_lead2_almg5_total",
            "exitID1_partial_orient_615AA_crss_2-3_lead2_almg10_total",
            "exitID1_partial_orient_615AA_crss_2-3_lead2_almg15_total"
        ],
        "leading_partial": 2
    }
    exit2_01_lead0 = {
        "dirs":[
            "exitID2_partial_orient_615AA_crss_0-1_almg5_total",
            "exitID2_partial_orient_615AA_crss_0-1_almg10_total",
            "exitID2_partial_orient_615AA_crss_0-1_almg15_total",
        ],
        "crss_dirs":[
            "exitID2_partial_orient_615AA_crss_0-1_almg5_total",
            "exitID2_partial_orient_615AA_crss_0-1_almg15_total",
            "exitID2_partial_orient_615AA_crss_0-1_almg10_total"
        ],
        "leading_partial": 0
    }
    exit2_23_lead3 = {
        "dirs":[
            "exitID2_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg15_total",
        ], 
        "crss_dirs":[
            "exitID2_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg15_total"
        ],
        "leading_partial": 3
    }
    exit2_23_lead2 = {
        "dirs":[
            "exitID2_partial_orient_615AA_crss_2-3_almg5_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg10_total",
            "exitID2_partial_orient_615AA_crss_2-3_almg15_total",
        ],
        "crss_dirs":[
            "exitID2_partial_orient_615AA_crss_2-3_lead2_almg5_total",
            "exitID2_partial_orient_615AA_crss_2-3_lead2_almg10_total",
            "exitID2_partial_orient_615AA_crss_2-3_lead2_almg15_total"
        ],
        "leading_partial": 2
    }
    data_dir = Path("../../production_runs")
    data_dir_other_dir = Path("../../production_runs_other_dir/")

    crss_data_dir = Path("../create_crss_table/output")
    fig_dir = Path("figures")

    data_frames = []
    data_columns = ["exitID", "sIDs", "leading_partial_id", "Mg", "stress", "separation_mean", "separation_ci95", "data_size"]

    for exit_pairs in [exit1_01_lead0, exit1_23_lead3, exit1_23_lead2, exit2_01_lead0, exit2_23_lead3, exit2_23_lead2]:
        print(exit_pairs)
        exit_pairs_dirs = exit_pairs["dirs"]
        leading_p = exit_pairs["leading_partial"]
        crss_paths = exit_pairs["crss_dirs"]
        #total_data_overall_all_almg = defaultdict(list)
        total_data_overall_all_almg = {}
        for gpath, crss_path in zip(exit_pairs_dirs, crss_paths):
            data_dir_local = data_dir/gpath/"generatedData"
            crss_data_path_local = crss_data_dir/crss_path/"all_crss_data.csv"
            crss_data_local = pd.read_csv(crss_data_path_local) 
            avg_crss = crss_data_local["crss"].mean()
            data_paths = [p for pattern in ("*MPa.txt", "MPa.txt") for p in data_dir_local.rglob(pattern)]
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
            alloy = f"AlMg{re.search(r"almg(\d+)", str(gpath)).group(1)}" 
            exit_id = f"exitID{re.search(r"exitID(\d+)", str(gpath)).group(1)}"
            sid_pair = re.search(r"crss_(\d)-(\d)", str(gpath))
            sid_pair = f"sID{sid_pair.group(1)}-{sid_pair.group(2)}"
            leading_par = f"lead{leading_p}"
            for seed, path in paths_per_seed.items():
                paths_stress_sorted = sorted(path, key=lambda p: int(re.search(r"(\d+)MPa.txt", str(p)).group(1)))
                # for overall average
                for data_path in paths_stress_sorted:
                    # extract stress value
                    stress = int(re.search(r"(\d+)MPa", str(data_path.name)).group(1))
                    # open data
                    separation_data = np.loadtxt(data_path)
                    if not separation_data.any():
                        continue
                    total_data_overall[stress].append(separation_data)

            for stress, sep_data in total_data_overall.items():
                sep_data = np.asarray(sep_data)
                sep_data = sep_data.reshape(-1)
                data_size = len(sep_data)
                #mean, error95 = zip(ci95(sep_data))
                mean, error95 = ci95(sep_data)
                mg = re.search(r"AlMg(\d+)", str(alloy)).group(1)
                sid_pair_tmp = re.search(r"sID(\d)-(\d)", str(sid_pair))
                new_entry = [(
                    re.search(r"exitID(\d+)", str(exit_id)).group(1),
                    f"{sid_pair_tmp.group(1)}-{sid_pair_tmp.group(2)}",
                    leading_p,
                    mg,
                    stress,
                    mean,
                    error95,
                    data_size
                )]

                data_frames.append(pd.DataFrame(new_entry, columns=data_columns))

    data_frames = pd.concat(data_frames, ignore_index=True)
    data_out_dir = Path("generated_data")
    os.makedirs(data_out_dir, exist_ok=True)
    data_frames.to_csv(data_out_dir/"separation_distance_total_data.csv", index=False)


if __name__ == "__main__":
    main()
