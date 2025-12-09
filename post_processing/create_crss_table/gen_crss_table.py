# standard lib
import os
import re
import sys
import json
import shutil
import tarfile
import tempfile

# 3rd party lib
import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict


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
    try:
        config_fname = sys.argv[1]
    except:
        exit("usage) gen_crss_table.py config.json")

    with open(config_fname, "r") as f:
        config = json.load(f)
    simulationDir = Path(config["data_path"])
    output_dir = Path(config["output_path"])

    data_paths = [p for pat in ("*MPa.tar.gz", "MPa") for p in simulationDir.rglob(pat)]
    # sort by seed number
    data_paths = sorted(
        data_paths, key=lambda p: int(re.search(r"seed(\d+)", str(p)).group(1))
    )

    path_seed_dict = defaultdict(list)  # seed -> list of matching Paths
    for p in data_paths:
        seed_num = int(re.search(r"seed(\d+)", str(p)).group(1))
        path_seed_dict[seed_num].append(p)
    path_seed_dict = dict(path_seed_dict)

    data_columns = [
        "d_type",
        "noise_type",
        "seed",
        "d_char",
        "alloy",
        "length",
        "stress",
        "dydx_mean",
    ]
    data_frames = []
    # for dat in data_paths:
    for seed_num, paths in path_seed_dict.items():
        # sort by stress
        paths = sorted(
            paths, key=lambda p: int(re.search(r"(\d+)MPa", str(p)).group(1))
        )
        for path in paths:
            # extract stress value
            stress = int(re.search(r"(\d+)MPa", str(path)).group(1))

            # open data
            src = path
            tmp = None  # will hold TemporaryDirectory handle
            work_dir = None  # will point to the folder we actually use
            if src.is_dir():  # in uncompressed format
                work_dir = src
            elif src.suffixes == [".tar", ".gz"]:  # if .tar.gz format
                tmp = tempfile.TemporaryDirectory()
                with tarfile.open(src, "r:gz") as tf:
                    tf.extractall(tmp.name)
                work_dir = Path(tmp.name) / src.name.removesuffix(".tar.gz")
            else:
                raise FileNotFoundError("neither directory nor .tar.gz found")

            alloy = f"{config["alloy"]}.txt"
            matDir = work_dir / "inputFiles"
            # extract material info
            b_SI = readValFromMaterialFile(matDir, alloy, "b_SI")
            mu0_SI = readValFromMaterialFile(matDir, alloy, "mu0_SI")
            rho_SI = readValFromMaterialFile(matDir, alloy, "rho_SI")
            cs = np.sqrt(mu0_SI / rho_SI)  # shear wave speed
            convertTimeUnit = b_SI / cs  # [sec]
            convertMPaToMu = 1 / (mu0_SI * 10 ** (-6))

            # read deformation gradient tensor from poly file
            fTensor = np.zeros((3, 3))
            with open(work_dir / "inputFiles" / "polycrystal.txt", "r") as fh:
                mat = []
                for line in fh:
                    if "F=" in line:  # found header
                        mat.append(line.strip("F=").strip("\n").strip())
                        mat.append(fh.readline().strip("\n").strip())
                        mat.append(fh.readline().strip().split(";")[0])
                        break
                fTensor = np.loadtxt(mat, dtype=np.float64)

            V = np.linalg.det(fTensor)
            fDiag = np.diag(fTensor)
            A = fDiag[0] * fDiag[1]

            # read labels
            with open(work_dir / "F/F_labels.txt", "r") as label:
                fLabels = label.read()
            # remove empty element and store it as a list
            fLabels = np.array([x for x in fLabels.split("\n") if x], dtype=str)

            # last 33 element indexes
            lastCols = np.arange(-33, 0)
            # open F file
            fData = defaultdict(list)
            with open(work_dir / "F/F_0.txt", "r") as f:
                for line in f:
                    line = [float(x) for x in line.split(" ") if x and x != "\n"]
                    # parse the first 14 elements
                    for i in range(len(fLabels[:14])):
                        fData[str(fLabels[i])].append(line[i])
                    # parse the last 33 elements
                    for j in lastCols:
                        fData[str(fLabels[j])].append(line[j])

            xAxisData = np.array(fData[config["xAxisData"]])
            yAxisData = np.array(fData[config["yAxisData"]])

            # if there is a data error, skip the data
            if len(xAxisData) != len(yAxisData):
                continue

            # filter the duplicate data (minimization creates duplicates on runID)
            # Find duplicates and their indices
            unique_values, counts = np.unique(xAxisData, return_counts=True)
            duplicates = unique_values[counts > 1]
            # remove the first duplicate data entry
            dup_idx = np.where(xAxisData == duplicates)[0][0]
            xAxisData = np.delete(xAxisData, dup_idx)
            yAxisData = np.delete(yAxisData, dup_idx)

            # assume the first step data is minimization
            min_id_idx = 1
            # calculate the absolute plastic distortion
            abs_betaP = np.abs(yAxisData - yAxisData[min_id_idx])

            dydx = np.gradient(abs_betaP, xAxisData)
            dydx_mean = np.mean(dydx[-3:])

            data_table = [
                (
                    config["dislocType"],
                    config["noiseType"],
                    seed_num,
                    config["dislocChar"],
                    alloy.strip(".txt"),
                    config["dislocLength"],
                    stress,
                    dydx_mean,
                )
            ]
            print(data_table)
            data_frames.append(
                pd.DataFrame(data_table, columns=data_columns)
            )  # note the extra [] to make it 2-D

            if tmp:  # automatically deletes the tmp tree
                tmp.cleanup()

    # generate output dir
    os.makedirs(output_dir, exist_ok=True)

    # save plastic strain rate data
    df = pd.concat(data_frames, ignore_index=True)
    df.to_csv(output_dir / "all_rate_data.csv", index=False)

    # generate CRSS table
    seeds = df["seed"].unique()
    threshold = 1e-10  # set your cut-off

    data_frames = []
    data_columns = [
        "d_type",
        "noise_type",
        "seed",
        "d_char",
        "alloy",
        "length",
        "str_step_size",
        "crss",
    ]
    for seed in seeds:
        local_df = df[df["seed"] == seed]
        high_df = local_df[local_df["dydx_mean"] > threshold]  # 1. apply threshold
        crss_above_row = high_df.loc[high_df["dydx_mean"].idxmin()]
        crss_above = crss_above_row["stress"]

        # rows whose stress is below the CRSS
        below = local_df[local_df["stress"] < crss_above]
        if not len(below):
            crss = crss_above
            crss_search_step = crss_above
        else:
            # pick the one with the highest such stress
            entry_below = below.loc[below["stress"].idxmax(), "stress"]
            crss_search_step = crss_above - entry_below
            crss = crss_above - (crss_search_step / 2)
        # pick the one with the highest such stress
        #entry_below = below.loc[below["stress"].idxmax(), "stress"]
        #crss_search_step = crss_above - entry_below
        #crss = crss_above - (crss_search_step / 2)
        new_entry = [
            (
                crss_above_row["d_type"],
                crss_above_row["noise_type"],
                crss_above_row["seed"],
                crss_above_row["d_char"],
                crss_above_row["alloy"],
                crss_above_row["length"],
                crss_search_step,
                crss,
            )
        ]
        data_frames.append(pd.DataFrame(new_entry, columns=data_columns))

    df_new = pd.concat(data_frames, ignore_index=True)
    df_new = df_new.sort_values(by=["seed"])
    df_new.to_csv(output_dir / "all_crss_data.csv", index=False)


if __name__ == "__main__":
    main()
