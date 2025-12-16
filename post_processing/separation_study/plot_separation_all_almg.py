# standard lib
import enum
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
from matplotlib import rcParams
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 200,
        # "figure.autolayout": True,  # Prevent label clipping
        "axes.grid": True,
        "grid.alpha": 0.8,
        "text.usetex": False,
        "font.size": 25,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
    }
)


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
    fig_dir = Path("figures")
    os.makedirs(fig_dir, exist_ok=True)

    data_out_dir = Path("generated_data")
    # total_separation_data = pd.read_csv(data_out_dir/"separation_distance_total_data.csv")
    total_separation_data = pd.read_pickle(
        data_out_dir / "separation_distance_total_data.pkl"
    )

    exit_ids = total_separation_data["exitID"].unique()
    slip_ids = total_separation_data["sIDs"].unique()
    alloys = total_separation_data["Mg"].unique()
    leading_partial_ids = total_separation_data["leading_partial_id"].unique()
    # n_colors = 10
    # colors = cm.tab10(np.linspace(0, 1, n_colors))

    # plot all
    for exit_id in exit_ids:
        for slip_id in slip_ids:
            for leading_partial_id in leading_partial_ids:
                for i, alloy in enumerate(alloys):
                    data = total_separation_data[
                        (total_separation_data["Mg"] == alloy)
                        & (total_separation_data["sIDs"] == slip_id)
                        & (total_separation_data["exitID"] == exit_id)
                        & (
                            total_separation_data["leading_partial_id"]
                            == leading_partial_id
                        )
                    ]
                    if data.empty:
                        continue
                    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
                    crss = data["CRSS"].unique()
                    ax.axvline(crss, linestyle="--", color="r")
                    # stresses = data['stress']
                    plt.errorbar(
                        data["stress"],
                        data["separation_mean"],
                        yerr=data["separation_ci95"],
                        fmt="o",
                        color="k",
                        capsize=5,
                        markersize=8,
                        label=f"AlMg{alloy}",
                    )
                    #ax.legend()
                    ax.set_xlabel("Stress [MPa]")
                    # ax.xlim(-10, 150)
                    ax.set_ylabel("Distance [b]")
                    plt.tight_layout()
                    fig.savefig(fig_dir/f"exitID{exit_id}_AlMg{alloy}_sID{slip_id}_lead{leading_partial_id}.png", transparent=True)
                    plt.close("all")


if __name__ == "__main__":
    main()
