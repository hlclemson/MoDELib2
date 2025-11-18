# standard lib
import re
import sys
import json
import math
import shutil
import tarfile
import pathlib
import tempfile

# 3rd party lib
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path
from matplotlib import rcParams
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def ci95(data):
    """return mean and 95 % margin for a 1-D array"""
    n = len(data)
    m = np.mean(data)
    se = stats.sem(data)
    h = se * stats.t.ppf((1 + 0.95) / 2, n - 1)
    return m, h


# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 200,
        "figure.autolayout": True,  # Prevent label clipping
        "axes.grid": True,
        "grid.alpha": 0.6,
        "text.usetex": False,
        "font.size": 16,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
    }
)


def main():
    ref_data = Path("reference_data/all_crss_data_final.csv")
    ref_df = pd.read_csv(ref_data)
    ref_df = ref_df.rename(columns={"d_char": "d_type", "d_type": "d_char"})


    step_size_limit = 30
    ref_edge_almg5_crss = ref_df[
        (ref_df["length"] == 215)
        & (ref_df["noise_type"] == "SF-noise")
        & (ref_df["d_char"] == "edge")
        & (ref_df["alloy"] == "AlMg5")
        & (ref_df["str_step_size"] < step_size_limit)
    ]
    ref_edge_almg10_crss = ref_df[
        (ref_df["length"] == 215)
        & (ref_df["noise_type"] == "SF-noise")
        & (ref_df["d_char"] == "edge")
        & (ref_df["alloy"] == "AlMg10")
        & (ref_df["str_step_size"] < step_size_limit)
    ]
    ref_edge_almg15_crss = ref_df[
        (ref_df["length"] == 215)
        & (ref_df["noise_type"] == "SF-noise")
        & (ref_df["d_char"] == "edge")
        & (ref_df["alloy"] == "AlMg15")
        & (ref_df["str_step_size"] < step_size_limit)
    ]

    # reference
    ref_edge_almg5_mu, ref_edge_almg5_err = zip(ci95(ref_edge_almg5_crss["crss"]))
    ref_edge_almg10_mu, ref_edge_almg10_err = zip(ci95(ref_edge_almg10_crss["crss"]))
    ref_edge_almg15_mu, ref_edge_almg15_err = zip(ci95(ref_edge_almg15_crss["crss"]))
    ref_x = np.array([5,10,15])
    ref_y = np.array([ref_edge_almg5_mu, ref_edge_almg10_mu, ref_edge_almg15_mu]).ravel()
    ref_y_err = np.array([ref_edge_almg5_err, ref_edge_almg10_err, ref_edge_almg15_err]).ravel()
    plt.errorbar(ref_x, ref_y, yerr=ref_y_err, fmt='o', capsize=5, label="sID0-sID1")

    # 2-3 system
    data5_23 = Path("output/partial_orient_615AA_crss_2-3_almg5_total/all_crss_data.csv")
    data10_23 = Path("output/partial_orient_615AA_crss_2-3_almg10_total/all_crss_data.csv")
    data15_23 = Path("output/partial_orient_615AA_crss_2-3_almg15_total/all_crss_data.csv")
    df5_23 = pd.read_csv(data5_23)
    df10_23 = pd.read_csv(data10_23)
    df15_23 = pd.read_csv(data15_23)
    edge_almg5_mu_23, edge_almg5_err_23 = zip(ci95(df5_23["crss"]))
    edge_almg10_mu_23, edge_almg10_err_23 = zip(ci95(df10_23["crss"]))
    edge_almg15_mu_23, edge_almg15_err_23 = zip(ci95(df15_23["crss"]))
    x_23 = np.array([5,10,15])
    y_23 = np.array([edge_almg5_mu_23, edge_almg10_mu_23, edge_almg15_mu_23]).ravel()
    y_23_err = np.array([edge_almg5_err_23, edge_almg10_err_23, edge_almg15_err_23]).ravel()
    plt.errorbar(x_23, y_23, yerr=y_23_err, fmt='o', capsize=5, label="sID2-sID3")

    # 4-5 system
    data5_45 = Path("output/partial_orient_615AA_crss_4-5_almg5_total/all_crss_data.csv")
    data10_45 = Path("output/partial_orient_615AA_crss_4-5_almg10_total/all_crss_data.csv")
    data15_45 = Path("output/partial_orient_615AA_crss_4-5_almg15_total/all_crss_data.csv")
    df5_45 = pd.read_csv(data5_45)
    df10_45 = pd.read_csv(data10_45)
    df15_45 = pd.read_csv(data15_45)
    edge_almg5_mu_45, edge_almg5_err_45 = zip(ci95(df5_45["crss"]))
    edge_almg10_mu_45, edge_almg10_err_45 = zip(ci95(df10_45["crss"]))
    edge_almg15_mu_45, edge_almg15_err_45 = zip(ci95(df15_45["crss"]))
    x_45 = np.array([5,10,15])
    y_45 = np.array([edge_almg5_mu_45, edge_almg10_mu_45, edge_almg15_mu_45]).ravel()
    y_45_err = np.array([edge_almg5_err_45, edge_almg10_err_45, edge_almg15_err_45]).ravel()
    plt.errorbar(x_45, y_45, yerr=y_45_err, fmt='o', capsize=5, label="sID4-sID5")
    plt.legend()

    plt.title("Edge Dislocation Partial, 215b")
    plt.xlabel("Mg Concentration [%]")
    plt.ylabel("CRSS [MPa]")

    plt.savefig("crss_comparison.png",transparent=True)


    return 0


if __name__ == "__main__":
    main()
