import re
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from scipy import stats
from pathlib import Path
from matplotlib import rcParams
from dataclasses import dataclass

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

@dataclass
class almg5_SI:
    mu: np.float64 = 28.595e+9
    alat: np.float64 = 4.069e-10
    nu: np.float64 = 0.3406
    isf: np.float64 = 110e-3

@dataclass
class almg10_SI:
    mu: np.float64 = 26.783697e9
    alat: np.float64 = 4.093e-10
    nu: np.float64 = 0.342741
    isf: np.float64 = 102e-3

@dataclass
class almg15_SI:
    mu: np.float64 = 24.971966e9
    alat: np.float64 = 4.119e-10
    nu: np.float64 = 0.344882
    isf: np.float64 = 98e-3


def ci95(data):
    """return mean and 95 % margin for a 1-D array"""
    n = len(data)
    m = np.mean(data)
    se = stats.sem(data)
    h = se * stats.t.ppf((1 + 0.95) / 2, n - 1)
    return m, h


def sfe_width_isotropic(mu, nu, b_mag, theta, isf):
    return (mu*b_mag**2)/(8*np.pi*isf)*((2-nu)/(1-nu))*(1-(2*nu*np.cos(np.deg2rad(2*theta))/(2-nu)))


# Number of colors
n_colors = 3
colors = sns.color_palette("tab10")

# fig filename
out_fname = Path("figures/relaxation_with_iso_width_comp.png")

fig, ax = plt.subplots(1, 1, figsize=(11, 6))
ths = np.linspace(0, 90)
aa_to_nm = 1e+9
#aa_to_nm = 1
b_p = np.array([-2,-1,1])*almg5_SI.alat/6
ax.plot(ths, sfe_width_isotropic(almg5_SI.mu, almg5_SI.nu, np.linalg.norm(b_p), ths, almg5_SI.isf)*aa_to_nm, color=colors[0], linewidth=2, label="AlMg5")
ax.plot(ths, sfe_width_isotropic(almg10_SI.mu, almg10_SI.nu, np.linalg.norm(b_p), ths, almg10_SI.isf)*aa_to_nm, color=colors[1], linewidth=2, label="AlMg10")
ax.plot(ths, sfe_width_isotropic(almg15_SI.mu, almg15_SI.nu, np.linalg.norm(b_p), ths, almg15_SI.isf)*aa_to_nm, color=colors[2], linewidth=2, label="AlMg15")

data_dir = Path("data/")
data_fnames = [
    "ex1_s01_mg5.jsonl",
    "ex1_s01_mg10.jsonl",
    "ex1_s01_mg15.jsonl",
    "ex1_s23_mg5.jsonl",
    "ex1_s23_mg10.jsonl",
    "ex1_s23_mg15.jsonl",
    "ex2_s01_mg5.jsonl",
    "ex2_s01_mg10.jsonl",
    "ex2_s01_mg15.jsonl",
    "ex2_s23_mg5.jsonl",
    "ex2_s23_mg10.jsonl",
    "ex2_s23_mg15.jsonl",
]

labels = {
    "ex1_s01_l0": "(a)",
    "ex2_s01_l0": "(b)",
    "ex1_s23_l2": "(c)",
    "ex1_s23_l3": "(d)",
    "ex2_s23_l3": "(e)",
    "ex2_s23_l2": "(f)",
}

# plot all
#angles_per_type = {
#    "e1_s0-1_l0": 90,
#    "e2_s0-1_l0": 0,
#    "e1_s2-3_l2": 30,
#    "e1_s2-3_l3": 30,
#    "e2_s2-3_l3": 60,
#    "e2_s2-3_l2": 60,
#}
angles_per_type = {
    "ex1_s01": 90,
    "ex2_s01": 0,
    "ex1_s23": 30,
    "ex2_s23": 60,
}

burger = 2.86 # angstrom
ang_to_nm = 1e-1
color_per_type = {
    "mg5": colors[0],
    "mg10": colors[1],
    "mg15": colors[2],
}
for idx, data_fname in enumerate(data_fnames):
    alloy = re.search(r"mg(\d+)", str(data_fname)).group(0)
    id = re.search(r"ex(\d+)_s(\d+)", str(data_fname)).group(0)
    df = pd.read_json(data_dir / data_fname, lines=True)

    sep_data = df["relaxed_sep_dist"]*burger*ang_to_nm
    mean, error95 = ci95(sep_data)
    angle = angles_per_type[id]
    ax.errorbar(
        angle,
        mean,
        yerr=error95,
        fmt="o",
        capsize=8,
        markersize=8,
        color=color_per_type[alloy],
        #label=labels[label_type],
    )

#ax.set_xticks(alloys.astype(int))
ax.set_xlabel("Angle [$\\theta$]")
ax.set_ylabel("Separation Width [nm]")
ax.set_ylim(0.4,1.05)
ax.set_xticks([0,30,60,90])
# legend fully outside the axes
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0., frameon=False)
fig.tight_layout() # keeps everything visible when saving/showing
#fig_dir = Path("figures")
os.makedirs(out_fname.parent, exist_ok=True)
fig.savefig(out_fname, transparent=True)
