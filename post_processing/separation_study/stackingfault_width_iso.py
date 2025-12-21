import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
from dataclasses import dataclass
import matplotlib.cm as cm
from pathlib import Path
import pandas as pd

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

def sfe_width_isotropic(mu, nu, b_mag, theta, isf):
    return (mu*b_mag**2)/(8*np.pi*isf)*((2-nu)/(1-nu))*(1-(2*nu*np.cos(np.deg2rad(2*theta))/(2-nu)))


fig, ax = plt.subplots(1, 1, figsize=(10, 6))
ths = np.linspace(0, 90)
aa_to_nm = 1e+9
#aa_to_nm = 1
b_p = np.array([-2,-1,1])*almg5_SI.alat/6
ax.plot(ths, sfe_width_isotropic(almg5_SI.mu, almg5_SI.nu, np.linalg.norm(b_p), ths, almg5_SI.isf)*aa_to_nm, linewidth=2, label="AlMg5")
ax.plot(ths, sfe_width_isotropic(almg10_SI.mu, almg10_SI.nu, np.linalg.norm(b_p), ths, almg10_SI.isf)*aa_to_nm, linewidth=2, label="AlMg10")
ax.plot(ths, sfe_width_isotropic(almg15_SI.mu, almg15_SI.nu, np.linalg.norm(b_p), ths, almg15_SI.isf)*aa_to_nm, linewidth=2, label="AlMg15")
#plt.tight_layout()
#fig.savefig("separation_iso.png")

data_out_dir = Path("generated_data")
# total_separation_data = pd.read_csv(data_out_dir/"separation_distance_total_data.csv")
total_separation_data = pd.read_pickle(
    data_out_dir / "separation_distance_total_data.pkl"
)

exit_ids = total_separation_data["exitID"].unique()
slip_ids = total_separation_data["sIDs"].unique()
alloys = total_separation_data["Mg"].unique()
leading_partial_ids = total_separation_data["leading_partial_id"].unique()
n_colors = 6
colors = cm.tab10(np.linspace(0, 1, n_colors))
labels = {
    "e1_s0-1_l0": "(a)",
    "e2_s0-1_l0": "(b)",
    "e1_s2-3_l2": "(c)",
    "e1_s2-3_l3": "(d)",
    "e2_s2-3_l3": "(e)",
    "e2_s2-3_l2": "(f)",
}
color_per_type = {key: colors[i] for i, key in enumerate(labels.keys())}

b_SI=2.86e-10 # m
# plot all
angles_per_type = {
    "e1_s0-1_l0": 90,
    "e2_s0-1_l0": 0,
    "e1_s2-3_l2": 30,
    "e1_s2-3_l3": 30,
    "e2_s2-3_l3": 60,
    "e2_s2-3_l2": 60,
}

#handles = []
#fig, ax = plt.subplots(1, 1, figsize=(9, 6))
for exit_id in exit_ids:
    for slip_id in slip_ids:
        for leading_partial_id in leading_partial_ids:
            data = total_separation_data[
                (total_separation_data["sIDs"] == slip_id)
                & (total_separation_data["exitID"] == exit_id)
                & (total_separation_data["stress"] == 1)
                & (
                    total_separation_data["leading_partial_id"]
                    == leading_partial_id
                )
            ]
            if data.empty:
                continue
            mg = data["Mg"].astype(int)
            # plot only 5 Mg
            data = data[data["Mg"]=="5"]
            label_type = f"e{exit_id}_s{slip_id}_l{leading_partial_id}"
            angle = angles_per_type[label_type]
            #angles = [angle, angle, angle]
            angles = angle
            hd = plt.errorbar(
                angles,
                data["separation_mean"]*b_SI*aa_to_nm,
                yerr=data["separation_ci95"]*b_SI*aa_to_nm,
                fmt="o",
                capsize=8,
                markersize=8,
                color=color_per_type[label_type],
                label=labels[label_type],
            )
            #handles.append(hd)
# reorder:
#labels = sorted([h.get_label() for h in handles])
#print(labels)
#exit()
#plt.legend(handles, labels)

#ax.set_xticks(alloys.astype(int))
ax.set_xlabel("$\\theta$")
ax.set_ylabel("Separation Width [nm]")
# legend fully outside the axes
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0., frameon=False)
fig.tight_layout() # keeps everything visible when saving/showing
fig_dir = Path("figures")
os.makedirs(fig_dir, exist_ok=True)
fig.savefig(fig_dir / f"relaxation_with_iso_width.png", transparent=True)
