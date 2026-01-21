import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from scipy import stats
from pathlib import Path
from matplotlib import rcParams


# Configure global plot settings (applies to all figures)
rcParams.update(
    {
        "figure.dpi": 200,
        "figure.autolayout": True,  # Prevent label clipping
        "axes.grid": True,
        "grid.alpha": 0.9,
        #'text.usetex': True,   # Enable LaTeX
        "text.usetex": False,  # Enable LaTeX
        "font.size": 20,  # Default font size for text
        "mathtext.fontset": "stix",  # Use STIX font for math text
        "font.family": "serif",  # Use serif font (matches LaTeX default)
        #'font.weight': 'bold',          # Base font weight
        #'text.latex.preamble': r'\usepackage{amsmath} \boldmath',  # Additional packages
        #'text.latex.preamble': r'\usepackage{amsmath}',  # Additional packages
    }
)

data_dir = Path("data/")
data_fnames = [
    "ex1_s01_mg5.jsonl",
    "ex1_s01_mg10.jsonl",
    "ex1_s01_mg15.jsonl",
    #"ex1_s23_mg5.jsonl",
    #"ex1_s23_mg10.jsonl",
    #"ex1_s23_mg15.jsonl",
    #"ex2_s01_mg5.jsonl",
    #"ex2_s01_mg10.jsonl",
    #"ex2_s01_mg15.jsonl",
    #"ex2_s23_mg5.jsonl",
    #"ex2_s23_mg10.jsonl",
    #"ex2_s23_mg15.jsonl",
]
out_fname = Path("figures/ex1_s01.png")
#out_fname = Path("figures/ex1_s23.png")
#out_fname = Path("figures/ex2_s01.png")
#out_fname = Path("figures/ex2_s23.png")

fig, ax = plt.subplots(1, 1, figsize=(12, 6), dpi=200)
# Number of colors
n_colors = 3
colors = sns.color_palette("tab10")

burger = 2.86 # angstrom
ang_to_nm = 1e-1
for idx, data_fname in enumerate(data_fnames):
    alloy = re.search(r"mg(\d+)", str(data_fname)).group(0)
    df = pd.read_json(data_dir / data_fname, lines=True)

    sep_data = df["relaxed_sep_dist"]*burger*ang_to_nm
    n, bins, patches = ax.hist(
        sep_data,
        bins=30,
        density=True,
        alpha=0.4,
        color=colors[idx],
        edgecolor=None,
        linewidth=0.5,
    )
    # Calculate Gaussian fit parameters
    mu, sigma = np.mean(sep_data), np.std(sep_data)
    # Generate x values for the Gaussian curve
    x = np.linspace(sep_data.min(), sep_data.max(), 1000)
    # Calculate Gaussian probability density function
    gaussian_fit = stats.norm.pdf(x, mu, sigma)
    # Plot the Gaussian fit
    ax.plot(
        x,
        gaussian_fit,
        linewidth=3,
        color=colors[idx],
        label=f"{alloy} $\\mu$={mu:.4f}, $\\sigma$={sigma:.4f}",
    )

# Add labels and legend
ax.set_xlabel("d [nm]")
ax.set_ylabel("Probability Density")
# ax.set_title('')
ax.legend(
    loc="center left",  # reference point inside the bbox
    bbox_to_anchor=(1.02, 0.5),  # (x, y) in *axes* coords
    frameon=False,
)  # no box, no background
plt.tight_layout()
#plt.show()
# save
os.makedirs(out_fname.parent, exist_ok=True)
fig.savefig(out_fname, transparent=True)
# plt.show()
