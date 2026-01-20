import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import stats
from pathlib import Path
from matplotlib import rcParams


# Configure global plot settings (applies to all figures)
rcParams.update({
    'figure.dpi': 200,
    'figure.autolayout': True,  # Prevent label clipping
    'axes.grid': True,
    'grid.alpha': 0.6,
    #'text.usetex': True,   # Enable LaTeX
    'text.usetex': False,   # Enable LaTeX
    'font.size': 20,     # Default font size for text
    'mathtext.fontset': 'stix',  # Use STIX font for math text
    'font.family': 'serif',  # Use serif font (matches LaTeX default)
    #'font.weight': 'bold',          # Base font weight
    #'text.latex.preamble': r'\usepackage{amsmath} \boldmath',  # Additional packages
    #'text.latex.preamble': r'\usepackage{amsmath}',  # Additional packages
})


data = Path("data/ex1_s23_mg5.jsonl")
out_fname = Path("figures/ex1_s23_mg5_hist.png")
df = pd.read_json(data, lines=True)

sep_data = df["relaxed_sep_dist"]
fig, ax = plt.subplots(1, 1, figsize=(14, 6), dpi=200)
n, bins, patches = ax.hist(sep_data, bins=30, density=True, alpha=0.7, color='skyblue', edgecolor='black', linewidth=0.5)
# Calculate Gaussian fit parameters
mu, sigma = np.mean(sep_data), np.std(sep_data)
# Generate x values for the Gaussian curve
x = np.linspace(sep_data.min(), sep_data.max(), 1000)
# Calculate Gaussian probability density function
gaussian_fit = stats.norm.pdf(x, mu, sigma)
# Plot the Gaussian fit
ax.plot(x, gaussian_fit, 'r-', linewidth=2, label=f'$\\mu$={mu:.2f}, $\\sigma$={sigma:.2f}')

# Add labels and legend
ax.set_xlabel('d')
ax.set_ylabel('Probability Density')
#ax.set_title('')
ax.legend(
    loc="center left",  # reference point inside the bbox
    bbox_to_anchor=(1.02, 0.5),  # (x, y) in *axes* coords
    frameon=False,
)  # no box, no background
plt.tight_layout()

# save
os.makedirs(out_fname.parent, exist_ok=True)
fig.savefig(out_fname, transparent=True)
#plt.show()

