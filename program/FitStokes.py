import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np
import sys

def plot_fitstokes(data_path, out_path, d_max=500):
    
    df = pd.read_csv(data_path, sep=r'\s+', header=None, comment='!', engine='python')
    
    labels = ['I', 'Q/I', 'U/I', 'V/I']
    y_cols = [(2, 4, 3), (5, 7, 6), (8, 10, 9), (11, 13, 12)]  # zero-based indexing
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), sharex='col')
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
    for i, (label, cols) in enumerate(zip(labels, y_cols)):
        d1 = df[cols[0]]
        d2 = df[cols[1]]
        d3 = df[cols[2]]
        x = df[0]

        # --- Top row: data + model
        ax = axes[0, i]
        ax.errorbar(x, d1, fmt='x', capsize=3, capthick=1, yerr=d2, label='DATA', color='red')
        ax.scatter(x, d3, label='MGMR3D', color='blue')
        ax.set_title(label)
        ax.set_xlim(0, d_max)
        ax.grid(True, which='both', linestyle=':', linewidth=0.5)
        ax.legend(loc='upper right')

        if label == 'I':
            ax.set_ylim(0, None)
        else:
            ax.set_ylim(-1, 1)

        # --- Bottom row: pull distribution
        ax2 = axes[1, i]
        with np.errstate(divide='ignore', invalid='ignore'):
            pull = np.where(d2 != 0, (d3 - d1) / d2, 0)
        ax2.plot(x, pull, 'k.')
        ax2.axhline(0, color='black', linestyle='-')
        ax2.set_xlim(0, d_max)
        ax2.set_ylim(-2, 2)
        ax2.set_xlabel("d [m]")
        ax2.grid(True, which='both', linestyle=':')

        if i == 0:
            ax2.set_ylabel("Pull")

    fig.suptitle("FitStokes Equivalent", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(f"{out_path}/FitStokes.jpg")
    plt.show()

if __name__ == "__main__":
    
    plot_fitstokes(sys.argv[1], sys.argv[2])