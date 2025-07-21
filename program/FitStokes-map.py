from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np
import sys

def plot_stokesmap(data_dir, out_path, stokes_type="I"):
    """
    Plot a Stokes map from data in polar coordinates.
    
    Parameters:
        data_dir (str): Directory containing FitResult.dat and grid_StX.grd files
        stokes_type (str): Which Stokes parameter to plot: 'I', 'Q/I', 'U/I', 'V/I', '\\Psi'
        save (bool): If True, saves the plot as PNG
    """
    # File paths
    fitresult_path = f"{data_dir}FitResult.dat"
    grid_path = f"{data_dir}grid_StI.grd"  # Always read from I for layout/meta

    # Read FitResult.dat
    df = pd.read_csv(fitresult_path, sep=r'\s+', header=None, comment='!', engine='python')
    R = df[0].astype(float)
    theta = df[1].astype(float)
    x = R * np.cos(theta)
    y = R * np.sin(theta)

    # Determine bounds
    LdistX, UdistX = x.min() - 50, x.max() + 50
    LdistY, UdistY = y.min() - 100, y.max() + 50

    # Read grid meta info
    with open(grid_path, 'r') as f:
        line = f.readline().strip().split()
        N_grid, min_grid, max_grid, StI_max, d_gdid, Norm_I, corex, corey = map(float, line[:8])
    cc_min = 0
    cc_max = StI_max * (Norm_I if stokes_type == "I" else 1.0)
    cc_Stp = (cc_max - cc_min) / 5

    # Map stokes_type to correct column index
    stokes_map = {
        "I": 3,
        "Q/I": 6,
        "U/I": 9,
        "V/I": 12,
        "\\Psi": 16,
    }
    if stokes_type not in stokes_map:
        raise ValueError(f"Unsupported stokes_type: {stokes_type}")
    stokes_idx = stokes_map[stokes_type]
    z = df[stokes_idx].astype(float)

    # Prepare plot
    fig, ax = plt.subplots(figsize=(9, 9))
    norm = Normalize(vmin=cc_min, vmax=cc_max)
    cmap = plt.get_cmap("plasma")

    sc = ax.scatter(x, y, c=z, cmap=cmap, norm=norm, s=10)
    ax.scatter(corex, corey, marker='x', color='black', s=100)
    ax.set_xlim(min_grid, max_grid)
    ax.set_ylim(min_grid, max_grid)
    ax.axis('off')

    # Colorbar
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)
    cbar.set_label(f"Stokes {stokes_type}")

    plt.tight_layout()
    fig.savefig(f"{out_path}Stokes_{stokes_type.replace('/', '')}.jpg") #FitStokes-map
    plt.show()

def plot_stokes_parameters(data_dir, out_path):
    """
    Plot all Stokes parameters in a 2x3 grid with pcolormesh and scatter overlays.

    Parameters:
        data_dir (str): Directory containing FitResult.dat and grid_StI.grd
        out_path (str): Path prefix to save output image
    """
    # Load grid metadata
    grid_path = f"{data_dir}/grid_StI.grd"
    with open(grid_path, 'r') as f:
        line = f.readline().strip().split()
        _, min_grid, max_grid, StI_max, _, Norm_I, corex, corey = map(float, line[:8])

    # Load data
    df = pd.read_csv(f"{data_dir}/FitResult.dat", sep=r'\s+', header=None, comment='!', engine='python')
    R = df[0].values
    theta = df[1].values
    X = R * np.cos(theta)
    Y = R * np.sin(theta)

    # Prepare grid
    n_grid = 500
    xi = np.linspace(min_grid, max_grid, n_grid)
    yi = np.linspace(min_grid, max_grid, n_grid)
    Xi, Yi = np.meshgrid(xi, yi)

    # Stokes parameters
    stokes_map = {
        "I": 3,
        "Q/I": 6,
        "U/I": 9,
        "V/I": 12,
        "Ψ": 16,
    }
    cmap_dict = {
        "I": "plasma",
        "Q/I": "bwr",
        "U/I": "bwr",
        "V/I": "bwr",
        "Ψ": "twilight_shifted",
    }
    norm_limits = {
        "I": (0, StI_max * Norm_I),
        "Q/I": (-1, 1),
        "U/I": (-1, 1),
        "V/I": (-1, 1),
        "Ψ": (-90, 90)
    }

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))

    for ax, (label, col_idx) in zip(axs.flat, stokes_map.items()):
        Z = df[col_idx].values
        cmap = cmap_dict[label]
        vmin, vmax = norm_limits[label]
        norm = Normalize(vmin=vmin, vmax=vmax)

        # Interpolate Z to grid
        Zi = griddata((X, Y), Z, (Xi, Yi), method='cubic')

        # Plot interpolated background with pcolormesh
        pc = ax.pcolormesh(Xi, Yi, Zi, shading='auto', cmap=cmap, norm=norm)

        # Overlay scatter points
        ax.scatter(X, Y, c=Z, cmap=cmap, norm=norm, edgecolors='k', linewidths=0.2, s=25, zorder=3)

        # Core location
        ax.plot(corex, corey, 'kx', markersize=8)

        # Axis formatting
        ax.set_aspect('equal')
        ax.grid()
        ax.set_xlim(min_grid, max_grid)
        ax.set_ylim(min_grid, max_grid)
        ax.set_xlabel('[m]')
        ax.set_ylabel('[m]')
        ax.set_title(f"Stokes {label}")

        # Colorbar
        cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax)

    axs[1, 2].axis('off')  # Hide unused subplot
    fig.suptitle("Stokes Parameters Map", fontsize=12)
    plt.tight_layout()
    fig.savefig(f"{out_path}/FitStokes-map.jpg")
    plt.show()

if __name__ == "__main__":
    plot_stokes_parameters(sys.argv[1], sys.argv[2])


