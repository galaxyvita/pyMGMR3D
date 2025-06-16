import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import numpy as np
import sys

def plot_sh_current(data_dir, out_path):
    """
    Plots the sh_Current.dat file with J_x, J_y, and J_Q*10.

    Parameters:
        data_dir (str): Directory where sh_Current.dat is located.
    """
    datfile = Path(data_dir) / "sh_Current.dat"

    # Read data
    try:
        df = pd.read_csv(datfile, delim_whitespace=True, header=None, comment="!")
    except Exception as e:
        print(f"Failed to read {datfile}: {e}")
        return

    # Extract columns
    x = df[0]
    J_x = df[3]
    J_y = df[4]
    J_Q = df[5] * 10  # scale by 10 as in GLE

    # Set up plot
    fig, ax = plt.subplots(figsize=(10, 5))

    ax.plot(x, J_x, label="Jₓ", color="green", linestyle='-')
    ax.plot(x, J_y, label="Jᵧ", color="red", linestyle='--')
    ax.plot(x, J_Q, label="J_Q × 10", color="blue", linestyle='-.')

    ax.set_xlim(0, 20)
    ax.set_xlabel("D to ground [km]")
    ax.set_ylabel("Current (a.u.)")  # Add appropriate unit if known
    ax.grid()

    ax.legend(loc="upper right", frameon=False)

    plt.tight_layout()
    fig.savefig(f"{out_path}sh_Current.jpg")
    plt.show()

if __name__ == "__main__":
    
    plot_sh_current(sys.argv[1], sys.argv[2])
