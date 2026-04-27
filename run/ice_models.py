"""
ice_models.py
=============
Ice density/refractivity profile generator for pyMGMR3D in-ice simulations.
Writes current_ice.dat, which is read by Initialize_Ice_Layer in Ice.f90.

Usage (standalone):
    python ice_models.py --input SIM000001.in

Usage (imported by pyMGMR3D.sh wrapper or another Python script):
    from ice_models import run_from_input_file
    run_from_input_file("SIM000001.in")

Column format of current_ice.dat:
    depth[m]  density[g/cm3]  xi=n-1  X_cumul_vertical[g/cm2]

Ice models:
    0 : no ice (file not created)
    1 : constant bulk ice (rho=0.917, n=1.78)
    2 : South Pole single-exponential firn
    3 : Greenland single-exponential firn

Density -> refractive index:  n(z) = 1 + 0.86 * rho(z)   [rho in g/cm3]
"""

import re
import numpy as np
import argparse
import sys
import os


# ---------------------------------------------------------------------------
# Ice model parameters
# ---------------------------------------------------------------------------

ICE_MODELS = {
    1: {
        "name": "Constant bulk ice",
        "rho_ice":  0.917,
        "rho_surf": 0.917,
        "scale_depth": None,      # constant density, no exponential
    },
    2: {
        "name": "South Pole single-exponential (Schytt/RICE/ARA)",
        "rho_ice":  0.917,        # g/cm3, asymptotic deep-ice density
        "rho_surf": 0.359,        # g/cm3, surface snow density at South Pole
        "scale_depth": 37.0,      # m, firn compaction scale depth
    },
    3: {
        "name": "Greenland single-exponential (RNO-G Summit Station)",
        "rho_ice":  0.917,
        "rho_surf": 0.400,        # g/cm3, surface snow at Summit Station
        "scale_depth": 26.0,      # m
    },
}

N_COEFF = 0.86   # n(z) = 1 + N_COEFF * rho(z)


# ---------------------------------------------------------------------------
# Input file reader  (matches style of get_val used for atmosphere)
# ---------------------------------------------------------------------------

def get_val(text, key):
    """
    Extract a value from an MGMR3D .in input file.
    Handles lines like:  key = value,  or  key = value ! comment
    Ignores inline comments after '!' and trailing commas.
    Matching is case-insensitive.
    """
    pattern = rf"^\s*{key}\s*=\s*([^,! \n]+)"
    match = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
    if match:
        return match.group(1).strip().replace("'", "").replace('"', "")
    return None


def read_ice_params(input_file):
    """
    Read ice-related parameters from a pyMGMR3D .in input file.

    Returns a dict with keys:
        ice_model  (int)   : 0=none, 1=constant, 2=S.Pole, 3=Greenland
        X_max_ice (float) : maximum ice depth to simulate [m]
        step_ice      (float) : depth step size [m]
        observer_z    (float) : antenna depth below surface [m], 0=surface, <0=in ice
        Zen_sh        (float) : shower zenith angle [deg], used for slant-depth conversion
    """
    with open(input_file, "r") as f:
        text = f.read()

    # Defaults (must match Fortran defaults in SetParams)
    params = {
        "ice_model":  int(get_val(text, "ice_model")  or 0),
        "X_max_ice": float(get_val(text, "X_max_ice") or 500.0),
        "step_ice":      float(get_val(text, "step_ice")      or 5.0),
        "observer_z":    float(get_val(text, "observer_z")    or 0.0),
        "Zen_sh":        float(get_val(text, "Zen_sh")        or 0.0),
    }
    return params

# ---------------------------------------------------------------------------
# Profile generation
# ---------------------------------------------------------------------------

def density_profile(model_id, depths):
    """Return density [g/cm3] as a function of depth [m]."""
    p = ICE_MODELS[model_id]
    if p["scale_depth"] is None:
        return np.full_like(depths, p["rho_ice"])
    C = 1.0 / p["scale_depth"]
    return p["rho_ice"] - (p["rho_ice"] - p["rho_surf"]) * np.exp(-C * depths)


def generate_ice_profile(ice_model_id, depth_max_m=500.0, dz_m=5.0,
                         output_file="current_ice.dat", zenith_deg=0.0,
                         verbose=True):
    """
    Generate current_ice.dat for a given ice model.

    Parameters
    ----------
    ice_model_id  : int    — 0=none, 1=constant, 2=S.Pole, 3=Greenland
    depth_max_m   : float  — maximum ice depth [m]
    dz_m          : float  — depth step size [m]
    output_file   : str    — output filename
    zenith_deg    : float  — shower zenith angle [deg] (for informational output only)
    verbose       : bool

    Returns
    -------
    dict with arrays: depth, density, refractivity, X_cumul_vertical
    (empty dict if ice_model_id == 0)
    """
    if ice_model_id == 0:
        if verbose:
            print("ice_model=0: no ice, current_ice.dat not created.")
        return {}

    if ice_model_id not in ICE_MODELS:
        raise ValueError(
            f"Unknown ice_model={ice_model_id}. "
            f"Valid values: {list(ICE_MODELS.keys())}"
        )

    depths  = np.arange(0.0, depth_max_m + dz_m * 0.5, dz_m)
    rho     = density_profile(ice_model_id, depths)
    n       = 1.0 + N_COEFF * rho
    xi      = n - 1.0   # refractivity; this is what the Fortran uses as IceRefrac

    # Cumulative vertical grammage from the surface downward [g/cm2]
    # rho [g/cm3] * dz [m] * 100 [cm/m] = dX [g/cm2]
    dX          = rho * dz_m * 100.0
    X_cumul     = np.zeros_like(depths)
    X_cumul[1:] = np.cumsum(dX[:-1])   # X_cumul[0] = 0 at surface

    with open(output_file, "w") as f:
        f.write(
            f"! pyMGMR3D ice profile  model={ice_model_id}: "
            f"{ICE_MODELS[ice_model_id]['name']}\n"
        )
        f.write("! depth[m]  density[g/cm3]  xi=n-1  X_cumul_vertical[g/cm2]\n")
        for z, r, x, X in zip(depths, rho, xi, X_cumul):
            f.write(f"{z:10.3f}  {r:12.6f}  {x:12.6f}  {X:14.4f}\n")

    if verbose:
        cos_zen = np.cos(np.radians(zenith_deg))
        print(f"Ice profile -> {output_file}")
        print(f"  Model:    {ICE_MODELS[ice_model_id]['name']}")
        print(f"  Steps:    {len(depths)}  (dz={dz_m} m, max={depth_max_m} m)")
        print(f"  Surface:  rho={rho[0]:.4f} g/cm3,  n={n[0]:.5f},  xi={xi[0]:.5f}")
        print(f"  Deep ice: rho={rho[-1]:.4f} g/cm3,  n={n[-1]:.5f},  xi={xi[-1]:.5f}")
        print(f"  Cherenkov angle deep ice: {np.degrees(np.arccos(1.0/n[-1])):.2f} deg")
        print(f"  Total vertical grammage:  {X_cumul[-1]:.1f} g/cm2")
        if zenith_deg > 0.0:
            print(f"  Slant grammage at zenith={zenith_deg:.1f} deg: "
                  f"{X_cumul[-1] / cos_zen:.1f} g/cm2")

    return {
        "depth": depths,
        "density": rho,
        "refractivity": xi,
        "X_cumul_vertical": X_cumul,
    }


# ---------------------------------------------------------------------------
# Top-level runner: reads .in file and generates current_ice.dat
# ---------------------------------------------------------------------------

def run_from_input_file(input_file, output_file="current_ice.dat", verbose=True):
    """
    Read ice parameters from a pyMGMR3D .in input file and generate
    current_ice.dat if ice_model_id > 0.

    This is the function called by pyMGMR3D.sh (or its Python equivalent)
    before launching the Fortran binary.
    """
    params = read_ice_params(input_file)

    if verbose:
        print(f"Ice parameters read from {input_file}:")
        for k, v in params.items():
            print(f"  {k} = {v}")

    if params["ice_model"] == 0:
        if verbose:
            print("ice_model=0: skipping ice profile generation.")
        return {}

    return generate_ice_profile(
        ice_model_id  = params["ice_model"],
        depth_max_m   = params["X_max_ice"],
        dz_m          = params["step_ice"],
        output_file   = output_file,
        zenith_deg    = params["Zen_sh"],
        verbose       = verbose,
    )
# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate current_ice.dat for pyMGMR3D in-ice simulations."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input", type=str, metavar="FILE",
        help="Read parameters from a pyMGMR3D .in input file (recommended)"
    )
    group.add_argument(
        "--model", type=int,
        help="Ice model ID directly: 1=constant, 2=S.Pole, 3=Greenland"
    )
    parser.add_argument("--depth",  type=float, default=500.0,
                        help="Max ice depth [m] (only used with --model, default 500)")
    parser.add_argument("--step",   type=float, default=5.0,
                        help="Depth step [m] (only used with --model, default 5)")
    parser.add_argument("--zenith", type=float, default=0.0,
                        help="Shower zenith [deg] (only used with --model, default 0)")
    parser.add_argument("--output", type=str,   default="current_ice.dat",
                        help="Output filename (default: current_ice.dat)")
    args = parser.parse_args()

    if args.input:
        run_from_input_file(args.input, output_file=args.output, verbose=True)
    else:
        generate_ice_profile(
            ice_model_id = args.model,
            depth_max_m  = args.depth,
            dz_m         = args.step,
            output_file  = args.output,
            zenith_deg   = args.zenith,
            verbose      = True,
        )
