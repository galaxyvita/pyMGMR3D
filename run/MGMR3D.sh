#!/bin/bash

# --- Validate input ---
# Checks if at least one argument (the input file) is provided.
# If not, it prints a usage message and exits with an error.
if [ $# -lt 1 ]; then
    echo "Usage: $0 <inputfile.in> [--plot] [--no-hdf5] [--plot-dir <directory>]"
    exit 1
fi

# Assign the first argument as the input file.
INPUT_FILE=$1
# Shift arguments so that $1 now refers to the next argument (flags).
shift

# Extract the base name of the input file (without extension).
# Get the current working directory, which is considered the run folder.
# Determine the program directory by going up one level from the run folder.
BASENAME=$(basename "$INPUT_FILE" .in)
RUN_FOLDER=$(pwd)
PROG_DIR=$(dirname "$RUN_FOLDER")/program

# Initialize flags for plotting and HDF5 conversion as enabled by default.
# Initialize plot directory.
# Initialize plot directory.
# Flag to track if a custom plot directory was specified.
# Initialize HDF5 filename, it will be the same as the input file with .hdf5 extension.
PLOT_ENABLED=false
HDF5_ENABLED=true
PLOT_DIR=""
CUSTOM_PLOT_DIR=false
HDF5_FILENAME="${RUN_FOLDER}/${BASENAME}.hdf5"

# --- Parse flags ---
# Loop through all provided command-line arguments to parse flags.
while [[ $# -gt 0 ]]; do
    case $1 in
        --plot)
            PLOT_ENABLED=true
            shift
            ;;
        --hdf5)
            HDF5_ENABLED=false
            shift
            ;;
        --plot-dir)
            PLOT_DIR=$2
            CUSTOM_PLOT_DIR=true
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Set the default plot directory to the run folder if no custom directory was specified.
if [ -z "$PLOT_DIR" ]; then
    PLOT_DIR="${RUN_FOLDER}"
fi

# Create the plot directory if it doesn't already exist.
mkdir -p "$PLOT_DIR"

# Print important paths for user information.
echo "Program folder = ${PROG_DIR}"
echo "Run folder     = ${RUN_FOLDER}"
echo "Input file     = ${INPUT_FILE}"
echo "Plot output to = ${PLOT_DIR}"
echo "HDF5 output to = ${HDF5_FILENAME}"

# --- Compile MGMR ---
# Change directory to the program directory.
# If changing directory fails.
# Compile the MGMR program using the specified makefile.
# If compilation fails, exit the script.
# Change back to the original run folder.
# If changing directory fails, exit the script.
cd "${PROG_DIR}" || exit 1
make -f "MGMR3D_fit-makefile-v5.mak" || exit 1
cd "${RUN_FOLDER}" || exit 1

# --- Run simulation ---
# Create a temporary 'plot' directory within the run folder for intermediate output.
# Execute the MGMR program, feeding it the input file's content.
mkdir plot
"${PROG_DIR}/MGMR3D_fit-v5" < "$INPUT_FILE"

# --- Plotting with GLE and Python ---
# Change directory to the temporary 'plot' folder to handle output files.
cd plot || exit 1

# Check if plotting is enabled.
if $PLOT_ENABLED; then
    # Check if GLE (Graphics Layout Engine) command is available.
    if command -v gle &> /dev/null; then
        gle -d pdf -o "${PLOT_DIR}/${BASENAME}_FitStokes.pdf" "${PROG_DIR}/FitStokes.GLE" "${RUN_FOLDER}/plot/FitResult"
        gle -d jpg -r 200 -o "${PLOT_DIR}/${BASENAME}_FitStokes-map.jpg" "${PROG_DIR}/FitStokes-map.GLE" "${RUN_FOLDER}/plot/"
        gle -d pdf -o "${PLOT_DIR}/${BASENAME}_sh-current.pdf" "${PROG_DIR}/sh-current.GLE" "${RUN_FOLDER}/plot/"
    else
        echo "GLE not found. Skipping GLE plots."
    fi
    # Check if python3 command is available and the required Python libraries.
    if command -v python3 &> /dev/null; then
        python3 - <<EOF
try:
    import pandas
    import matplotlib
    import numpy
except ImportError as e:
    print(f"Missing library: {e.name}")
    exit(1)
EOF
        if [ $? -eq 0 ]; then
            python3 "${PROG_DIR}/FitStokes.py" "${RUN_FOLDER}/plot/FitResult.dat" "${PLOT_DIR}"
            python3 "${PROG_DIR}/FitStokes-map.py" "${RUN_FOLDER}/plot/" "${PLOT_DIR}"
            python3 "${PROG_DIR}/sh-current.py" "${RUN_FOLDER}/plot/" "${PLOT_DIR}"
        else
            echo "Skipping Python plotting due to missing libraries."
        fi
    else
        echo "Python3 not installed. Skipping plot generation."
    fi
fi

# --- Convert to HDF5 ---
# Check if HDF5 conversion is enabled and the required Python libraries.
if $HDF5_ENABLED; then
    if command -v python3 &> /dev/null; then
        python3 - <<EOF
try:
    import pandas as pd
    import numpy as np
    import h5py
    from pathlib import Path
    import re
except ImportError as e:
    print(f"Missing library: {e.name}")
    exit(1)
EOF
        if [ $? -eq 0 ]; then
            python3 "${RUN_FOLDER}/MGMR_HDF5.py" "${RUN_FOLDER}/${BASENAME}.hdf5" "${RUN_FOLDER}/${INPUT_FILE}" "${RUN_FOLDER}/plot/FitResult.dat" "${RUN_FOLDER}/plot/sh_Current.dat" "${RUN_FOLDER}/plot/"
        else
            echo "Skipping HDF5 conversion due to missing libraries."
        fi
    else
        echo "Python3 not installed. Skipping HDF5 generation."
    fi
fi

# --- Cleanup only if no custom output ---
# Change back to the original run folder.
cd "${RUN_FOLDER}" || exit 1
# Remove the temporary 'plot' directory ONLY if a custom plot directory was NOT specified.
# This prevents deleting desired output if the user explicitly set a plot-dir.
if ! $CUSTOM_PLOT_DIR; then
    rm -r plot
fi

# Inform the user that the script has completed.
# Exit with a success status.
echo "Run complete."
exit 0
