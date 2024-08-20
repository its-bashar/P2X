#!/bin/bash
#SBATCH --job-name=h5ad_files_p2x_filtering
#SBATCH --cpus-per-task=50
#SBATCH --mem=239G
#SBATCH --time=01:00:00
#SBATCH --output=h5ad_files_p2x_filtering_%j.out
#SBATCH --error=h5ad_files_p2x_filtering_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=USERNAME@student.le.ac.uk
#SBATCH --export=NONE

# Load the Python module
module load python/3.10.12-tggsi7t

# Create a virtual environment in a directory where you have write permissions
VENV_DIR="$HOME/venv_$SLURM_JOB_ID"
python -m venv $VENV_DIR

# Activate the virtual environment
source $VENV_DIR/bin/activate

# Upgrade pip and install scanpy within the virtual environment
pip install --upgrade pip
pip install scanpy
pip install numpy
pip install pandas

# Run the Python script
python 4.h5ad_files_p2x_filtering_converting_to_csvs.py

# Deactivate the virtual environment
deactivate

# Clean up the virtual environment
rm -rf $VENV_DIR

