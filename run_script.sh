#!/bin/bash

# Runscript for the BlueCloud pipeline
# Author: Urs Hofmann Elizondo, 28/02/2024

##################################################################
# MAKE SURE TO SWITCH TO THE NEW SOFTWARE STACK
##################################################################
source /cluster/apps/local/env2lmod.sh

# Remove modules
module purge

# Load all needed modules
module load gcc/8.2.0
module load gdal/3.5.3
module load proj/8.2.1
module load geos/3.9.1
module load udunits2/2.2.28
module load r/4.2.2

# --vanilla means R starts with a clean slate, without saving or restoring anything from previous sessions
# --no-echo supresses the printing commands

sbatch -n 8 --time=3:00:00 --mem-per-cpu=15G --wrap "R --vanilla < master.R > ./data/EULER_cache/beta_div.out" --output=./data/EULER_cache/slurm_%j.out

