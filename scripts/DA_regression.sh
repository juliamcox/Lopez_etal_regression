#!/bin/bash

#SBATCH --account=p31438
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=julia.cox@northwestern.edu
#SBATCH -o logs/%N.%j_%A_%a.out       # STDOUT
#SBATCH -e logs/%N.%j_%A_%a.err       # STDERR


module load matlab/r2022b

matlab -nodisplay -r "addpath(genpath('~/Documents/MATLAB/')); rng('shuffle'); DA_activeAvoid_script_noPar;  exit"