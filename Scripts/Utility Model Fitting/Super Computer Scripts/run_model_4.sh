#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=04:00:00
#SBATCH --partition=amilan
#SBATCH --output=sample-matlab-model4-%j.out

module purge
module load matlab
matlab -nodesktop -nodisplay -r 'clear;utilitymodeltypes = 4;run_modelfits;'
