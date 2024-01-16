#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
#SBATCH --partition=amilan
#SBATCH --output=sample-matlab-model2-%j.out

module purge
module load matlab
matlab -nodesktop -nodisplay -r 'clear;utilitymodeltypes = 2;run_modelfits;'
