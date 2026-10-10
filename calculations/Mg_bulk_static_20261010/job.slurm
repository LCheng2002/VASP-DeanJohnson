#!/usr/bin/env bash
#SBATCH --job-name=MgStaticRef
#SBATCH --partition=64c512g
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:40:00
#SBATCH --mem=16G
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail
module purge
module load vasp/6.3.2-vtst-intel-2021.4.0 || module load vasp/6.3.2-intel-2021.4.0
command -v vasp_std
ulimit -s unlimited
mpirun -np "$SLURM_NTASKS" vasp_std
