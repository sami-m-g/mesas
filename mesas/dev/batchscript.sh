#!/bin/bash
#SBATCH --job-name=mesas
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=end
#SBATCH --mail-user=charman1@jhu.edu

#### load and unload modules you may need
module load intel
module load python/3.7-anaconda
. /software/apps/anaconda/5.2/python/3.7/etc/profile.d/conda.sh
conda activate
conda activate py3
mkdir ~/scratch/$SLURM_JOBID
mkdir ~/scratch/$SLURM_JOBID/junk
mkdir ~/scratch/$SLURM_JOBID/junk/plots
mkdir ~/scratch/$SLURM_JOBID/dev
mkdir ~/scratch/$SLURM_JOBID/data
cp ./lowerhafren.py ~/scratch/$SLURM_JOBID/dev
cp ./plots.py ~/scratch/$SLURM_JOBID/dev
cp ./recursive_split.py ~/scratch/$SLURM_JOBID/dev
cp ../data/lower_hafren.csv ~/scratch/$SLURM_JOBID/data
cd ~/scratch/$SLURM_JOBID/dev
python ~/scratch/$SLURM_JOBID/dev/lowerhafren.py
mv ~/scratch/$SLURM_JOBID ~/data/charman1/mesas


