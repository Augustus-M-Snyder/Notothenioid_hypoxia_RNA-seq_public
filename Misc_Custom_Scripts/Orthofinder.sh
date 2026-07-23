#!/bin/bash
#SBATCH --job-name=Orthofinder_orthogroups
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --mem=128G
#SBATCH --mail-user=amsnyder7
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j_Orthofinder_orthogroups.txt

# Load necessary modules
module purge
module load GCCcore/12.3.0
module load slurm

# Load Conda Environment
conda init
source ~/.bashrc
conda activate NOTO_hyp

# Set Python path to NOTO_hyp environment
export PATH=~/miniconda3/envs/NOTO_hyp/bin:$PATH

# Run Orthofinder
~/miniconda3/envs/NOTO_hyp/bin/orthofinder -f ./prot_fastas
