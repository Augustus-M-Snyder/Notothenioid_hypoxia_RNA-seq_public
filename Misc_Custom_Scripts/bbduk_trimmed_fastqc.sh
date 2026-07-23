#!/bin/bash
#SBATCH --job-name=Fastqc_bbduk_Emac
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=amsnyder7@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j_%x.txt


# Load necessary modules
module purge
module load Python/3.11.3
module load Java
module load slurm

# Activate conda environment
conda init
source ~/.bashrc
source activate NOTO_hyp

# Script for Fastqc
mkdir -p bbduk_trimmed_fastqc_output

echo "Evaluating read count files for each species..."
# Array of species prefixes
species=("MAC")

# Loop through each species
for sp in "${species[@]}"; do
    echo ${sp}
        for file in ${sp}*bb.fq.gz; do
        fastqc -o ./bbduk_trimmed_fastqc_output "$file";
        done
done
echo "Done!"

# Deactivate conda environment (optional)
conda deactivate

# End of script
