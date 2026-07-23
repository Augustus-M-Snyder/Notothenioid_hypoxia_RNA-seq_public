#!/bin/bash
#SBATCH --job-name=bbduk_trim
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --mail-user=amsnyder7@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j.txt


# Load necessary modules
module purge
module load GCCcore/12.3.0
module load Python/3.11.3
module load Java
module load slurm

# Activate conda environment
conda init
source ~/.bashrc
conda activate NOTO_hyp

# Script for BBduk:
ls *_R1_001.fastq.gz | while read file; do
    name=$(echo "${file}" | cut -d "/" -f 8 | cut -d "." -f 1 | cut -d "_" -f 1)
    echo "${name}"                                                                  
    if [ -f "./${name}.1.bb.fq.gz" ]; then
        echo "Output files exist for $name; skipping..."
    else
        /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/bbduk.sh -Xmx1g \
            in1="${name}_R1_001.fastq.gz" \
            in2="${name}_R2_001.fastq.gz" \
            out1="./${name}.1.bb.fq.gz" \
            out2="./${name}.2.bb.fq.gz" \
            ref=/home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/adapters.fa \
            ktrim=r k=23 mink=11 hdist=1 tpe tbo
    fi
done
echo "Done!"

# End of script
