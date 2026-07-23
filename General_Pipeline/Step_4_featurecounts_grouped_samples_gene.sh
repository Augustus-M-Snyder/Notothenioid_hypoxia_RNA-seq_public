#!/bin/bash
#SBATCH --job-name=FeatureCounts
#SBATCH --partition=t1small
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --mail-user=amsnyder7@alaska.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j_%x.txt

# Load necessary modules
module purge
module load GCCcore/12.3.0
module load Python/3.11.3
module load Java
module load slurm

conda init
source ~/.bashrc
conda activate NOTO_hyp

# Change to the directory containing the BAM files and annotation files
cd /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed

# Array of species prefixes
species=("ACE" "GEO" "COR" "MAC")

echo "Starting gene counts..."
# Script for generating gene count data with featureCounts:
for sp in "${species[@]}"; do
    echo "Processing species: $sp"

    # Run featureCounts with all .bam files for the current species using wildcard and enable verbose reporting
    /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/featureCounts \
        -p --countReadPairs -O -T 12 -t gene -F GTF -g gene_id \
        -a Annotations/${sp}_annotation.gtf \
        -o Trimmed_reads/feature_counts_test/genes/${sp}_counts.txt \
        bam/*${sp}_align.bam

    # Capture the exit status of featureCounts
    status=$?
    if [[ $status -ne 0 ]]; then
        echo "featureCounts failed for species $sp with exit status $status."
        echo "Check the log file Trimmed_reads/feature_counts/${sp}_featureCounts.log for details."
    else
        echo "featureCounts completed successfully for species $sp."
    fi
done
echo "Done!"


