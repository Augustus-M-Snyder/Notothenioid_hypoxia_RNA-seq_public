#!/bin/bash
#SBATCH --job-name=FeatureCounts_samps_sep
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
species=("ACE" "GEO" "COR")

echo "Starting gene counts..."
# Script for generating gene count data with featureCounts:
for sp in "${species[@]}"; do
    echo "Processing species: $sp"

    # Loop through each BAM file in the directory
    for bam_file in bam/*${sp}_align.bam; do
        echo "Processing BAM file: $bam_file"

        # Run featureCounts for the current BAM file
        /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/featureCounts \
            -p --countReadPairs -O -T 12 -t exon -F GTF -g gene_id \
            -a Annotations/${sp}_genomic.gtf \
            -o Trimmed_reads/feature_counts/separated_samples/"${sp}_$(basename "$bam_file" .bam)_counts.txt" \
            "$bam_file"

        # Capture the exit status of featureCounts
        status=$?
        if [[ $status -ne 0 ]]; then
            echo "featureCounts failed for BAM file $bam_file with exit status $status."
            echo "Check the log file Trimmed_reads/feature_counts/${sp}_$(basename "$bam_file" .bam)_featureCounts.log for details."
        else
            echo "featureCounts completed successfully for BAM file $bam_file."
        fi
    done
done
echo "Done!"

