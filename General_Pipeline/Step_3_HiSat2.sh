#!/bin/bash
#SBATCH --job-name=Alignment_via_HiSAT2
#SBATCH --partition=t1standard
#SBATCH --nodes=3
#SBATCH --ntasks=36
#SBATCH --ntasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --mail-user=amsnyder7
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output=output_%j_Alignment_via_HISAT2.txt

# Load necessary modules
module purge
module load GCCcore/12.3.0
module load Python/3.11.3
module load Java
module load slurm

conda init
source ~/.bashrc
conda activate samtools_env

# Script for Aligning reads to genome with HiSat2:

echo "Aligning reads to genomes for each species..."
# Array of species prefixes
species=("ACE" "GEO" "COR" "MAC")

# Loop through each species
for sp in "${species[@]}"; do
    echo ${sp}
    # List all trimmed files for the current species
    ls ${sp}*.trimmed_1P.fq.gz | while read file; do
        name=$(echo "${file}" | cut -d "." -f 1)
        echo "${name}"
        # Run HiSat2 to align the trimmed files to the indexed genome file:
        /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/hisat2 -p 48 --mp 2,0 -q -x /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/indexes/${sp}_genome_tran \
	-1 "${name}.fastp.trimmed_1P.fq.gz" \
	-2 "${name}.fastp.trimmed_2P.fq.gz" \
	-U "${name}.fastp.trimmed_1U.fq.gz" \
	--summary-file "${name}_align.summ.txt" \
	--un "${name}.unmapped.fq" | \
	/home/amsnyder7/miniconda3/envs/samtools_env/bin/samtools view -S -b > /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/bam/${name}_${sp}_align.bam
	/home/amsnyder7/miniconda3/envs/samtools_env/bin/samtools sort /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/bam/${name}_${sp}_align.bam -o /center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/bam/${name}_${sp}_align.s.bam
    done
done

echo "Done!"
