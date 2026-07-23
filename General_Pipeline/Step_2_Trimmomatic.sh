#!/bin/bash
#SBATCH --job-name=Fastqc_&_Trimming_RNAseq
#SBATCH --partition=t1standard
#SBATCH --nodes=3
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
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

# Function to calculate percentage
calculate_percentage() {
    local current=$1
    local total=$2
    echo $(awk "BEGIN { pc=100*${current}/${total}; i=int(pc); print (pc-i<0.5)?i:i+1 }")
}

# Function to display the progress bar
show_progress() {
    local progress=$1
    local total=$2
    local percent=$(( 100 * progress / total ))
    local filled=$(( percent / 2 ))
    local empty=$(( 50 - filled ))
    local bar=$(printf "%0.s#" $(seq 1 $filled))
    bar+=$(printf "%0.s." $(seq 1 $empty))
    printf "\r[%-50s] %d%%" "$bar" "$percent"
}

# Count total number of files
total_files=$(ls *.1.bb.fq.gz | wc -l)
processed_files=0

# Initialize progress bar
show_progress $processed_files $total_files

# Loop through input files
for file in *.1.bb.fq.gz; do
    name=$(echo "${file}" | cut -d "." -f 1)
    echo "Processing file: $file"

    # Check if trimmed output file already exists
    if [ -f "${name}.fastp.trimmed_1P.fq.gz" ]; then
        echo "Trimmed file ${name}.fastp.trimmed_1P.fq.gz exists; skipping..."
    else
        # Run Trimmomatic
        /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/trimmomatic PE -threads 8 -phred33 -trimlog "${name}.trimmomatic.log" \
            "${file}" "${name}.2.bb.fq.gz" \
            -baseout "${name}.fastp.trimmed.fq.gz" \
            ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:36

        # Check if Trimmomatic produced the expected output files
        if [ -f "${name}.fastp.trimmed_1P.fq.gz" ]; then
            echo "Trimming completed successfully for $file"

            # Clean up intermediate files if necessary
            echo "Deleting bbduk files for disk space..."
            rm -f "${file}" "${name}.2.bb.fq.gz"
        else
            echo "Error: Trimming failed for $file"
        fi
    fi

    # Update progress
    processed_files=$((processed_files + 1))
    show_progress $processed_files $total_files
done

echo -e "\nTrimming and processing complete!"
echo "Done!"

# Deactivate conda environment (optional)
conda deactivate

# End of script

