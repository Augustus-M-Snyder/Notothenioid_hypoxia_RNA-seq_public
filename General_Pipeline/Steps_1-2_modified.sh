#!/bin/bash
#SBATCH --job-name=Fastqc_&_Trimming_RNAseq
#SBATCH --partition=t1standard
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
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
source activate NOTO_hyp # This didn't seem to work before with source because of the word activate...

# Script for BBduk:
# Loop over files starting from GIB628
ls *_R1_001.fastq.gz | while read file; do
    name=$(echo "${file}" | cut -d "/" -f 8 | cut -d "." -f 1 | cut -d "_" -f 1)
    if [ "$name" \< "GIB631" ]; then
        echo "Skipping $name"
        continue
    fi
    echo "${name}"                                                                  
    if [ -f "fastq_trimmed/${name}.1.bb.fq.gz" ]; then
        echo "Output files exist for $name; skipping..."
    else
        /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/bbduk.sh -Xmx1g \
            in1="${name}_R1_001.fastq.gz" \
            in2="${name}_R2_001.fastq.gz" \
            out1="fastq_trimmed/${name}.1.bb.fq.gz" \
            out2="fastq_trimmed/${name}.2.bb.fq.gz" \
            ref=/home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/adapters.fa \
            ktrim=r k=23 mink=11 hdist=1 tpe tbo
    fi
done
echo "Done!"

# Script for Trimmomatic (this didn't get done with previous job)
ls *.1.bb.fq.gz | while read file; do 
        name=$(echo "${file}" | cut -d "." -f 1) 		
        echo "${name}" 
        if [ -f "${name}.fastp.trimmed_1P.fq.gz" ]; then 		
                echo "exists; skipping..."				
        else

                /home/amsnyder7/miniconda3/envs/NOTO_hyp/bin/trimmomatic PE -threads 8 -phred33 -trimlog \
                	"${name}.trimmomatic.log" "${file}" "${name}.2.bb.fq.gz" \
                	-baseout "${name}.fastp.trimmed.fq.gz" \
                	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:36 
        fi
done
echo "Done!"

# Deactivate conda environment (optional)
source deactivate

# End of script
