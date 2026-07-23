#!/usr/bin/env bash
set -euo pipefail

# ---- config you may tweak ----
CPUS="${SLURM_CPUS_PER_TASK:-8}"
EMAPPER="${EMAPPER:-emapper.py}"     # in your eggnog-mapper conda env
EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR:-$HOME/eggnog_data}"  # point to your DB
OUTDIR="${OUTDIR:-emapper_out}"
SENS="${SENS:-diamond}"              # diamond|mmseqs
# --------------------------------

mkdir -p "$OUTDIR"

# Map array index -> species, or accept -s SPECIES
declare -a SPECIES=(ACE COR GEO GIB MAC)
SPEC=""; while getopts "s:" o; do case $o in s) SPEC="$OPTARG";; esac; done

if [[ -z "${SPEC}" ]]; then
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    SPEC="${SPECIES[$SLURM_ARRAY_TASK_ID]}"
  else
    echo "Usage: $0 -s {ACE|COR|GEO|GIB|MAC}   (or: sbatch --array=0-4 $0)"; exit 1
  fi
fi

IN="${SPEC}_cds_longest_per_gene.prot.fa"
[[ -s "$IN" ]] || { echo "[ERR] missing $IN"; exit 2; }

PREFIX="${OUTDIR}/${SPEC}"
LOG="${PREFIX}.log"

echo "[INFO] $(date) SPEC=$SPEC CPUS=$CPUS IN=$IN" | tee "$LOG"

# If needed, uncomment and adapt to your environment:
# module load miniconda; source activate eggnog-mapper
# or: mamba activate eggnog-mapper

"$EMAPPER" \
  -i "$IN" --itype proteins \
  --output "$PREFIX" --output_dir "$OUTDIR" \
  --cpu "$CPUS" --data_dir "$EGGNOG_DATA_DIR" \
  --override --temp_dir "${TMPDIR:-/tmp}" \
  --decorate_gff no --sensmode "$SENS" \
  2>&1 | tee -a "$LOG"

echo "[DONE] $(date) -> ${PREFIX}.emapper.annotations"

