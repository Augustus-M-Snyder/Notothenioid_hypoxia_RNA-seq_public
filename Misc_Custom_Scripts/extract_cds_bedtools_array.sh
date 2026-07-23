#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Step 2: Extract CDS transcript FASTA with bedtools, then pick
#         the longest CDS per gene.
# Works in two modes:
#   - Single-species:   ./script -s GIB
#   - SLURM array:      sbatch --array=0-4 script
#                       (ACE=0, COR=1, GEO=2, MAC=3, GIB=4)
# Requires: bedtools, python3, samtools (for .fai if needed)
# Expects (from Step 1):
#   ./<SP>_cds.bed12
#   ../../Genomes/<SP>_genome.fa  (indexed optional: .fai)
#   optional ./<SP>_cds_index.tsv (tx -> gene map)
# Produces:
#   ./<SP>_cds_transcripts.fa
#   ./<SP>_cds_longest_per_gene.fa
#   ./<SP>_bedtools_getfasta.warn.log
# ============================================================

usage() {
  cat <<EOF
Usage:
  $0 -s <ACE|COR|GEO|MAC|GIB>     # run single species (local test)
  sbatch --array=0-4 $0           # SLURM array: 0=ACE,1=COR,2=GEO,3=MAC,4=GIB

Optional flags:
  -g <genome.fa>   override genome path (default ../../Genomes/<SP>_genome.fa)
  -b <cds.bed12>   override CDS BED12 (default ./<SP>_cds.bed12)
  -o <outdir>      outputs directory (default .)
  -h               help
EOF
}

# -------- arg parsing (single-species mode only) ----------
SP=""
GENOME=""
BED12=""
OUTDIR="."
while getopts "s:g:b:o:h" opt; do
  case "$opt" in
    s) SP="$OPTARG" ;;
    g) GENOME="$OPTARG" ;;
    b) BED12="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# -------- figure out species (single vs array) -------------
declare -a SPECIES_LIST=("ACE" "COR" "GEO" "MAC" "GIB")

if [[ -z "${SP}" ]]; then
  # No -s given → try SLURM array
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    idx="${SLURM_ARRAY_TASK_ID}"
    if (( idx < 0 || idx >= ${#SPECIES_LIST[@]} )); then
      echo "[ERROR] SLURM_ARRAY_TASK_ID=$idx out of range 0..$(( ${#SPECIES_LIST[@]}-1 ))" >&2
      exit 2
    fi
    SP="${SPECIES_LIST[$idx]}"
  else
    echo "[ERROR] No -s given and SLURM_ARRAY_TASK_ID is not set." >&2
    usage
    exit 2
  fi
fi

# -------- paths & checks -----------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME="${GENOME:-${SCRIPT_DIR}/../../Genomes/${SP}_genome.fa}"
BED12="${BED12:-${SCRIPT_DIR}/${SP}_cds.bed12}"
OUTDIR="${OUTDIR:-.}"

mkdir -p "$OUTDIR"
TX_FA="${OUTDIR}/${SP}_cds_transcripts.fa"
LP_FA="${OUTDIR}/${SP}_cds_longest_per_gene.fa"
BT_LOG="${OUTDIR}/${SP}_bedtools_getfasta.warn.log"
IDX_TSV="${SCRIPT_DIR}/${SP}_cds_index.tsv"   # optional fallback

echo "========== [${SP}] $(date) =========="
echo "[INFO] genome: $GENOME"
echo "[INFO] bed12 : $BED12"
echo "[INFO] outdir: $OUTDIR"

for f in "$GENOME" "$BED12"; do
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] Required file missing or empty: $f" >&2
    exit 3
  fi
done

# -------- Step 2.1: bedtools getfasta ----------------------
SORTED_BED="${OUTDIR}/.${SP}_cds.sorted.bed12.tmp"
if [[ -s "${GENOME}.fai" ]]; then
  bedtools sort -faidx "${GENOME}.fai" -i "$BED12" > "$SORTED_BED"
else
  bedtools sort -i "$BED12" > "$SORTED_BED"
fi

set +e
bedtools getfasta \
  -fi "$GENOME" \
  -bed "$SORTED_BED" \
  -s -split -name \
  -fo "$TX_FA" 2> "$BT_LOG"
rc=$?
set -e

if [[ $rc -ne 0 ]]; then
  echo "[ERROR] bedtools getfasta exit code $rc (see $BT_LOG)" >&2
  exit $rc
fi
if [[ ! -s "$TX_FA" ]]; then
  echo "[ERROR] ${TX_FA} is empty (see $BT_LOG)" >&2
  exit 4
fi

# -------- Step 2.2: longest CDS per gene -------------------
export SP OUTDIR
python3 - << 'PY'
import os, re, sys

sp       = os.environ.get("SP")
outdir   = os.environ.get("OUTDIR","." )
tx_fa    = os.path.join(outdir, f"{sp}_cds_transcripts.fa")
lp_fa    = os.path.join(outdir, f"{sp}_cds_longest_per_gene.fa")
idx_tsv  = os.path.join(os.path.dirname(__file__), f"{sp}_cds_index.tsv")

def itfa(path):
    h=None; seq=[]
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(seq)
                h=line.strip()[1:]
                seq=[]
            else:
                seq.append(line.strip())
    if h is not None:
        yield h, "".join(seq)

# optional transcript -> gene mapping
t2g = {}
if os.path.exists(idx_tsv) and os.path.getsize(idx_tsv) > 0:
    with open(idx_tsv) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 2:
                t2g[p[0]] = p[1]

def gene_from_header(h):
    m = re.search(r'\|gene:([^|]+)\|', h)
    if m:
        return m.group(1)
    t = h.split("|",1)[0]
    return t2g.get(t, t)

best = {}
n_tx = 0
for hdr, seq in itfa(tx_fa):
    n_tx += 1
    gid = gene_from_header(hdr)
    L   = len(seq)
    cur = best.get(gid)
    if (cur is None) or (L > cur[0]):
        best[gid] = (L, hdr, seq)

with open(lp_fa, "w") as fo:
    for L, hdr, seq in best.values():
        fo.write(f">{hdr}|len:{L}\n")
        for i in range(0, L, 60):
            fo.write(seq[i:i+60] + "\n")

print(f"[OK] {sp}: transcripts={n_tx} genes={len(best)} -> {lp_fa}")
PY

# -------- Summary ------------------------------------------
N_TX=$(grep -c '^>' "$TX_FA" || true)
N_GN=$(grep -c '^>' "$LP_FA" || true)
N_WARN=$(wc -l < "$BT_LOG" || echo 0)

echo "[INFO] Transcript CDS FASTA : ${TX_FA}  (entries: ${N_TX})"
echo "[INFO] Longest CDS per gene : ${LP_FA}  (entries: ${N_GN})"
echo "[INFO] bedtools warnings    : ${BT_LOG} (lines: ${N_WARN})"
echo "[DONE] [${SP}] $(date)"

