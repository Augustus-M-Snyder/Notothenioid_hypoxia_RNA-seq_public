#!/usr/bin/env bash
#SBATCH -J GIBacc2gene
#SBATCH -p standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o GIBacc2gene.%j.out
#SBATCH -e GIBacc2gene.%j.err

set -euo pipefail

# -------------------- CONFIG (edit if needed) --------------------
BASE="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
ANN_DIR="${BASE}/Annotations"
GEN_DIR="${BASE}/Genomes"
PROT_DIR="${ANN_DIR}/protein_mappings"
OUT_DIR="${ANN_DIR}/Symbol_Lists"

GIB_GTF="${ANN_DIR}/GIB_annotation.gtf"
GIB_GENOME="${GEN_DIR}/GIB_genome.fa"

ACC_TSV="${PROT_DIR}/GIB_accessions.tsv"   # OR provide TXT:
ACC_TXT="${PROT_DIR}/GIB_accessions.txt"

MIN_PIDENT=60     # % identity threshold for best hit
MIN_QCOV=0.70     # query coverage (0-1)

CONDA_EXE="$HOME/miniconda3/bin/conda"
ENV_GFFREAD="NOTO_hyp"
ENV_EDIRECT="edirect"
ENV_EMAPPER="emapper"

mkdir -p "$PROT_DIR" "$OUT_DIR"

echo "== GIB accession → gene mapper =="
echo "[CONFIG] GTF:     $GIB_GTF"
echo "[CONFIG] GENOME:  $GIB_GENOME"
echo "[CONFIG] OUT_DIR: $OUT_DIR"

# -------------------- STEP 0: pre-clean GTF ----------------------
CLEAN_GTF="${OUT_DIR}/GIB_annotation.noTx.gtf"

echo "[clean] Stripping only transcript/mRNA feature rows (duplicates cause gffread crashes)..."
awk 'BEGIN{OFS="\t"} 
     /^#/ {next} 
     { if ($3!="transcript" && $3!="mRNA") print $0 }' \
     "$GIB_GTF" > "$CLEAN_GTF"

echo "[clean] Cleaned GTF → $CLEAN_GTF"

# -------------------- STEP 1: gffread proteins/CDS ----------------
if [[ ! -s "$OUT_DIR/GIB_proteins.fa" || ! -s "$OUT_DIR/GIB_cds.fa" ]]; then
  echo "[gffread] Attempt 1 on original GTF (env: $ENV_GFFREAD)..."
  if ! "$CONDA_EXE" run -n "$ENV_GFFREAD" gffread -g "$GIB_GENOME" \
        -y "$OUT_DIR/GIB_proteins.fa" \
        -x "$OUT_DIR/GIB_cds.fa" \
        "$GIB_GTF" ; then
    echo "[gffread] Attempt 1 failed; trying cleaned GTF..."
    "$CONDA_EXE" run -n "$ENV_GFFREAD" gffread -g "$GIB_GENOME" \
        -y "$OUT_DIR/GIB_proteins.fa" \
        -x "$OUT_DIR/GIB_cds.fa" \
        "$CLEAN_GTF"
  fi
else
  echo "[gffread] Found existing proteins/CDS; skipping."
fi

# -------------------- STEP 2: make transcript→gene & gene→symbol ---
echo "[maps] Extracting transcript→gene and gene→symbol tables..."
python - << 'PY'
import re, csv, os
OUT=os.environ["OUT_DIR"]; GTF=os.environ["GIB_GTF"]
CLEAN=os.environ["CLEAN_GTF"]
path = GTF if os.path.getsize(GTF)>0 else CLEAN
# Prefer cleaned file for mapping too (avoids dup IDs in parent chains)
path = CLEAN if os.path.getsize(CLEAN)>0 else GTF

def attrs(s):
    d={}
    for k,v in re.findall(r'(\S+)\s+"([^"]+)"', s): d[k]=v
    for k,v in re.findall(r'(\S+)=([^;]+)', s): d[k]=v
    return d

gene_name_by_gene={}; tx_to_gene={}
with open(path, "r", errors="ignore") as fh:
    for ln in fh:
        if not ln or ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t"); 
        if len(f)<9: continue
        ft, a = f[2], attrs(f[8])
        gid = a.get("gene_id") or a.get("gene") or a.get("ID")
        tid = a.get("transcript_id") or (a.get("ID") if ft in ("transcript","mRNA") else None)
        par = a.get("Parent")
        if a.get("gene_name") and gid: gene_name_by_gene.setdefault(gid, a["gene_name"])
        if ft in ("transcript","mRNA") and tid and (gid or par):
            tx_to_gene[tid] = gid or par

os.makedirs(OUT, exist_ok=True)
with open(f"{OUT}/GIB_transcript_to_gene.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["transcript_id","gene_id"])
    for k,v in sorted(tx_to_gene.items()): w.writerow([k,v])
with open(f"{OUT}/GIB_gene_id_to_symbol.tsv","w",newline="") as f:
    w=csv.writer(f,delimiter="\t"); w.writerow(["gene_id","gene_name"])
    for k,v in sorted(gene_name_by_gene.items()): w.writerow([k,v])
print("[maps] Wrote:", f"{OUT}/GIB_transcript_to_gene.tsv", f"{OUT}/GIB_gene_id_to_symbol.tsv")
PY

# -------------------- STEP 3: accession FASTA via EDirect ----------
# Normalize accession list to one-per-line
if [[ -s "$ACC_TSV" ]]; then
  echo "[accessions] Using TSV → column 2"
  cut -f2 "$ACC_TSV" | grep -v -E '^\s*$' | sort -u > "$OUT_DIR/GIB_accessions.txt"
elif [[ -s "$ACC_TXT" ]]; then
  echo "[accessions] Using TXT list"
  grep -v -E '^\s*$' "$ACC_TXT" | sort -u > "$OUT_DIR/GIB_accessions.txt"
else
  echo "[ERROR] No accession list at $ACC_TSV or $ACC_TXT" >&2
  exit 2
fi

echo "[efetch] Downloading protein sequences (env: $ENV_EDIRECT)..."
> "$OUT_DIR/GIB_accessions.faa"
> "$OUT_DIR/GIB_accessions_unfetched.txt"
BATCH=200
mapfile -t IDS < <(cat "$OUT_DIR/GIB_accessions.txt")
for ((i=0; i<${#IDS[@]}; i+=BATCH)); do
  ids_str=$(IFS=, ; echo "${IDS[@]:i:BATCH}")
  if ! "$CONDA_EXE" run -n "$ENV_EDIRECT" efetch -db protein -id "$ids_str" -format fasta >> "$OUT_DIR/GIB_accessions.faa" 2>/dev/null; then
    printf "%s\n" "${IDS[@]:i:BATCH}" >> "$OUT_DIR/GIB_accessions_unfetched.txt"
  fi
done
ACC_N=$(grep -c '^>' "$OUT_DIR/GIB_accessions.faa" || true)
echo "[efetch] Wrote $ACC_N sequences → $OUT_DIR/GIB_accessions.faa"

# -------------------- STEP 4: DIAMOND align accessions → proteins ---
echo "[diamond] Build DB from GIB_proteins.fa (env: $ENV_EMAPPER)..."
"$CONDA_EXE" run -n "$ENV_EMAPPER" diamond makedb \
  --in "$OUT_DIR/GIB_proteins.fa" -d "$OUT_DIR/GIB_proteins.dmnd"

echo "[diamond] BLASTP accessions vs GIB proteins..."
"$CONDA_EXE" run -n "$ENV_EMAPPER" diamond blastp \
  -q "$OUT_DIR/GIB_accessions.faa" \
  -d "$OUT_DIR/GIB_proteins.dmnd" \
  -o "$OUT_DIR/GIB_accession_vs_GIB_proteins.tsv" \
  --max-target-seqs 5 --evalue 1e-6 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore qlen slen

# -------------------- STEP 5: best hits + join to gene/symbol -------
if [[ ! -s "$OUT_DIR/GIB_accession_vs_GIB_proteins.tsv" ]]; then
  echo "[diamond] No hits file created; aborting." >&2
  exit 0
fi

awk -v P="$MIN_PIDENT" -v Q="$MIN_QCOV" 'BEGIN{OFS="\t"}
{
  qcov = ($4>0 && $7>0) ? $4/$7 : 0;
  scov = ($4>0 && $8>0) ? $4/$8 : 0;
  if ($3+0 >= P && qcov >= Q) print $1,$2,$3,qcov,scov,$5,$6
}' "$OUT_DIR/GIB_accession_vs_GIB_proteins.tsv" \
| sort -k1,1 -k7,7nr -k3,3nr \
| awk -v OFS="\t" '!seen[$1]++' \
> "$OUT_DIR/GIB_besthits.filtered.tsv"

# Load maps and emit final table
awk -v OFS="\t" -v T2G="$OUT_DIR/GIB_transcript_to_gene.tsv" -v G2S="$OUT_DIR/GIB_gene_id_to_symbol.tsv" '
  BEGIN{
    # transcript -> gene
    while ((getline < T2G) > 0) { if (++n1==1) continue; split($0,a,"\t"); t2g[a[1]]=a[2] }
    close(T2G)
    # gene -> symbol
    while ((getline < G2S) > 0) { if (++n2==1) continue; split($0,b,"\t"); g2s[b[1]]=b[2] }
    close(G2S)
    print "protein_accession","transcript_id","gene_id","gene_name","pident","qcov","scov","evalue","bitscore"
  }
  {
    acc=$1; tid=$2; pident=$3; qcov=$4; scov=$5; evalue=$6; bits=$7
    # gffread uses transcript IDs as sequence IDs; keep the token before any spaces/pipes
    split(tid,parts,/[\|\t ]/); tid=parts[1]
    gid = (tid in t2g) ? t2g[tid] : ""
    gnm = (gid in g2s) ? g2s[gid] : ""
    print acc, tid, gid, gnm, pident, qcov, scov, evalue, bits
  }
' "$OUT_DIR/GIB_besthits.filtered.tsv" > "$OUT_DIR/GIB_accession_to_gene_map.tsv"

echo "[DONE] → $OUT_DIR/GIB_accession_to_gene_map.tsv"

