#!/usr/bin/env bash
#SBATCH --partition=t1small
#SBATCH --job-name=map2zf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --export=ALL

# Submit with:
#   sbatch map_to_zebrafish_and_upset.sh -d /center1/NOTOHYP/amsnyder7/eggnog_data -s ACE,COR,GEO,GIB,MAC -p 8 -o emapper_zf
# Or as an array (one species per task):
#   sbatch --array=0-4 map_to_zebrafish_and_upset.sh -d /center1/.../eggnog_data -s ACE,COR,GEO,GIB,MAC -o emapper_zf

# -------------------- your original options & defaults --------------------
# Map longest CDS/proteins to zebrafish with eggNOG-mapper,
# standardize gene IDs, and build an overlap matrix.

# Defaults:
DATA_DIR="${EGGNOG_DATA_DIR:-$(pwd)/eggnog_data}"
SPECIES_CSV="ACE,COR,GEO,GIB,MAC"
CPUS=48
OUTDIR="emapper_zf"
PLOT="no"

usage() {
  echo "Usage: $0 [-d DATA_DIR] [-s SPECIES_CSV] [-p CPUS] [-o OUTDIR] [--plot]"
  echo "  -d  eggNOG database dir (must contain eggnog.db, eggnog_proteins.dmnd, eggnog.taxa.db)"
  echo "  -s  species list (comma-separated), default: ${SPECIES_CSV}"
  echo "  -p  threads for emapper (diamond), default: ${CPUS}"
  echo "  -o  output directory, default: ${OUTDIR}"
  echo "  --plot  also try to produce an UpSet plot (needs pandas, upsetplot, matplotlib)"
  exit 1
}

# -------------------- parse args --------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d) DATA_DIR="$2"; shift 2;;
    -s) SPECIES_CSV="$2"; shift 2;;
    -p) CPUS="$2"; shift 2;;
    -o) OUTDIR="$2"; shift 2;;
    --plot) PLOT="yes"; shift;;
    -h|--help) usage;;
    *) echo "[ERR] Unknown arg: $1"; usage;;
  esac
done

# If Slurm set cpus-per-task, use that (overrides -p).
if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
  CPUS="${SLURM_CPUS_PER_TASK}"
fi

# Load conda environment
conda init
source ~/.bashrc
conda activate NOTO_hyp

echo "[INFO] job_id=${SLURM_JOB_ID:-NA} node=${SLURM_JOB_NODELIST:-NA} cpus=${CPUS}"

mkdir -p "$OUTDIR"

# Check emapper
if ! command -v emapper.py >/dev/null 2>&1; then
  echo "[ERR] emapper.py not found in PATH (activate NOTO_hyp conda env?)"
  exit 1
fi

# -------------------- DB checks --------------------
pushd "$DATA_DIR" >/dev/null
for f in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db; do
  if [[ ! -s "$f" && -s "${f}.gz" ]]; then
    echo "[INFO] Decompressing ${f}.gz ..."
    gunzip -f "${f}.gz"
  fi
done
for must in eggnog.db eggnog_proteins.dmnd eggnog.taxa.db; do
  [[ -s "$must" ]] || { echo "[ERR] Missing $must in $DATA_DIR"; exit 1; }
done
popd >/dev/null

# Species list (and handle array mode)
IFS=',' read -r -a SPECIES <<< "$SPECIES_CSV"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx="${SLURM_ARRAY_TASK_ID}"
  if (( idx < 0 || idx >= ${#SPECIES[@]} )); then
    echo "[ERR] SLURM_ARRAY_TASK_ID=$idx is out of range for species list (${#SPECIES[@]} items)."
    exit 1
  fi
  SPECIES=( "${SPECIES[$idx]}" )
fi

echo "========== Mapping to zebrafish (taxid 7955) =========="
echo "[INFO] DATA_DIR=$DATA_DIR  CPUS=$CPUS  OUTDIR=$OUTDIR"
echo "[INFO] Species: ${SPECIES[*]}"

# -------------------- run emapper per species --------------------
for sp in "${SPECIES[@]}"; do
  in=""
  itype="proteins"
  for cand in "${sp}_cds_longest_per_gene.prot.fa" "${sp}_cds_longest_per_gene.faa" "${sp}_cds_longest_per_gene.fa"; do
    if [[ -s "$cand" ]]; then
      in="$cand"
      break
    fi
  done
  if [[ -z "$in" ]]; then
    echo "[WARN] No input FASTA for $sp; skipping (expected ${sp}_cds_longest_per_gene.{prot.fa,faa,fa})"
    continue
  fi
  case "$in" in
    *.prot.fa|*.faa) itype="proteins" ;;
    *.fa|*.fna|*.fas) itype="CDS" ;;
  esac

  echo "[RUN] $sp  (in=$in itype=$itype) -> ${OUTDIR}/${sp}_to_DR.emapper.annotations"
  emapper.py \
    -i "$in" \
    --itype "$itype" \
    -m diamond --sensmode sensitive \
    --cpu "$CPUS" \
    --data_dir "$DATA_DIR" \
    --target_taxa 7955 \
    --target_orthologs one2one \
    --report_orthologs \
    --no_file_comments \
    --override \
    --output_dir "$OUTDIR" \
    --output "${sp}_to_DR"
done

# -------------------- Parsing & Matrix (Python) --------------------
export OUTDIR SPECIES_CSV  # so Python can read them

python <<'PY'
import os, re, sys, gzip
from collections import defaultdict

outdir = os.environ.get("OUTDIR","emapper_zf")
species_csv = os.environ.get("SPECIES_CSV","ACE,COR,GEO,GIB,MAC")
species = [s for s in species_csv.split(",") if s]

os.makedirs(outdir, exist_ok=True)

def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

ensg_re = re.compile(r"(ENSDARG\d+)")
def parse_annotations(sp):
    ann = os.path.join(outdir, f"{sp}_to_DR.emapper.annotations")
    if not os.path.exists(ann) and os.path.exists(ann + ".gz"):
        ann = ann + ".gz"
    if not os.path.exists(ann):
        sys.stderr.write(f"[WARN] Missing annotations for {sp}: {ann}\n")
        return {}

    header = None
    rows = []
    with open_any(ann) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                h = line.lstrip("#").strip()
                if "\t" in h and "query_name" in h:
                    header = [c.strip() for c in h.split("\t")]
                continue
            rows.append(line.rstrip("\n"))
    if header is None:
        header = ["query_name","seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl","COG_category","Description","Preferred_name"]
    idx = {c:i for i,c in enumerate(header)}

    zf_map = {}
    mapped_any = 0
    for r in rows:
        parts = r.split("\t")
        if len(parts) < 2:
            continue
        q = parts[idx.get("query_name",0)]
        dr = None
        for p in parts:
            m = ensg_re.search(p)
            if m:
                dr = m.group(1)
                break
        if not dr:
            pn = parts[idx["Preferred_name"]] if "Preferred_name" in idx and idx["Preferred_name"] < len(parts) else ""
            pn = pn.strip()
            if pn and pn != "-":
                dr = pn
        if dr:
            zf_map[q] = dr
            mapped_any += 1
    sys.stderr.write(f"[INFO] {sp}: mapped {mapped_any} entries\n")

    pairs_path = os.path.join(outdir, f"{sp}_to_DR.ids.tsv")
    with open(pairs_path, "w") as fw:
        fw.write("zfish_key\tspecies_gene\n")
        for q,dr in zf_map.items():
            fw.write(f"{dr}\t{q}\n")

    keys_path = os.path.join(outdir, f"{sp}_zfish_keys.txt")
    with open(keys_path, "w") as fw:
        for dr in sorted(set(zf_map.values())):
            fw.write(f"{dr}\n")

    return zf_map

per_sp = {sp: parse_annotations(sp) for sp in species}

all_keys = set()
for m in per_sp.values():
    all_keys.update(m.values())

ordered_keys = sorted(all_keys)
matrix_path = os.path.join(outdir, "presence_matrix.tsv")
with open(matrix_path, "w") as fw:
    fw.write("zfish_key\t" + "\t".join(species) + "\n")
    for dr in ordered_keys:
        row = [dr]
        for sp in species:
            row.append("1" if dr in set(per_sp.get(sp,{}).values()) else "0")
        fw.write("\t".join(row) + "\n")

print(f"[OK] Wrote matrix: {matrix_path}")
PY

# -------------------- Optional UpSet plot --------------------
if [[ "$PLOT" == "yes" ]]; then
python <<'PY'
import os
outdir = os.environ.get("OUTDIR","emapper_zf")
path = os.path.join(outdir, "presence_matrix.tsv")

try:
    import matplotlib
    matplotlib.use("Agg")
    import pandas as pd
    import matplotlib.pyplot as plt
    from upsetplot import from_memberships, UpSet
except Exception:
    print("[WARN] Skipping plot (need pandas, matplotlib, upsetplot).")
    raise SystemExit(0)

df = pd.read_csv(path, sep="\t")
species_cols = [c for c in df.columns if c != "zfish_key"]

memberships = []
for _, r in df.iterrows():
    members = [sp for sp in species_cols if str(r[sp]) == "1"]
    if members:
        memberships.append(members)

if not memberships:
    print("[WARN] No intersections to plot.")
    raise SystemExit(0)

data = from_memberships(memberships)
plt.figure()
UpSet(data, subset_size='count', show_counts=True).plot()
plt.suptitle("Zebrafish ortholog overlap (eggNOG-mapper)")
png = os.path.join(outdir, "upset_all_species.png")
plt.savefig(png, dpi=200, bbox_inches="tight")
print(f"[OK] Wrote plot: {png}")
PY
fi

echo "[DONE] Outputs in: $OUTDIR"
