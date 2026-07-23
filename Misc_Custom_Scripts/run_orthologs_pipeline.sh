#!/usr/bin/env bash

# Load conda environment
conda init
source ~/.bashrc
conda activate NOTO_hyp
set -euo pipefail

# ===== User paths (adjust if needed) =====
BASE_DIR="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
PMAP_DIR="$BASE_DIR/Annotations/protein_mappings"
GENOME_DIR="$BASE_DIR/Genomes"
ANN_DIR="$BASE_DIR/Annotations"

cd "$PMAP_DIR"

# ===== Species =====
SPECIES=("ACE" "COR" "GEO" "MAC" "GIB")

# ===== Tools (module load if needed) =====
# module load DIAMOND/2.x
# module load gffread/0.12.x
# module load python/3.10

# ===== STEP A. Ensure per-species protein FASTAs exist; build gene<->protein maps when missing =====
# For each species:
#   - If <spp>_protein_sequences.fa exists and non-empty, keep it.
#   - Else, generate CDS/protein from GTF+genome using gffread.
#   - Build/refresh <spp>_protein_mapping.tsv: gene_id, transcript_id, protein_id, gene_name_raw, gene_name_std
python3 - <<'PY'
import os, re, sys, subprocess, pandas as pd

BASE_DIR = "/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed"
PMAP_DIR = f"{BASE_DIR}/Annotations/protein_mappings"
GENOME_DIR = f"{BASE_DIR}/Genomes"
ANN_DIR = f"{BASE_DIR}/Annotations"
species = ["ACE","COR","GEO","MAC","GIB"]

os.chdir(PMAP_DIR)

def std_symbol(s):
    if s is None or pd.isna(s) or str(s).strip()=="":
        return None
    # normalize whitespace and case; keep a reproducible version for joining
    t = str(s).strip()
    return t.upper()  # store uppercase as a standardized join key; keep original separately

def ensure_proteins(sp):
    prot_fa = f"{sp}_protein_sequences.fa"
    if os.path.exists(prot_fa) and os.path.getsize(prot_fa) > 0:
        return prot_fa
    # Build from genome + GTF
    ann = f"{ANN_DIR}/{sp}_annotation.gtf"
    genome = f"{GENOME_DIR}/{sp}_genome.fa"
    cds_fa = f"{sp}_cds.fa"

    if not os.path.exists(ann):
        print(f"[WARN] Missing annotation for {sp}: {ann}", file=sys.stderr)
    if not os.path.exists(genome):
        print(f"[WARN] Missing genome for {sp}: {genome}", file=sys.stderr)

    # gffread to get CDS and proteins
    # -x CDS, -y proteins
    cmd = ["gffread", "-g", genome, "-x", cds_fa, "-y", prot_fa, ann]
    subprocess.run(cmd, check=True)
    return prot_fa

def parse_gtf_for_map(gtf_path):
    # Return dataframe with gene_id, transcript_id, gene_name_raw
    # We’ll read only lines with 'transcript' or 'mRNA' features to bind transcript->gene and gene_name
    rows = []
    attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
    with open(gtf_path, 'r') as fh:
        for line in fh:
            if not line or line.startswith("#"): 
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9: 
                continue
            feature = parts[2]
            if feature not in ("transcript","mRNA"): 
                continue
            attrs = parts[8]
            d = {}
            for m in attr_re.finditer(attrs):
                key, val = m.group(1), m.group(2)
                d[key] = val
            gid = d.get("gene_id")
            tid = d.get("transcript_id")
            gname = d.get("gene_name") or d.get("gene_symbol") or d.get("Name")
            rows.append((gid, tid, gname))
    if not rows:
        return pd.DataFrame(columns=["gene_id","transcript_id","gene_name_raw"])
    df = pd.DataFrame(rows, columns=["gene_id","transcript_id","gene_name_raw"]).drop_duplicates()
    # Deduplicate preferring non-null names
    df = df.sort_values(["gene_id","transcript_id","gene_name_raw"], na_position="last").drop_duplicates(["gene_id","transcript_id"], keep="first")
    return df

def parse_fasta_headers_for_protein_ids(fa_path):
    # gffread writes headers like:
    # >transcript_id pep|someID|... or just >transcript_id
    # We’ll try to capture a protein-like token if present; otherwise use transcript_id as protein_id surrogate.
    recs = []
    with open(fa_path, 'r') as fh:
        tid = None
        prot_id = None
        for line in fh:
            if line.startswith(">"):
                hdr = line[1:].strip()
                # transcript id is first token until space
                tid = hdr.split()[0]
                prot_id = None
                # look for something that looks like a protein/accession-esque token in the rest of header
                # fallback to tid
                prot_id = tid
                recs.append((tid, prot_id))
    df = pd.DataFrame(recs, columns=["transcript_id","protein_id"]).drop_duplicates()
    return df

for sp in species:
    gtf = f"{ANN_DIR}/{sp}_annotation.gtf"
    prot_fa = ensure_proteins(sp)
    # Build map table
    df_t = parse_gtf_for_map(gtf) if os.path.exists(gtf) else pd.DataFrame(columns=["gene_id","transcript_id","gene_name_raw"])
    df_p = parse_fasta_headers_for_protein_ids(prot_fa)
    df = df_p.merge(df_t, how="left", on="transcript_id")
    df["gene_name_std"] = df["gene_name_raw"].apply(std_symbol)
    out_tsv = f"{sp}_protein_mapping.tsv"
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[OK] {sp}: proteins={prot_fa}, mapping={out_tsv}  (n={len(df)})")
PY

# ===== STEP B. Download zebrafish reference proteome + GTF and make DIAMOND DB =====
# Using Ensembl current (adjust release if you prefer pinning)
DR_DIR="drerio_ref"
mkdir -p "$DR_DIR"
cd "$DR_DIR"

if [ ! -s "Danio_rerio.GRCz11.pep.all.fa.gz" ]; then
  # Try Ensembl/current or the latest mirror (adjust if needed)
  wget -q --show-progress https://ftp.ensembl.org/pub/current/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz || \
  wget -q --show-progress https://ftp.ensembl.org/pub/release-111/fasta/danio_rerio/pep/Danio_rerio.GRCz11.pep.all.fa.gz
fi

if [ ! -s "Danio_rerio.GRCz11.112.gtf.gz" ]; then
  wget -q --show-progress https://ftp.ensembl.org/pub/current/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz || \
  wget -q --show-progress https://ftp.ensembl.org/pub/release-112/gtf/danio_rerio/Danio_rerio.GRCz11.112.gtf.gz
fi

if [ ! -s "drerio.pep.fa" ]; then
  pigz -dc Danio_rerio.GRCz11.pep.all.fa.gz > drerio.pep.fa
fi

if [ ! -s "drerio.gtf" ]; then
  pigz -dc Danio_rerio.GRCz11.112.gtf.gz > drerio.gtf
fi

# Build DIAMOND DB (once)
if [ ! -s "drerio.dmnd" ]; then
  diamond makedb --in drerio.pep.fa -d drerio
fi

# Build a protein->(gene_id,gene_symbol) map from GTF + FASTA headers
python3 - <<'PY'
import re, pandas as pd

pep = "drerio.pep.fa"
gtf = "drerio.gtf"

# Map ENSEMBL transcript->gene_id, gene_name from GTF
rows=[]
attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
with open(gtf,'r') as fh:
    for line in fh:
        if line.startswith("#"): continue
        parts=line.rstrip("\n").split("\t")
        if len(parts)<9: continue
        if parts[2] not in ("transcript","mRNA"): continue
        d={}
        for m in attr_re.finditer(parts[8]):
            d[m.group(1)] = m.group(2)
        tid = d.get("transcript_id")
        gid = d.get("gene_id")
        gnm = d.get("gene_name") or d.get("Name")
        if tid:
            rows.append((tid,gid,gnm))
df_t = pd.DataFrame(rows, columns=["transcript_id","dr_gene_id","dr_gene_symbol"]).drop_duplicates()

# FASTA headers have ENSEMBL protein IDs and linked transcript IDs
# Example header: >ENSDARP00000123456 pep:known chromosome:... transcript:ENSDART00000... gene:ENSDARG0000... gene_symbol:actb1
records=[]
with open(pep,'r') as fh:
    pid=None; tid=None; gid=None; gsym=None
    for line in fh:
        if line.startswith(">"):
            hdr=line[1:].strip()
            pid=hdr.split()[0]
            tid=None; gid=None; gsym=None
            # try to parse transcript: and gene_symbol: tokens if present
            for tok in hdr.split():
                if tok.startswith("transcript:"):
                    tid=tok.split(":",1)[1]
                elif tok.startswith("gene:"):
                    gid=tok.split(":",1)[1]
                elif tok.startswith("gene_symbol:"):
                    gsym=tok.split(":",1)[1]
            records.append((pid, tid, gid, gsym))
df_p = pd.DataFrame(records, columns=["dr_protein_id","transcript_id","dr_gene_id_hdr","dr_gene_symbol_hdr"]).drop_duplicates()

# Join preference: use gene_symbol from header if present else from GTF transcript map
df = df_p.merge(df_t, how="left", on="transcript_id")
df["dr_gene_id"] = df["dr_gene_id_hdr"].fillna(df["dr_gene_id"])
df["dr_gene_symbol"] = df["dr_gene_symbol_hdr"].fillna(df["dr_gene_symbol"])
df = df[["dr_protein_id","dr_gene_id","dr_gene_symbol"]].drop_duplicates()
df.to_csv("drerio_protein_to_gene.tsv", sep="\t", index=False)
print(f"[OK] drerio_protein_to_gene.tsv    n={len(df)}")
PY

cd "$PMAP_DIR"

# ===== STEP C. DIAMOND BLASTP forward (species->drerio) and reverse (drerio->species) =====
mkdir -p diamond_out

DR_DB="$PMAP_DIR/drerio_ref/drerio.dmnd"

for sp in "${SPECIES[@]}"; do
  SP_PEP="${sp}_protein_sequences.fa"
  if [ ! -s "$SP_PEP" ]; then
    echo "[WARN] Missing $SP_PEP, skipping $sp"
    continue
  fi

  # Forward
  FWD="diamond_out/${sp}_vs_dr.tsv"
  if [ ! -s "$FWD" ]; then
    diamond blastp -q "$SP_PEP" -d "$DR_DB" -o "$FWD" \
      --max-target-seqs 1 --evalue 1e-5 --very-sensitive --outfmt 6 qseqid sseqid pident length evalue bitscore
  fi

  # Reverse DB creation (species)
  if [ ! -s "diamond_out/${sp}.dmnd" ]; then
    diamond makedb --in "$SP_PEP" -d "diamond_out/${sp}"
  fi

  # Reverse
  REV="diamond_out/dr_vs_${sp}.tsv"
  if [ ! -s "$REV" ]; then
    diamond blastp -q "drerio_ref/drerio.pep.fa" -d "diamond_out/${sp}.dmnd" -o "$REV" \
      --max-target-seqs 1 --evalue 1e-5 --very-sensitive --outfmt 6 qseqid sseqid pident length evalue bitscore
  fi
done

# ===== STEP D. Compute RBH 1:1 orthologs per species and build clean mapping tables =====
python3 - <<'PY'
import os, pandas as pd, glob

PMAP_DIR = "."
DRMAP = "drerio_ref/drerio_protein_to_gene.tsv"
dr = pd.read_csv(DRMAP, sep="\t")

def load_map(sp):
    m = pd.read_csv(f"{sp}_protein_mapping.tsv", sep="\t")
    # ensure columns
    for c in ["gene_id","transcript_id","protein_id","gene_name_raw","gene_name_std"]:
        if c not in m.columns: m[c] = None
    return m

def top1(df, qcol, scol, keep="first"):
    # df: qseqid sseqid ... choose single top hit per q based on best bitscore (and then lowest evalue)
    # Input from DIAMOND is already --max-target-seqs 1, so it’s unique. But keep safety anyway.
    return df.drop_duplicates([qcol], keep=keep)

all_out = []
for fwd in glob.glob("diamond_out/*_vs_dr.tsv"):
    sp = os.path.basename(fwd).split("_vs_")[0]
    rev = f"diamond_out/dr_vs_{sp}.tsv"
    if not os.path.exists(rev): 
        continue

    df_f = pd.read_csv(fwd, sep="\t", header=None,
                       names=["qseqid","sseqid","pident","length","evalue","bitscore"])
    df_r = pd.read_csv(rev, sep="\t", header=None,
                       names=["qseqid","sseqid","pident","length","evalue","bitscore"])

    df_f = top1(df_f, "qseqid","sseqid")
    df_r = top1(df_r, "qseqid","sseqid")

    # RBH condition: species_protein <-> dr_protein
    # forward: species_q (protein_id) -> dr_s (dr_protein_id)
    # reverse: dr_q (dr_protein_id)   -> sp_s (species protein_id)
    rbh = df_f.merge(df_r, left_on="sseqid", right_on="qseqid", suffixes=("_fwd","_rev"))
    rbh = rbh[ rbh["qseqid_rev"] == rbh["qseqid_fwd"] ]  # dr_protein back to same species protein

    # Map to gene IDs / symbols (species & zebrafish)
    spmap = load_map(sp)
    # qseqid_fwd is species protein_id; link to transcript -> gene
    rbh = rbh.merge(spmap[["transcript_id","protein_id","gene_id","gene_name_raw","gene_name_std"]], 
                    left_on="qseqid_fwd", right_on="protein_id", how="left")
    rbh = rbh.merge(dr, left_on="sseqid_fwd", right_on="dr_protein_id", how="left")

    # Keep a clean table
    keep = ["gene_id","transcript_id","protein_id",
            "gene_name_raw","gene_name_std",
            "dr_protein_id","dr_gene_id","dr_gene_symbol",
            "pident_fwd","length_fwd","evalue_fwd","bitscore_fwd"]
    rbh = rbh.rename(columns={"pident_fwd":"pident","length_fwd":"aln_len",
                              "evalue_fwd":"evalue","bitscore_fwd":"bitscore"})
    rbh = rbh[keep].drop_duplicates()

    # One species gene -> one dr_gene_id: pick best bitscore per species gene_id
    rbh = rbh.sort_values(["gene_id","bitscore","-pident".replace("-","")], ascending=[True,False,False]).drop_duplicates(["gene_id"], keep="first")

    rbh.insert(0, "Species", sp)
    out = f"RBH_{sp}_to_DR.tsv"
    rbh.to_csv(out, sep="\t", index=False)
    print(f"[OK] {out}   n={len(rbh)}")
    all_out.append(rbh)

if all_out:
    allcat = pd.concat(all_out, ignore_index=True)
    allcat.to_csv("RBH_all_species_to_DR.tsv", sep="\t", index=False)
    # Wide matrix: rows = dr_gene_id, columns species -> 1 if mapped
    wide = allcat.assign(val=1).pivot_table(index="dr_gene_id", columns="Species", values="val", fill_value=0, aggfunc="max")
    wide.to_csv("RBH_presence_matrix.tsv", sep="\t")
PY

echo "[DONE] All outputs in: $PMAP_DIR"

