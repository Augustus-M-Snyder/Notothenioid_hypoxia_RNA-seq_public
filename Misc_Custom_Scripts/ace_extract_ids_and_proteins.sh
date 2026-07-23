#!/usr/bin/env bash
set -euo pipefail

ACE_GTF="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Annotations/ACE_annotation.gtf"
OUT_DIR="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Annotations/Symbol_Lists"
mkdir -p "$OUT_DIR"

read_in() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then gzip -cd "$f"; else cat "$f"; fi
}

# --- 1A. Unique gene_ids from 'gene' features ---
read_in "$ACE_GTF" | awk -F'\t' '
function get(attr, key,   m){
  if (match(attr, key" \"([^\"]*)\"", m)) return m[1];    # GTF key "value"
  if (match(attr, key"=([^;]+)", m))     return m[1];    # GFF key=value
  return ""
}
$0 !~ /^#/ && $3=="gene" {
  gid = get($9, "gene_id")
  if (gid != "" && gid != "." && gid != "-") print gid
}' | sort -u > "${OUT_DIR}/ACE_gene_ids.txt"

# --- 1B. gene_id <tab> protein_id from CDS features (many-to-one allowed) ---
read_in "$ACE_GTF" | awk -F'\t' '
function get(attr, key,   m){
  if (match(attr, key" \"([^\"]*)\"", m)) return m[1];
  if (match(attr, key"=([^;]+)", m))     return m[1];
  return ""
}
$0 !~ /^#/ && $3=="CDS" {
  gid = get($9, "gene_id")
  pid = get($9, "protein_id")
  if (gid != "" && pid != "" && pid != "-" && pid != ".") {
    print gid "\t" pid
  }
}' | sort -u > "${OUT_DIR}/ACE_gene_to_protein.tsv"

wc -l "${OUT_DIR}/ACE_gene_ids.txt"
wc -l "${OUT_DIR}/ACE_gene_to_protein.tsv"
echo "Wrote:"
echo "  ${OUT_DIR}/ACE_gene_ids.txt"
echo "  ${OUT_DIR}/ACE_gene_to_protein.tsv"

