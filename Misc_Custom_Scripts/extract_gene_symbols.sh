#!/usr/bin/env bash
set -euo pipefail

# ======= USER SETTINGS (edit as needed) =======
# Species tags you use throughout the project
species=("ACE" "GEO" "MAC" "COR" "GIB")

# Where your annotation files live. One file per species, e.g. ACE_annotation.gtf(.gz)
ANN_DIR="/center1/NOTOHYP/amsnyder7/RNA_seq/Noto_hyp_fastq/fastq_trimmed/Annotations"
# Output folder for symbol lists
OUT_DIR="${ANN_DIR}/Symbol_Lists"
mkdir -p "$OUT_DIR"

# If your files don’t follow the pattern below, set explicit paths in the case block.
# Expected filenames:
#   ACE_annotation.gtf(.gz), GEO_annotation.gtf(.gz), ... (adjust if different)
# ==============================================

# Function to cat or zcat transparently
read_in() {
  local f="$1"
  if [[ "$f" == *.gz ]]; then gzip -cd "$f"; else cat "$f"; fi
}

# AWK program:
# - Looks for common keys in priority order: gene_name, gene_symbol, Name, gene
# - Accepts either GTF-like (key "value") or GFF-like (key=value)
# - Prints non-empty values; we later sort -u
awk_program='
BEGIN { FS = "\t" }
function get_attr(attr, keys,   n,i,k,pat,m) {
  n = split(keys, K, ",")
  for (i=1; i<=n; i++) {
    k = K[i]
    # GTF style: key "value"
    pat = k " \"([^\"]+)\""
    if (match(attr, pat, m)) return m[1]
    # GFF style: key=value
    pat = k "=([^;]+)"
    if (match(attr, pat, m)) return m[1]
  }
  return ""
}
{
  # Column 9 is attributes
  sym = get_attr($9, "gene_name,gene_symbol,Name,gene")
  if (sym != "" && sym != "-" && sym != ".") {
    print sym
  }
}
'

for sp in "${species[@]}"; do
  # Try gz then plain; adjust pattern here if your filenames differ
  cand1="${ANN_DIR}/${sp}_annotation.gtf.gz"
  cand2="${ANN_DIR}/${sp}_annotation.gtf"
  if [[ -f "$cand1" ]]; then
    ann="$cand1"
  elif [[ -f "$cand2" ]]; then
    ann="$cand2"
  else
    echo "WARNING: No annotation file found for ${sp} (looked for ${cand1} and ${cand2})" >&2
    continue
  fi

  echo "[${sp}] Extracting gene symbols from: $ann"
  read_in "$ann" \
    | awk "$awk_program" \
    | sed 's/^"//; s/"$//' \
    | sort -u \
    > "${OUT_DIR}/${sp}_gene_symbols.txt"

  # Quick report
  n=$(wc -l < "${OUT_DIR}/${sp}_gene_symbols.txt")
  echo "[${sp}] Wrote ${n} unique symbols -> ${OUT_DIR}/${sp}_gene_symbols.txt"
done

echo "Done."

