#!/usr/bin/env bash
set -euo pipefail

# --- config you can tweak ---
ENV_NAME="NOTO_hyp"
: "${EGGNOG_DATA_DIR:=/center1/NOTOHYP/amsnyder7/eggnog_data}"
mkdir -p "$EGGNOG_DATA_DIR"

# Mirrors that still host emapperdb-5.0.2
MIRRORS=(
  "https://eggnog6.embl.de/download/emapperdb-5.0.2"
  "https://eggnog5.embl.de/download/emapperdb-5.0.2"
)

# Files we need (minimal set for DIAMOND mode)
FILES=(
  "eggnog.db.gz"
  "eggnog_proteins.dmnd.gz"
  # taxonomic helper DB (small; nice to have)
  "eggnog.taxa.tar.gz"
)

# --- conda activate ---
if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda not on PATH"; exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "[INFO] Using env: $CONDA_DEFAULT_ENV"
echo "[INFO] Data dir: $EGGNOG_DATA_DIR"
echo "[INFO] emapper:  $(emapper.py --version || true)"
echo "[INFO] diamond:  $(diamond version 2>/dev/null | head -n1 || true)"
echo

fetch_one() {
  local fname="$1"
  local ok=1
  for base in "${MIRRORS[@]}"; do
    url="${base}/${fname}"
    echo "[DL] ${url}"
    if wget -O "${EGGNOG_DATA_DIR}/${fname}" -q --show-progress --progress=dot:giga "$url"; then
      ok=0; break
    else
      echo "[WARN] failed: ${url}"
    fi
  done
  return $ok
}

# download each file (skip if already non-empty)
for f in "${FILES[@]}"; do
  if [ -s "${EGGNOG_DATA_DIR}/${f}" ]; then
    echo "[SKIP] ${f} already present (size: $(stat -c%s "${EGGNOG_DATA_DIR}/${f}") bytes)"
    continue
  fi
  fetch_one "$f" || { echo "[ERROR] Could not download $f from any mirror"; exit 2; }
done

# decompress / unpack safely
cd "$EGGNOG_DATA_DIR"

if [ -s eggnog.db.gz ]; then
  echo "[INFO] Decompressing eggnog.db.gz ..."
  gunzip -f eggnog.db.gz
fi

if [ -s eggnog_proteins.dmnd.gz ]; then
  echo "[INFO] Decompressing eggnog_proteins.dmnd.gz ..."
  gunzip -f eggnog_proteins.dmnd.gz
fi

if [ -s eggnog.taxa.tar.gz ]; then
  echo "[INFO] Unpacking eggnog.taxa.tar.gz ..."
  tar -xzf eggnog.taxa.tar.gz
  # keep the tarball to avoid re-downloads
fi

# sanity checks
fail=0
for must in eggnog.db eggnog_proteins.dmnd; do
  if [ ! -s "$must" ]; then
    echo "[ERROR] Missing or empty: $EGGNOG_DATA_DIR/$must"
    fail=1
  fi
done
[ $fail -eq 0 ] || exit 3

echo
echo "[OK] eggNOG DB ready in: $EGGNOG_DATA_DIR"
ls -lh "$EGGNOG_DATA_DIR" | sed -n '1,40p'
echo
echo "[TIP] When running emapper, pass:  --data_dir \"$EGGNOG_DATA_DIR\""

