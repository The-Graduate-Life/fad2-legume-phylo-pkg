#!/usr/bin/env bash
# End-to-end FAD2 phylogenetics pipeline.
#
# Steps:
#   1. MAFFT   -- multiple sequence alignment (L-INS-i, --auto)
#   2. FastTree -- approximate maximum-likelihood tree (JTT+CAT, SH-like support)
#   3. pairwise_identity.py -- % identity matrix from the alignment
#   4. plot_tree.py -- render the Newick tree as a PNG
#
# Requires on PATH: mafft, FastTree, python3 (with matplotlib).
# Runtime: well under a minute for the 3-sequence example set; still only
# a few minutes even with a few dozen sequences on a laptop.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="${ROOT}/data/fad2_seed.fasta"
RESULTS="${ROOT}/results"
mkdir -p "$RESULTS"

if [[ ! -s "$DATA" ]]; then
  echo "ERROR: ${DATA} not found or empty. Run scripts/fetch_data.sh first." >&2
  exit 1
fi

echo "[1/4] Aligning $(grep -c '^>' "$DATA") sequences with MAFFT..."
mafft --auto --thread 2 "$DATA" > "${RESULTS}/fad2_aligned.fasta" 2> "${RESULTS}/mafft.log"

echo "[2/4] Building tree with FastTree..."
FastTree "${RESULTS}/fad2_aligned.fasta" > "${RESULTS}/fad2_tree.nwk" 2> "${RESULTS}/fasttree.log"

echo "[3/4] Computing pairwise percent identity..."
python3 "${ROOT}/scripts/pairwise_identity.py" "${RESULTS}/fad2_aligned.fasta" \
  > "${RESULTS}/pairwise_identity.tsv"

echo "[4/4] Rendering tree figure..."
python3 "${ROOT}/scripts/plot_tree.py" "${RESULTS}/fad2_tree.nwk" "${RESULTS}/fad2_tree.png"

echo
echo "Done. Outputs in ${RESULTS}:"
ls -1 "${RESULTS}"
