#!/usr/bin/env bash
# End-to-end FAD2 phylogenetics pipeline. Both tiers are always rooted on
# the Arabidopsis outgroup (see fetch_data.sh) -- there is no unrooted mode.
#
# Steps (per input set):
#   1. MAFFT    -- multiple sequence alignment (L-INS-i, --auto)
#   2. FastTree -- approximate maximum-likelihood tree (JTT+CAT, SH-like support)
#   3. reroot_tree.py -- reroot FastTree's tree on the Arabidopsis outgroup
#              (FastTree has no built-in outgroup option; its raw output
#               is an unrooted tree with an arbitrary root)
#   4. IQ-TREE -- ML tree with ModelFinder + 1000 ultrafast bootstrap
#              replicates, then also rerooted with reroot_tree.py (IQ-TREE's
#              own -o flag only affects which taxon is *drawn* at the root --
#              its .treefile output stays genuinely unrooted, confirmed by
#              its own report: "NOTE: Tree is UNROOTED although outgroup
#              taxon ... is drawn at root" -- so -o alone isn't enough)
#
#   IQ-TREE is run with a fixed random seed (-seed) AND forced to a single
#   thread (-nt 1, not -nt AUTO) for reproducibility: IQ-TREE's own docs
#   note that even with a fixed seed, a multi-threaded run isn't guaranteed
#   to produce identical output, since thread scheduling can affect the
#   order of operations in its stochastic tree search. -nt AUTO would also
#   pick a different thread count on a different machine's core count,
#   which is exactly the situation an independent reproduction run is in.
#   MAFFT and FastTree have no equivalent randomness to control here (no
#   bootstrap/stochastic search involved in the steps used).
#   5. pairwise_identity.py -- % identity matrix from the alignment
#   6. plot_tree.py -- render both rooted trees as PNGs (FastTree + IQ-TREE)
#
# Requires on PATH: mafft, FastTree, iqtree (the bioconda `iqtree` package's
# binary is named `iqtree` even for IQ-TREE 2.x/3.x -- there is no `iqtree2`
# binary), python3 (with matplotlib).
# Runtime: well under a minute for the 4-sequence core set; still only a
# few minutes for the 11-sequence legume expansion on a laptop. IQ-TREE's
# 1000-replicate ultrafast bootstrap is the slowest single step and still
# finishes in seconds to low minutes at this sequence count/length, even
# pinned to a single thread for reproducibility (see above).
#
# Usage:
#   bash scripts/run_pipeline.sh
#       Runs on data/fad2_seed.fasta -> results/ (whatever fetch_data.sh
#       last wrote there -- core+outgroup, or core+legumes+outgroup if
#       fetch_data.sh was last run with -l/--legumes).
#
#   bash scripts/run_pipeline.sh --both
#       Runs BOTH tiers in one invocation: the 4-sequence core set
#       (data/fad2_seed_core.fasta, core + outgroup) into results/, then
#       the 11-sequence legume expansion (data/fad2_seed_legumes.fasta,
#       core + legume expansion + outgroup) into results/legume_expansion/.
#       Requires having run `bash scripts/fetch_data.sh --legumes` at
#       least once first (that writes both seed files; see fetch_data.sh).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${ROOT}/data"

# Substring used to find the outgroup leaf in aligned/tree headers -- see
# fetch_data.sh for why Arabidopsis (UniProt P46313) is the outgroup.
OUTGROUP="P46313"

# Fixed for reproducibility -- see the header comment above. Change
# IQTREE_SEED if you deliberately want a different (still reproducible)
# run; the specific value doesn't matter, only that it's pinned.
IQTREE_SEED="42"
IQTREE_THREADS="1"

run_stage() {
  local DATA="$1" RESULTS="$2" LABEL="$3"
  mkdir -p "$RESULTS"

  if [[ ! -s "$DATA" ]]; then
    echo "ERROR: ${DATA} not found or empty. Run scripts/fetch_data.sh first." >&2
    exit 1
  fi

  local N_SEQ
  N_SEQ=$(grep -c '^>' "$DATA")
  local NSTEPS=6
  echo "=== ${LABEL} (${N_SEQ} sequences) -> ${RESULTS} ==="

  echo "[1/${NSTEPS}] Aligning ${N_SEQ} sequences with MAFFT..."
  mafft --auto --thread 2 "$DATA" > "${RESULTS}/fad2_aligned.fasta" 2> "${RESULTS}/mafft.log"

  echo "[2/${NSTEPS}] Building tree with FastTree..."
  FastTree "${RESULTS}/fad2_aligned.fasta" > "${RESULTS}/fad2_tree.nwk" 2> "${RESULTS}/fasttree.log"

  echo "[3/${NSTEPS}] Rerooting FastTree's tree on the outgroup..."
  python3 "${ROOT}/scripts/reroot_tree.py" "${RESULTS}/fad2_tree.nwk" \
    "${RESULTS}/fad2_tree_rooted.nwk" "$OUTGROUP"

  local OUTGROUP_LABEL
  OUTGROUP_LABEL=$(grep "^>" "${RESULTS}/fad2_aligned.fasta" | grep "$OUTGROUP" | head -1 | sed 's/^>//; s/ .*//')
  local IQTREE_OUTGROUP_ARGS=()
  if [[ -z "$OUTGROUP_LABEL" ]]; then
    echo "  WARNING: could not find outgroup '${OUTGROUP}' in aligned headers;" >&2
    echo "  IQ-TREE's report will not have an outgroup taxon to note." >&2
  else
    IQTREE_OUTGROUP_ARGS=(-o "$OUTGROUP_LABEL")
  fi

  local HAVE_IQTREE=0
  if [[ "$N_SEQ" -ge 4 ]]; then
    echo "[4/${NSTEPS}] Building bootstrap-supported tree with IQ-TREE (ModelFinder + 1000 UFBoot)..."
    rm -f "${RESULTS}"/fad2_iqtree.*
    iqtree -s "${RESULTS}/fad2_aligned.fasta" -st AA -m MFP -bb 1000 \
      -nt "$IQTREE_THREADS" -seed "$IQTREE_SEED" \
      "${IQTREE_OUTGROUP_ARGS[@]}" \
      -pre "${RESULTS}/fad2_iqtree" > "${RESULTS}/iqtree.log" 2>&1
    cp "${RESULTS}/fad2_iqtree.treefile" "${RESULTS}/fad2_tree_iqtree.nwk"
    # -o above only affects which taxon IQ-TREE's own report draws at the
    # root; the .treefile itself is still unrooted (see header comment).
    # Reroot it for real with the same tool used for FastTree's tree.
    python3 "${ROOT}/scripts/reroot_tree.py" "${RESULTS}/fad2_iqtree.treefile" \
      "${RESULTS}/fad2_tree_iqtree.nwk" "$OUTGROUP"
    HAVE_IQTREE=1
  else
    echo "[4/${NSTEPS}] Skipping IQ-TREE bootstrap: only ${N_SEQ} sequences (bootstrap needs >= 4)."
  fi

  echo "[5/${NSTEPS}] Computing pairwise percent identity..."
  python3 "${ROOT}/scripts/pairwise_identity.py" "${RESULTS}/fad2_aligned.fasta" \
    > "${RESULTS}/pairwise_identity.tsv"

  echo "[6/${NSTEPS}] Rendering tree figure(s)..."
  python3 "${ROOT}/scripts/plot_tree.py" "${RESULTS}/fad2_tree_rooted.nwk" "${RESULTS}/fad2_tree_rooted.png" \
    "FAD2 protein tree (FastTree, JTT+CAT, rooted on Arabidopsis outgroup)" "SH-like local support"
  if [[ "$HAVE_IQTREE" -eq 1 ]]; then
    python3 "${ROOT}/scripts/plot_tree.py" "${RESULTS}/fad2_tree_iqtree.nwk" "${RESULTS}/fad2_tree_iqtree.png" \
      "FAD2 protein tree (IQ-TREE, ModelFinder + UFBoot), rooted on Arabidopsis outgroup" "ultrafast bootstrap % (1000 reps)"
  fi

  echo
  echo "Done. Outputs in ${RESULTS}:"
  ls -1 "${RESULTS}"
  echo
}

if [[ "${1:-}" == "--both" ]]; then
  if [[ ! -s "${DATA_DIR}/fad2_seed_legumes.fasta" ]]; then
    echo "ERROR: ${DATA_DIR}/fad2_seed_legumes.fasta not found." >&2
    echo "Run 'bash scripts/fetch_data.sh --legumes' first -- it writes both" >&2
    echo "the core-only and legume-expansion seed files that --both needs." >&2
    exit 1
  fi
  run_stage "${DATA_DIR}/fad2_seed_core.fasta"    "${ROOT}/results"                    "Core set (rooted)"
  run_stage "${DATA_DIR}/fad2_seed_legumes.fasta" "${ROOT}/results/legume_expansion"   "Legume expansion (rooted)"
  echo "Both tiers complete: results/ (core set) and results/legume_expansion/ (legume expansion)."
else
  run_stage "${DATA_DIR}/fad2_seed.fasta" "${ROOT}/results" "Default run (rooted)"
fi
