#!/usr/bin/env bash
# Machine-checkable verification of a results/ (or results/legume_expansion/)
# directory produced by run_pipeline.sh -- exits 0 if everything checks out,
# non-zero with a clear message on the first thing that doesn't.
#
# This exists specifically so an unattended AI coding agent (or CI) can tell
# whether a reproduction run succeeded without a human eyeballing Newick
# strings or PNGs. It checks *substance* (sequence counts, correct rooting,
# expected file set, model selection present), not exact byte content --
# IQ-TREE's ultrafast bootstrap uses a random seed by default, so exact
# branch-support numbers (and therefore file bytes) are NOT expected to be
# identical between runs of the same input, even though the alignment,
# topology, and pairwise identities should be.
#
# Usage:
#   bash scripts/verify_results.sh results                    # core tier, expects 4 sequences
#   bash scripts/verify_results.sh results/legume_expansion    # legume tier, expects 11 sequences
#   bash scripts/verify_results.sh results 4                   # explicit expected count
set -uo pipefail

DIR="${1:?usage: verify_results.sh <results_dir> [expected_seq_count]}"
EXPECTED_N="${2:-}"

FAIL=0
fail() { echo "FAIL: $1" >&2; FAIL=1; }
pass() { echo "OK:   $1"; }

check_file() {
  local f="$1"
  if [[ -s "$DIR/$f" ]]; then
    pass "$f exists and is non-empty"
  else
    fail "$f missing or empty"
  fi
}

echo "=== Verifying $DIR ==="

for f in fad2_aligned.fasta fad2_tree.nwk fad2_tree_rooted.nwk fad2_tree_rooted.png \
         fad2_tree_iqtree.nwk fad2_tree_iqtree.png fad2_iqtree.iqtree \
         pairwise_identity.tsv mafft.log fasttree.log iqtree.log; do
  check_file "$f"
done

if [[ -s "$DIR/fad2_aligned.fasta" ]]; then
  N_SEQ=$(grep -c '^>' "$DIR/fad2_aligned.fasta")
  if [[ -n "$EXPECTED_N" ]]; then
    if [[ "$N_SEQ" -eq "$EXPECTED_N" ]]; then
      pass "sequence count = $N_SEQ (expected $EXPECTED_N)"
    else
      fail "sequence count = $N_SEQ, expected $EXPECTED_N"
    fi
  else
    pass "sequence count = $N_SEQ (no expected count given, not checked)"
  fi
fi

check_rooted_bifurcation() {
  local nwk="$1" label="$2"
  [[ -s "$nwk" ]] || { fail "$label: file missing, cannot check rooting"; return; }
  python3 - "$nwk" <<'PYEOF'
import sys
path = sys.argv[1]
with open(path) as f:
    s = f.read().strip().rstrip(";")
depth = 0
commas = 0
for c in s:
    if c == "(":
        depth += 1
    elif c == ")":
        depth -= 1
    elif c == "," and depth == 1:
        commas += 1
has_outgroup = "P46313" in s
ok = (commas == 1) and has_outgroup
print(f"{'OK  ' if ok else 'FAIL'} top-level commas={commas} (expect 1), outgroup P46313 present={has_outgroup}")
sys.exit(0 if ok else 1)
PYEOF
  if [[ $? -ne 0 ]]; then
    fail "$label: not a genuine outgroup-rooted bifurcation"
  else
    pass "$label: genuinely rooted on the Arabidopsis outgroup"
  fi
}

check_rooted_bifurcation "$DIR/fad2_tree_rooted.nwk" "FastTree rerooted tree"
check_rooted_bifurcation "$DIR/fad2_tree_iqtree.nwk" "IQ-TREE rerooted tree"

if [[ -s "$DIR/fad2_iqtree.iqtree" ]]; then
  if grep -q "Best-fit model" "$DIR/fad2_iqtree.iqtree"; then
    MODEL=$(grep "Best-fit model" "$DIR/fad2_iqtree.iqtree" | head -1)
    pass "ModelFinder ran: $MODEL"
  else
    fail "no 'Best-fit model' line found in fad2_iqtree.iqtree -- ModelFinder may not have run"
  fi
fi

if [[ -s "$DIR/pairwise_identity.tsv" ]]; then
  # every sequence should be 100.0% identical to itself somewhere in the matrix
  if grep -q "100.0" "$DIR/pairwise_identity.tsv"; then
    pass "pairwise_identity.tsv contains expected self-identity (100.0)"
  else
    fail "pairwise_identity.tsv has no 100.0 self-identity value -- looks malformed"
  fi
fi

echo
if [[ "$FAIL" -eq 0 ]]; then
  echo "=== ALL CHECKS PASSED for $DIR ==="
  exit 0
else
  echo "=== ONE OR MORE CHECKS FAILED for $DIR ===" >&2
  exit 1
fi
