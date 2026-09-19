#!/usr/bin/env bash
# Fetch FAD2 (omega-6 fatty acid desaturase) protein sequences.
#
# Two tiers, both always outgroup-rooted, all using pinned accessions
# fetched by plain curl -- no live search queries, no extra tooling
# beyond curl itself:
#   1. Core set (default, no flags) -- 3 legume sequences + 1 outgroup (4 total).
#   2. Legume expansion (-l/--legumes) -- adds 7 more legume species,
#      still + the same outgroup (11 total).
#
# The tree-building tools (MAFFT/FastTree/IQ-TREE) all produce UNROOTED
# trees by default -- without a taxon known to have diverged before every
# legume in the set, the root drawn in a figure is arbitrary, not a
# biological claim. A single non-legume outgroup (Arabidopsis) is
# therefore always included, so both tiers this pipeline produces are
# genuinely rooted. See run_pipeline.sh (and reroot_tree.py) for how the
# rooting itself is done.
#
# Requires: curl only.
#
# Runtime: a few seconds total for either tier.
#
# Usage:
#   bash scripts/fetch_data.sh                    # core + outgroup (4 seqs)
#   bash scripts/fetch_data.sh --legumes           # core + legume expansion + outgroup (11 seqs)
#
# Core set (used in the committed example results/ in this repo):
#   G8GUK7  AhFAD2     Arachis hypogaea (peanut)      [UniProtKB]
#   Q5FBA1  GmFAD2-1A  Glycine max (soybean), paralog A [UniProtKB]
#   Q19AK8  GmFAD2-1B  Glycine max (soybean), paralog B [UniProtKB]
#
# Note: the original pinned peanut accession (E9M5E3) was retired from
# UniProtKB after this project was first built ("not part of a reference
# proteome" -- confirmed via the UniProt search API, which reported it as
# an inactive/deleted entry with no direct successor accession). G8GUK7
# was selected as a replacement: full-length (~379 aa), annotated omega-6
# fatty acid desaturase (EC 1.3.1.35), with all three canonical membrane
# desaturase histidine-box motifs (HECGHH / HHSNT / HVAHH) intact. This is
# a real example of accession churn even with "pinned" identifiers --
# UniProtKB periodically retires entries that fall out of reference
# proteomes, so pinning still benefits from periodic re-verification.
#
# Legume expansion set (-l/--legumes):
#   A0A1S2Y1S8   Cicer arietinum      chickpea            [UniProtKB]
#   A0A4D6M8B4   Vigna unguiculata     cowpea              [UniProtKB]
#   XP_013467668.1  Medicago truncatula  barrel medic       [NCBI, by accession]
#   A0A151U3S7   Cajanus cajan          pigeon pea          [UniProtKB]
#   E2JFD6       Phaseolus vulgaris     common bean         [UniProtKB]
#   A0A9D4W617   Pisum sativum          garden pea          [UniProtKB]
#   A0A1S3V1H8   Vigna radiata           mung bean          [UniProtKB]
#
# All seven are pinned by explicit accession and fetched with a direct
# curl GET (for the NCBI entry, a single-accession efetch call, not a
# search) -- no query-based lookup, so results are exactly reproducible.
#
# Several other legume species were considered (lentil, alfalfa, red
# clover, faba bean, Lotus japonicus) but dropped: their UniProt hits were
# either absent or only "uncharacterized"/generically-named fatty acid
# desaturase domain proteins with no confident omega-6/FAD2-specific
# annotation, indistinguishable by search alone from the omega-3, acyl-ACP,
# or other paralogs also present in these genomes. Rather than guess,
# those species were left out -- see the topology caveat in the README's
# Results section for why picking the wrong paralog for even one species
# is enough to visibly distort a small gene tree.
#
# An earlier version of this script used NCBI Entrez Direct (esearch |
# efetch) to live-query each organism by name; that two-step search/
# history-server flow proved unreliable on at least one real network
# (repeated "SSL_read: unexpected eof" failures on 3 of 4 organism
# queries, reproducible across runs), so it was replaced with fixed
# accessions found once via manual UniProt/NCBI searches and pinned here.
# Selection criteria for each: full-length (~370-400 aa, matching the
# canonical microsomal FAD2 length), annotated as the omega-6/delta-12
# desaturase (explicitly excluding the omega-3/FAD3-family paralogs and
# the longer, chloroplast-targeted omega-6 paralogs that also turn up in
# these searches -- both are related but functionally distinct enzymes).
#
# Outgroup (always fetched, both tiers):
#   P46313   AtFAD2   Arabidopsis thaliana (thale cress)   [UniProtKB]
#
# Arabidopsis is not a legume, so including it breaks the "legumes only"
# scope on purpose: it diverged from all the legume species before they
# diverged from each other, making it a valid root for this comparison.
# Okuley et al. (1994, Plant Cell 6:147) is the original characterization
# of AtFAD2, which is why this specific accession is the standard choice
# as an outgroup for FAD2 comparisons in the literature.
set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data"
mkdir -p "$OUTDIR"

WANT_LEGUMES=0
for arg in "$@"; do
  case "$arg" in
    -l|--legumes) WANT_LEGUMES=1 ;;
    *) echo "WARNING: unrecognized argument '${arg}' ignored" >&2 ;;
  esac
done

CORE_ACCESSIONS=(G8GUK7 Q5FBA1 Q19AK8)
UNIPROT_LEGUME_ACCESSIONS=(A0A1S2Y1S8 A0A4D6M8B4 A0A151U3S7 E2JFD6 A0A9D4W617 A0A1S3V1H8)
NCBI_ACC="XP_013467668.1"
OUTGROUP_ACC="P46313"

echo "Fetching core accession set (${#CORE_ACCESSIONS[@]} sequences) from UniProt..."
for acc in "${CORE_ACCESSIONS[@]}"; do
  echo "  fetching ${acc} ..."
  curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
    || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
done

echo "Fetching outgroup (${OUTGROUP_ACC}, AtFAD2, Arabidopsis thaliana) from UniProt..."
curl -sS -f "https://rest.uniprot.org/uniprotkb/${OUTGROUP_ACC}.fasta" -o "${OUTDIR}/${OUTGROUP_ACC}.fasta" \
  || { echo "  WARNING: failed to fetch ${OUTGROUP_ACC}" >&2; rm -f "${OUTDIR}/${OUTGROUP_ACC}.fasta"; }

if [[ "$WANT_LEGUMES" -eq 1 ]]; then
  echo "Fetching legume expansion set (pinned accessions, UniProt + NCBI)..."

  for acc in "${UNIPROT_LEGUME_ACCESSIONS[@]}"; do
    echo "  fetching ${acc} (UniProt) ..."
    curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
      || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
  done

  echo "  fetching ${NCBI_ACC} (NCBI, by accession) ..."
  curl -sS -f "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${NCBI_ACC}&rettype=fasta&retmode=text" \
    -o "${OUTDIR}/${NCBI_ACC}.fasta" \
    || { echo "  WARNING: failed to fetch ${NCBI_ACC}" >&2; rm -f "${OUTDIR}/${NCBI_ACC}.fasta"; }

  echo
  echo "Legume expansion headers fetched:"
  grep "^>" "${OUTDIR}"/{A0A1S2Y1S8,A0A4D6M8B4,A0A151U3S7,E2JFD6,A0A9D4W617,A0A1S3V1H8,XP_013467668.1}.fasta 2>/dev/null || echo "  (some fetches may have failed -- see warnings above)"
fi

# Build seed files explicitly from the accession lists above -- NOT by
# globbing "${OUTDIR}/*.fasta". A glob-based cat is a trap here: once a
# seed file already exists in the same directory, it matches *.fasta too
# and gets folded back into the next concatenation, silently duplicating
# sequences. (This bit us once already -- a legume seed that should have
# had 10 sequences came out with 13, exactly 3 too many: the just-written
# core seed got re-included.)

# Core seed = core 3 + outgroup (4 sequences), always (re)written
# regardless of other flags, so run_pipeline.sh --both has a clean input
# on hand without requiring separate fetch_data.sh invocations.
: > "${OUTDIR}/fad2_seed_core.fasta"
for acc in "${CORE_ACCESSIONS[@]}"; do
  [[ -f "${OUTDIR}/${acc}.fasta" ]] && cat "${OUTDIR}/${acc}.fasta" >> "${OUTDIR}/fad2_seed_core.fasta"
done
[[ -f "${OUTDIR}/${OUTGROUP_ACC}.fasta" ]] && cat "${OUTDIR}/${OUTGROUP_ACC}.fasta" >> "${OUTDIR}/fad2_seed_core.fasta"

# Legacy/default combined seed (core+outgroup, or core+legumes+outgroup if
# this invocation also fetched the legume expansion). Kept for backwards
# compatibility with a plain `bash scripts/run_pipeline.sh` (no --both).
cp "${OUTDIR}/fad2_seed_core.fasta" "${OUTDIR}/fad2_seed.fasta"

if [[ "$WANT_LEGUMES" -eq 1 ]]; then
  : > "${OUTDIR}/fad2_seed_legumes.fasta"
  for acc in "${CORE_ACCESSIONS[@]}" "${UNIPROT_LEGUME_ACCESSIONS[@]}" "$NCBI_ACC"; do
    [[ -f "${OUTDIR}/${acc}.fasta" ]] && cat "${OUTDIR}/${acc}.fasta" >> "${OUTDIR}/fad2_seed_legumes.fasta"
  done
  [[ -f "${OUTDIR}/${OUTGROUP_ACC}.fasta" ]] && cat "${OUTDIR}/${OUTGROUP_ACC}.fasta" >> "${OUTDIR}/fad2_seed_legumes.fasta"
  cp "${OUTDIR}/fad2_seed_legumes.fasta" "${OUTDIR}/fad2_seed.fasta"
fi

echo
echo "Wrote ${OUTDIR}/fad2_seed_core.fasta ($(grep -c '^>' "${OUTDIR}/fad2_seed_core.fasta") sequences)"
if [[ -f "${OUTDIR}/fad2_seed_legumes.fasta" ]]; then
  echo "Wrote ${OUTDIR}/fad2_seed_legumes.fasta ($(grep -c '^>' "${OUTDIR}/fad2_seed_legumes.fasta") sequences)"
fi
