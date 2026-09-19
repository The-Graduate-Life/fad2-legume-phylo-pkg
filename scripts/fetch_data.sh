#!/usr/bin/env bash
# Fetch FAD2 (omega-6 fatty acid desaturase) protein sequences.
#
# Two tiers:
#   1. Core set (default, no flags) -- 3 verified UniProtKB accessions,
#      fetched by curl. This is what results/ in this repo was built from.
#   2. Legume expansion (-l/--legumes) -- adds 4 more legume species,
#      fetched dynamically from NCBI via Entrez Direct (esearch/efetch).
#      These are NOT pinned accessions: NCBI's annotation naming isn't
#      fully standardized, so a query can return zero, one, or several
#      hits per species, and occasionally a partial/fragment sequence.
#      ALWAYS inspect the header list this script prints at the end
#      before running the pipeline on the result -- see the note below.
#
# Requires:
#   - curl (core set)
#   - NCBI Entrez Direct, i.e. esearch/efetch on PATH (legume expansion
#     only). Install reproducibly via conda/mamba -- see environment.yml
#     (bioconda::entrez-direct) -- rather than system apt/sudo, so the
#     whole environment stays declared in one file.
#
# Runtime: a few seconds for the core set; the legume expansion adds a
# handful of NCBI round-trips, typically well under a minute total.
#
# Core set (used in the committed example results/ in this repo):
#   E9M5E3  AhFAD2     Arachis hypogaea (peanut)
#   Q5FBA1  GmFAD2-1A  Glycine max (soybean), seed-specific paralog A
#   Q19AK8  GmFAD2-1B  Glycine max (soybean), seed-specific paralog B
#
# Legume expansion set (-l/--legumes):
#   Cicer arietinum      chickpea
#   Vigna unguiculata     cowpea
#   Medicago truncatula   barrel medic (legume model species)
#   Cajanus cajan          pigeon pea
#
# An earlier version of this script also offered an Arabidopsis (non-legume)
# outgroup via -e/--extended, for rooting the tree. That is intentionally
# dropped here -- this project's scope is a legume-only, unrooted
# comparison. To add an outgroup back, fetch UniProt accession P46313
# (AtFAD2, Arabidopsis thaliana; Okuley et al. 1994, Plant Cell 6:147)
# the same way the core set is fetched below.
set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data"
mkdir -p "$OUTDIR"

CORE_ACCESSIONS=(E9M5E3 Q5FBA1 Q19AK8)

echo "Fetching core accession set (${#CORE_ACCESSIONS[@]} sequences) from UniProt..."
for acc in "${CORE_ACCESSIONS[@]}"; do
  echo "  fetching ${acc} ..."
  curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
    || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
done

if [[ "${1:-}" == "-l" || "${1:-}" == "--legumes" ]]; then
  if ! command -v esearch >/dev/null 2>&1 || ! command -v efetch >/dev/null 2>&1; then
    echo "ERROR: esearch/efetch not found on PATH." >&2
    echo "Install NCBI Entrez Direct reproducibly via conda:" >&2
    echo "  conda env update -f environment.yml && conda activate fad2-legume-phylo" >&2
    exit 1
  fi

  LEGUME_ORGANISMS=(
    "Cicer arietinum"
    "Vigna unguiculata"
    "Medicago truncatula"
    "Cajanus cajan"
  )

  echo "Fetching legume expansion set from NCBI (live query, not pinned accessions)..."
  LEGUME_FASTA="${OUTDIR}/legume_expansion.fasta"
  : > "$LEGUME_FASTA"
  for organism in "${LEGUME_ORGANISMS[@]}"; do
    echo "  querying FAD2 for ${organism} ..."
    esearch -db protein -query "FAD2[Title] AND \"${organism}\"[Organism]" \
      | efetch -format fasta >> "$LEGUME_FASTA" \
      || echo "  WARNING: query failed for ${organism}" >&2
  done

  echo
  echo "Legume expansion headers fetched -- REVIEW THESE before running the pipeline:"
  grep "^>" "$LEGUME_FASTA" || echo "  (none returned)"
  echo
  echo "If any entry above looks like a duplicate, a fragment, or an"
  echo "unrelated/mis-annotated protein, edit ${LEGUME_FASTA} to remove it"
  echo "before proceeding."
fi

cat "${OUTDIR}"/*.fasta > "${OUTDIR}/fad2_seed.fasta" 2>/dev/null || true
echo
echo "Wrote ${OUTDIR}/fad2_seed.fasta ($(grep -c '^>' "${OUTDIR}/fad2_seed.fasta") sequences)"
