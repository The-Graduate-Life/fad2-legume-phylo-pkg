#!/usr/bin/env bash
# Fetch FAD2 (omega-6 fatty acid desaturase) protein sequences from UniProtKB.
#
# Requires: curl, internet access to rest.uniprot.org (no API key needed).
# Runtime: a few seconds.
#
# Core set (used in the committed example results/ in this repo):
#   E9M5E3  AhFAD2  Arachis hypogaea (peanut), omega-6 fatty acid desaturase
#   Q5FBA1  GmFAD2-1A  Glycine max (soybean), seed-specific FAD2 paralog A
#   Q19AK8  GmFAD2-1B  Glycine max (soybean), seed-specific FAD2 paralog B
#
# Optional expansion set (NOT fetched by default -- pass -e/--extended to
# also pull this, which adds an outgroup outside the legumes):
#   P46313  AtFAD2   Arabidopsis thaliana, single-copy FAD2 (Okuley et al.
#                    1994, Plant Cell 6:147) -- classic outgroup for legume
#                    FAD2 comparisons.
#
# To broaden further, search UniProt (https://www.uniprot.org) for "FAD2"
# restricted to the species of interest, or use NCBI Entrez Direct on a
# machine with unrestricted internet access, e.g.:
#   esearch -db protein -query "FAD2[Gene] AND Cicer arietinum[Organism]" \
#     | efetch -format fasta >> data/fad2_seed.fasta
# Always record the accession(s) you add and re-run the pipeline from
# scratch rather than hand-editing the alignment/tree outputs.
#
set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data"
mkdir -p "$OUTDIR"

CORE_ACCESSIONS=(E9M5E3 Q5FBA1 Q19AK8)
EXTENDED_ACCESSIONS=(P46313)

ACCESSIONS=("${CORE_ACCESSIONS[@]}")
if [[ "${1:-}" == "-e" || "${1:-}" == "--extended" ]]; then
  ACCESSIONS+=("${EXTENDED_ACCESSIONS[@]}")
  echo "Fetching core + extended accession set (${#ACCESSIONS[@]} sequences)"
else
  echo "Fetching core accession set (${#ACCESSIONS[@]} sequences)."
  echo "Pass -e/--extended to also fetch an Arabidopsis outgroup sequence."
fi

for acc in "${ACCESSIONS[@]}"; do
  echo "  fetching ${acc} ..."
  curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
    || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
done

cat "${OUTDIR}"/*.fasta > "${OUTDIR}/fad2_seed.fasta"
echo "Wrote ${OUTDIR}/fad2_seed.fasta ($(grep -c '^>' "${OUTDIR}/fad2_seed.fasta") sequences)"
