#!/usr/bin/env bash
# Fetch FAD2 (omega-6 fatty acid desaturase) protein sequences.
#
# Two tiers, both using pinned accessions fetched by plain curl -- no
# live search queries, no extra tooling beyond curl itself:
#   1. Core set (default, no flags) -- 3 sequences.
#   2. Legume expansion (-l/--legumes) -- adds 4 more legume species.
#
# Requires: curl only.
#
# Runtime: a few seconds total for either tier.
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
#
# All four are pinned by explicit accession and fetched with a direct
# curl GET (for the NCBI entry, a single-accession efetch call, not a
# search) -- no query-based lookup, so results are exactly reproducible.
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
# An earlier version of this script also offered an Arabidopsis (non-legume)
# outgroup via -e/--extended, for rooting the tree. That is intentionally
# dropped here -- this project's scope is a legume-only, unrooted
# comparison. To add an outgroup back, fetch UniProt accession P46313
# (AtFAD2, Arabidopsis thaliana; Okuley et al. 1994, Plant Cell 6:147)
# the same way the core set is fetched below.
set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/data"
mkdir -p "$OUTDIR"

CORE_ACCESSIONS=(G8GUK7 Q5FBA1 Q19AK8)

echo "Fetching core accession set (${#CORE_ACCESSIONS[@]} sequences) from UniProt..."
for acc in "${CORE_ACCESSIONS[@]}"; do
  echo "  fetching ${acc} ..."
  curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
    || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
done

if [[ "${1:-}" == "-l" || "${1:-}" == "--legumes" ]]; then
  echo "Fetching legume expansion set (pinned accessions, UniProt + NCBI)..."

  UNIPROT_LEGUME_ACCESSIONS=(A0A1S2Y1S8 A0A4D6M8B4 A0A151U3S7)
  for acc in "${UNIPROT_LEGUME_ACCESSIONS[@]}"; do
    echo "  fetching ${acc} (UniProt) ..."
    curl -sS -f "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "${OUTDIR}/${acc}.fasta" \
      || { echo "  WARNING: failed to fetch ${acc}" >&2; rm -f "${OUTDIR}/${acc}.fasta"; }
  done

  NCBI_ACC="XP_013467668.1"
  echo "  fetching ${NCBI_ACC} (NCBI, by accession) ..."
  curl -sS -f "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${NCBI_ACC}&rettype=fasta&retmode=text" \
    -o "${OUTDIR}/${NCBI_ACC}.fasta" \
    || { echo "  WARNING: failed to fetch ${NCBI_ACC}" >&2; rm -f "${OUTDIR}/${NCBI_ACC}.fasta"; }

  echo
  echo "Legume expansion headers fetched:"
  grep "^>" "${OUTDIR}"/{A0A1S2Y1S8,A0A4D6M8B4,A0A151U3S7,XP_013467668.1}.fasta 2>/dev/null || echo "  (some fetches may have failed -- see warnings above)"
fi

cat "${OUTDIR}"/*.fasta > "${OUTDIR}/fad2_seed.fasta" 2>/dev/null || true
echo
echo "Wrote ${OUTDIR}/fad2_seed.fasta ($(grep -c '^>' "${OUTDIR}/fad2_seed.fasta") sequences)"
