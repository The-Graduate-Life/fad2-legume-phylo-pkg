# FAD2 legume phylogenetics: a small, reproducible external-tool pipeline

A ~1-minute-to-run demonstration pipeline that fetches public protein
sequences for the **FAD2 (omega-6 / delta-12 fatty acid desaturase)** gene
family, aligns them with **MAFFT**, and builds a tree with **FastTree** --
producing a real (if small) phylogenetic result, not a toy dataset.

FAD2 converts oleic acid (18:1) to linoleic acid (18:2) and is the classic
target for "high-oleic" trait breeding in peanut and soybean (it's the same
gene targeted in peanut high-oleic marker-assisted selection and soybean
FAD2-1A/1B knockouts).

## What this run does

This is a **legume-only** comparison (no outgroup) of FAD2 protein
sequences, in two tiers:

**Core set** -- 3 sequences, pinned UniProt accessions, verified and
committed to `results/`:

| Sequence | UniProt accession | Organism | Role |
|---|---|---|---|
| AhFAD2 | [E9M5E3](https://www.uniprot.org/uniprotkb/E9M5E3) | *Arachis hypogaea* (peanut) | omega-6 FAD2 |
| GmFAD2-1A | [Q5FBA1](https://www.uniprot.org/uniprotkb/Q5FBA1) | *Glycine max* (soybean) | seed-specific FAD2 paralog A |
| GmFAD2-1B | [Q19AK8](https://www.uniprot.org/uniprotkb/Q19AK8) | *Glycine max* (soybean) | seed-specific FAD2 paralog B |

**Legume expansion** (optional, `-l`/`--legumes` flag) -- adds FAD2 from 4
more legume species, fetched live from NCBI rather than pinned by
accession (see [Data fetch](#data-fetch) for why that matters):

| Organism | Common name |
|---|---|
| *Cicer arietinum* | chickpea |
| *Vigna unguiculata* | cowpea |
| *Medicago truncatula* | barrel medic (legume model species) |
| *Cajanus cajan* | pigeon pea |

The pipeline itself places no ceiling on the number of sequences -- swap
in any FAD2 (or other single-gene-family) FASTA and re-run.

**Note on scale:** the core 3-sequence set is intentionally minimal
("light data"), chosen to keep the whole run well under a minute and to
make every step of the Methods below auditable by hand. The result is
still real and biologically sensible (see Results): the two soybean
paralogs are far more similar to each other than either is to the peanut
ortholog, consistent with FAD2-1A/1B arising from a relatively recent
within-lineage gene duplication in soybean's paleopolyploid genome, while
the peanut/soybean split reflects much deeper legume divergence. The
legume expansion scales this up to a small within-family comparison
without changing runtime meaningfully (see [Running](#running)).

## Repository layout

```
.
├── README.md
├── environment.yml          # conda environment spec
├── .gitignore                # data/ is fetched, not committed
├── data/                     # (empty in git; populated by fetch_data.sh)
├── scripts/
│   ├── fetch_data.sh          # core: curl -> UniProt; -l/--legumes: esearch/efetch -> NCBI
│   ├── run_pipeline.sh        # MAFFT -> FastTree -> stats -> figure
│   ├── pairwise_identity.py   # % identity matrix from the alignment
│   └── plot_tree.py           # renders the Newick tree as a PNG
└── results/                   # example outputs from this exact run
    ├── fad2_aligned.fasta
    ├── fad2_tree.nwk
    ├── fad2_tree.png
    ├── pairwise_identity.tsv
    ├── mafft.log
    └── fasttree.log
```

## Environment and installs

Everything runs from the command line. Core dependencies are **MAFFT** and
**FastTree**; the legume expansion additionally needs **NCBI Entrez
Direct** (`esearch`/`efetch`).

**conda (recommended -- this is what `environment.yml` declares, including
Entrez Direct, so the whole toolchain is reproducible from one file with
no `sudo`/system-package step):**

```bash
conda env create -f environment.yml   # first time
# or, to add entrez-direct to an existing env:
conda env update -f environment.yml
conda activate fad2-legume-phylo
```

If `conda` is slow to resolve dependencies, `mamba env update -f
environment.yml` is a drop-in faster solver.

Versions used to produce the committed `results/`: MAFFT v7.505, FastTree
v2.1.11 (Double precision), Python 3.11 + matplotlib 3.10. (Entrez Direct
is only exercised by the optional legume-expansion fetch, so it isn't
reflected in the committed core-set results.)

## Data fetch

**Core set** (3 sequences, pinned UniProt accessions, `curl` only):

```bash
bash scripts/fetch_data.sh
```

This writes `data/E9M5E3.fasta`, `data/Q5FBA1.fasta`, `data/Q19AK8.fasta`
and concatenates them into `data/fad2_seed.fasta`, which the pipeline
consumes. The exact accessions are documented both here and inline in
`fetch_data.sh`, so this half of the fetch is fully reproducible byte-for-byte.

**Legume expansion** (+4 species, live NCBI query, needs Entrez Direct):

```bash
bash scripts/fetch_data.sh --legumes   # or -l
```

This re-fetches the core 3, then queries NCBI's protein database for
`FAD2[Title] AND "<organism>"[Organism]` for each of the 4 additional
legumes and appends whatever comes back to `data/fad2_seed.fasta`. Unlike
the core set, **these are not pinned accessions** -- NCBI annotation
naming isn't fully standardized across species, so a query can return
zero, one, or several hits, and occasionally a partial/fragment sequence
or an off-target hit. The script prints every header it fetched for the
expansion; **read that list before running the pipeline**, and manually
remove (from `data/fad2_seed.fasta`) any entry that looks like a
duplicate, a fragment, or an unrelated protein. Record whatever you keep
or remove in your own lab notes / commit message, since a live query is
not by itself a reproducible citation the way a pinned accession is.

Raw sequence data is **not** committed to this repository (see
`.gitignore`) -- re-run the fetch script(s) to regenerate it.

To broaden further still, search UniProt directly or adjust the
`LEGUME_ORGANISMS` list in `scripts/fetch_data.sh`.

## Running

Core set only:

```bash
bash scripts/fetch_data.sh
bash scripts/run_pipeline.sh
```

With the legume expansion:

```bash
bash scripts/fetch_data.sh --legumes
# review the printed header list, edit data/fad2_seed.fasta if needed
bash scripts/run_pipeline.sh
```

Total wall-clock time for the 3-sequence core set: a few seconds. The
legume expansion adds a handful of NCBI round-trips (typically well under
a minute) plus a few more seconds of MAFFT/FastTree on 7 sequences instead
of 3 -- still trivial next to the 30-minute budget. FastTree in particular
scales near-linearly and is explicitly designed for alignments with
thousands of sequences, so this pipeline has plenty of headroom to grow
well past a full plant FAD2 gene family without approaching that budget.

## Methods (as you'd write it up)

Protein sequences for the FAD2 (omega-6 fatty acid desaturase) gene family
were retrieved from UniProtKB (accessions E9M5E3, Q5FBA1, Q19AK8; UniProt
Consortium) via the UniProt REST API. Sequences were aligned with MAFFT
v7.505 (Katoh & Standley, 2013) using the `--auto` strategy, which selected
the L-INS-i algorithm (iterative refinement with local pairwise alignment
information) appropriate for this sequence set. An approximate
maximum-likelihood phylogeny was inferred from the resulting alignment
with FastTree v2.1.11 (Price et al., 2010) under the Jones-Taylor-Thornton
(JTT) amino acid substitution model with CAT-approximated among-site rate
variation (20 categories); branch lengths were further optimized under
this model and support was estimated with FastTree's SH-like local
support test. Pairwise percent sequence identity was computed directly
from the MAFFT alignment (identical residues / aligned, non-gap columns).
The resulting tree was rendered with a minimal custom Newick parser and
matplotlib (no third-party tree-plotting library required).

## Results

**Core set (committed in `results/`, this repo's checked-in example run):**

Pairwise percent identity (`results/pairwise_identity.tsv`):

|  | AhFAD2 (peanut) | GmFAD2-1A | GmFAD2-1B |
|---|---|---|---|
| AhFAD2 (peanut) | 100.0 | 83.6 | 83.3 |
| GmFAD2-1A | 83.6 | 100.0 | 93.8 |
| GmFAD2-1B | 83.3 | 93.8 | 100.0 |

FastTree branch lengths (`results/fad2_tree.nwk`, substitutions/site):
AhFAD2 = 0.171, GmFAD2-1A = 0.029, GmFAD2-1B = 0.037 -- i.e. the two
soybean paralogs sit on much shorter branches from their common node than
either does from the peanut sequence, matching the ~94% vs ~83% identity
split above. See `results/fad2_tree.png` for the rendered tree.

**Legume expansion:** not yet run in this checked-in example -- the
sequences depend on a live NCBI query (see [Data fetch](#data-fetch)), so
they aren't pinned/reproducible enough to commit as a fixed "expected"
result the way the core set is. Run `bash scripts/fetch_data.sh --legumes
&& bash scripts/run_pipeline.sh` to generate it locally; if you want that
expanded tree/table committed to this repo as a second example, add it
under `results/` (e.g. `results/legume_expansion/`) alongside a note of
exactly which accessions/headers ended up in `data/fad2_seed.fasta`, so
the Methods section can cite it precisely.

## References

- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment
  software version 7. *Mol Biol Evol* 30(4):772-780.
- Price, M.N., Dehal, P.S., Arkin, A.P. (2010). FastTree 2 -- approximately
  maximum-likelihood trees for large alignments. *PLoS ONE* 5(3):e9490.
- UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase.
  *Nucleic Acids Res* 51(D1):D523-D531.
- Okuley, J. et al. (1994). Arabidopsis FAD2 gene encodes the enzyme
  essential for polyunsaturated lipid synthesis. *Plant Cell* 6(1):147-158.
