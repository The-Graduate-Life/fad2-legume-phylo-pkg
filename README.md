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

Compares three real FAD2 protein sequences:

| Sequence | UniProt accession | Organism | Role |
|---|---|---|---|
| AhFAD2 | [E9M5E3](https://www.uniprot.org/uniprotkb/E9M5E3) | *Arachis hypogaea* (peanut) | omega-6 FAD2 |
| GmFAD2-1A | [Q5FBA1](https://www.uniprot.org/uniprotkb/Q5FBA1) | *Glycine max* (soybean) | seed-specific FAD2 paralog A |
| GmFAD2-1B | [Q19AK8](https://www.uniprot.org/uniprotkb/Q19AK8) | *Glycine max* (soybean) | seed-specific FAD2 paralog B |

An optional Arabidopsis outgroup (UniProt P46313) can be added with one
flag (see below) to root the comparison; the pipeline itself places no
ceiling on the number of sequences -- swap in any FAD2 (or other
single-gene-family) FASTA and re-run.

**Note on scale:** three sequences is intentionally minimal ("light data"),
chosen to keep the whole run well under a minute and to make every step of
the Methods below auditable by hand. The result is still real and
biologically sensible (see Results): the two soybean paralogs are far more
similar to each other than either is to the peanut ortholog, consistent
with FAD2-1A/1B arising from a relatively recent within-lineage gene
duplication in soybean's paleopolyploid genome, while the peanut/soybean
split reflects much deeper legume divergence.

## Repository layout

```
.
├── README.md
├── environment.yml          # conda environment spec
├── .gitignore                # data/ is fetched, not committed
├── data/                     # (empty in git; populated by fetch_data.sh)
├── scripts/
│   ├── fetch_data.sh          # downloads sequences from UniProt (curl)
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

Everything runs from the command line; the only non-Python/R dependencies
are **MAFFT** and **FastTree**, both standard, actively-maintained
bioinformatics binaries.

**Option A -- conda (recommended, matches `environment.yml`):**

```bash
conda env create -f environment.yml
conda activate fad2-legume-phylo
```

**Option B -- apt (Debian/Ubuntu):**

```bash
sudo apt-get update
sudo apt-get install -y mafft fasttree python3-matplotlib curl
```

**Option C -- Homebrew (macOS):**

```bash
brew install mafft fasttree
pip install matplotlib
```

Versions used to produce the committed `results/`: MAFFT v7.505, FastTree
v2.1.11 (Double precision), Python 3.11 + matplotlib 3.10.

## Data fetch

Sequences are pulled from UniProtKB's REST API (no account or API key
needed):

```bash
bash scripts/fetch_data.sh            # core 3-sequence set
bash scripts/fetch_data.sh --extended # + Arabidopsis outgroup (4 sequences)
```

This writes `data/E9M5E3.fasta`, `data/Q5FBA1.fasta`, `data/Q19AK8.fasta`
(and optionally `data/P46313.fasta`), and concatenates them into
`data/fad2_seed.fasta`, which is what the pipeline consumes. Raw sequence
data is **not** committed to this repository (see `.gitignore`) -- re-run
the fetch script to regenerate it; the exact accessions are documented
both here and inline in `fetch_data.sh`, so the fetch is fully
reproducible.

To go beyond FAD2 in these three species, search UniProt or use NCBI
Entrez Direct (`esearch`/`efetch`) for additional accessions and append
them to `data/fad2_seed.fasta` before running the pipeline -- see the
comments in `scripts/fetch_data.sh` for a worked example command.

## Running

```bash
bash scripts/fetch_data.sh
bash scripts/run_pipeline.sh
```

Total wall-clock time for the 3-sequence set: a few seconds. Even scaled
up to a few dozen sequences (e.g. a full plant FAD2 gene family), MAFFT +
FastTree on a laptop should stay well within the 30-minute budget --
FastTree in particular scales near-linearly and is explicitly designed for
alignments with thousands of sequences.

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

## Results (this run)

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

## References

- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment
  software version 7. *Mol Biol Evol* 30(4):772-780.
- Price, M.N., Dehal, P.S., Arkin, A.P. (2010). FastTree 2 -- approximately
  maximum-likelihood trees for large alignments. *PLoS ONE* 5(3):e9490.
- UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase.
  *Nucleic Acids Res* 51(D1):D523-D531.
- Okuley, J. et al. (1994). Arabidopsis FAD2 gene encodes the enzyme
  essential for polyunsaturated lipid synthesis. *Plant Cell* 6(1):147-158.
