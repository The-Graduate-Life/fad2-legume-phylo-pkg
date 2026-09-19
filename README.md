# FAD2 legume phylogenetics: a small, reproducible external-tool pipeline

A demonstration pipeline that fetches public protein sequences for the
**FAD2 (omega-6 / delta-12 fatty acid desaturase)** gene family, aligns
them with **MAFFT**, and builds trees with **FastTree** and **IQ-TREE**
(ModelFinder + 1000 ultrafast bootstrap replicates) -- producing real,
statistically supported phylogenetic results, not a toy dataset.

FAD2 converts oleic acid (18:1) to linoleic acid (18:2) and is the classic
target for "high-oleic" trait breeding in peanut and soybean (it's the same
gene targeted in peanut high-oleic marker-assisted selection and soybean
FAD2-1A/1B knockouts).

## What this run does

This is a **legume-focused** comparison of FAD2 protein sequences, in two
tiers. Both are rooted on the same non-legume outgroup, Arabidopsis: MAFFT
/FastTree/IQ-TREE all build **unrooted** trees by default, and without a
taxon known to have diverged before every legume in the set, the "root"
drawn in a figure is just wherever the algorithm happened to place it, not
a claim about ancestry. Arabidopsis is not a legume, so it diverged from
every legume species here before they diverged from each other, making it
a valid, biologically meaningful root.

**Core set** -- 4 sequences (3 legumes + 1 outgroup), pinned UniProt
accessions, verified and committed to `results/`:

| Sequence | UniProt accession | Organism | Role |
|---|---|---|---|
| AhFAD2 | [G8GUK7](https://www.uniprot.org/uniprotkb/G8GUK7) | *Arachis hypogaea* (peanut) | omega-6 FAD2 |
| GmFAD2-1A | [Q5FBA1](https://www.uniprot.org/uniprotkb/Q5FBA1) | *Glycine max* (soybean) | seed-specific FAD2 paralog A |
| GmFAD2-1B | [Q19AK8](https://www.uniprot.org/uniprotkb/Q19AK8) | *Glycine max* (soybean) | seed-specific FAD2 paralog B |
| AtFAD2 | [P46313](https://www.uniprot.org/uniprotkb/P46313) | *Arabidopsis thaliana* (thale cress) | outgroup (not a legume) |

*P46313 is the standard reference outgroup choice for FAD2 comparisons: Okuley et al.
(1994) is the original functional characterization of Arabidopsis FAD2.)*

**Legume expansion** (optional, `-l`/`--legumes` flag) -- adds FAD2 from 7
more legume species, each pinned to a specific accession and fetched with
plain `curl`, alongside the same core set and outgroup (11 sequences
total; see [Data fetch](#data-fetch) for how these were chosen):

| Organism | Common name | Accession | Source |
|---|---|---|---|
| *Cicer arietinum* | chickpea | [A0A1S2Y1S8](https://www.uniprot.org/uniprotkb/A0A1S2Y1S8) | UniProtKB |
| *Vigna unguiculata* | cowpea | [A0A4D6M8B4](https://www.uniprot.org/uniprotkb/A0A4D6M8B4) | UniProtKB |
| *Medicago truncatula* | barrel medic (legume model species) | [XP_013467668.1](https://www.ncbi.nlm.nih.gov/protein/XP_013467668.1) | NCBI (by accession) |
| *Cajanus cajan* | pigeon pea | [A0A151U3S7](https://www.uniprot.org/uniprotkb/A0A151U3S7) | UniProtKB |
| *Phaseolus vulgaris* | common bean | [E2JFD6](https://www.uniprot.org/uniprotkb/E2JFD6) | UniProtKB |
| *Pisum sativum* | garden pea | [A0A9D4W617](https://www.uniprot.org/uniprotkb/A0A9D4W617) | UniProtKB |
| *Vigna radiata* | mung bean | [A0A1S3V1H8](https://www.uniprot.org/uniprotkb/A0A1S3V1H8) | UniProtKB |

The pipeline itself places no ceiling on the number of sequences -- swap
in any FAD2 (or other single-gene-family) FASTA and re-run.

**Note on scale:** the core 4-sequence set is intentionally minimal
("light data"), chosen to keep the whole run well under a minute and to
make every step of the Methods below auditable by hand. The result is
still real and biologically sensible (see Results): the two soybean
paralogs are far more similar to each other than either is to the peanut
ortholog, consistent with FAD2-1A/1B arising from a relatively recent
within-lineage gene duplication in soybean's paleopolyploid genome, while
the peanut/soybean split reflects much deeper legume divergence, and the
Arabidopsis outgroup is more diverged still. The legume expansion scales
this up to an 11-taxon comparison with real bootstrap support, still
trivial against the 30-minute budget (see [Running](#running)).

## Repository layout

```
.
├── README.md
├── environment.yml          # conda environment spec
├── .gitignore                # data/ is fetched, not committed
├── data/                     # (empty in git; populated by fetch_data.sh)
├── scripts/
│   ├── fetch_data.sh          # curl -> pinned accessions (UniProt + one NCBI); always includes the outgroup; -l/--legumes adds 7 more species
│   ├── run_pipeline.sh        # MAFFT -> FastTree + IQ-TREE -> reroot on outgroup -> stats -> figures; --both runs both tiers
│   ├── reroot_tree.py         # outgroup rerooting for a Newick tree (splits the outgroup's branch, reverses edges back to the old root)
│   ├── pairwise_identity.py   # % identity matrix from the alignment
│   └── plot_tree.py           # renders a Newick tree (with support labels) as a PNG
└── results/                   # example outputs: core set (4 seqs, rooted)
    ├── fad2_aligned.fasta
    ├── fad2_tree.nwk             # FastTree, raw (unrooted) output
    ├── fad2_tree_rooted.nwk      # FastTree, rerooted on the outgroup
    ├── fad2_tree_rooted.png
    ├── fad2_tree_iqtree.nwk      # IQ-TREE, rerooted on the outgroup (see gotcha below)
    ├── fad2_tree_iqtree.png
    ├── fad2_iqtree.iqtree        # full IQ-TREE report (model selection, etc.)
    ├── pairwise_identity.tsv
    ├── mafft.log
    ├── fasttree.log
    ├── iqtree.log
    └── legume_expansion/       # second example: full 11-sequence run
        ├── fad2_aligned.fasta
        ├── fad2_tree.nwk
        ├── fad2_tree_rooted.nwk
        ├── fad2_tree_rooted.png
        ├── fad2_tree_iqtree.nwk
        ├── fad2_tree_iqtree.png
        ├── fad2_iqtree.iqtree
        ├── pairwise_identity.tsv
        ├── mafft.log
        ├── fasttree.log
        └── iqtree.log
```

## Environment and installs

Everything runs from the command line. Dependencies are **MAFFT**,
**FastTree**, **IQ-TREE**, and `curl` -- no extra tooling is needed
beyond these for either data fetch or tree-building.

**conda (recommended -- this is what `environment.yml` declares, so the
whole toolchain is reproducible from one file with no `sudo`/system-package
step):**

```bash
conda env create -f environment.yml   # first time
# or, to update an existing env after a change to environment.yml:
conda env update -f environment.yml
conda activate fad2-legume-phylo
```

If `conda` is slow to resolve dependencies, `mamba env update -f
environment.yml` is a drop-in faster solver.

Versions used to produce the committed `results/`: MAFFT v7.526, FastTree
v2.2.0 (Double precision), IQ-TREE v3.1.3, Python 3.11.16 + matplotlib.
These are exactly what `environment.yml` pins (`mafft=7.526`,
`fasttree=2.2.0`, `iqtree=3.1.3`, `python=3.11`) -- confirmed by running
`mafft --version`, `fasttree -help`, and `iqtree --version` in the active
environment before generating the committed results.

**A naming gotcha worth knowing:** the bioconda `iqtree` package's binary
is called `iqtree`, not `iqtree2` -- true even for IQ-TREE 2.x/3.x. There
is no `iqtree2` binary on PATH after installing this package. The scripts
here already call `iqtree` correctly; this note exists because the first
draft of `run_pipeline.sh` called `iqtree2` and failed with "command not
found."


## Data fetch

**Core set** (4 sequences: 3 legumes + the Arabidopsis outgroup, pinned
UniProt accessions, `curl` only):

```bash
bash scripts/fetch_data.sh
```

This writes `data/G8GUK7.fasta`, `data/Q5FBA1.fasta`, `data/Q19AK8.fasta`,
and `data/P46313.fasta` (the outgroup), concatenating the first three plus
the outgroup into `data/fad2_seed_core.fasta` (and, for backwards
compatibility, `data/fad2_seed.fasta`). The exact accessions are
documented both here and inline in `fetch_data.sh`, so this fetch is fully
reproducible byte-for-byte.

The outgroup accession, P46313 (AtFAD2, *Arabidopsis thaliana*), was
chosen because it's the standard reference for this exact purpose: Okuley
et al. (1994, *Plant Cell* 6:147) is the original functional
characterization of Arabidopsis FAD2, and this accession recurs across the
literature whenever an outgroup is needed for FAD2 comparisons.

**Legume expansion** (+7 species, pinned accessions, `curl` only):

```bash
bash scripts/fetch_data.sh --legumes   # or -l
```

This fetches the core 3 legumes, the outgroup, and 7 additional legume
FAD2 sequences by explicit accession (6 from UniProtKB, 1 from NCBI by
direct accession lookup -- see the table in
[What this run does](#what-this-run-does)), writing
`data/fad2_seed_legumes.fasta` (11 sequences) alongside the
core-only `data/fad2_seed_core.fasta` (4 sequences).

Raw sequence data is **not** committed to this repository (see
`.gitignore`) -- re-run the fetch script(s) to regenerate it.

To broaden further still, add more pinned accessions to
`scripts/fetch_data.sh` the same way.

## Running

Core set only (default):

```bash
bash scripts/fetch_data.sh
bash scripts/run_pipeline.sh
```

Both tiers in one command:

```bash
bash scripts/fetch_data.sh --legumes
bash scripts/run_pipeline.sh --both
```

This builds the 4-sequence core set into `results/` and the 11-sequence
legume expansion into `results/legume_expansion/`, one after the other.

For each input set, `run_pipeline.sh` builds **two rooted trees**:
FastTree (fast, approximate, SH-like support) and IQ-TREE (ModelFinder
model selection + 1000 ultrafast bootstrap replicates), both rerooted on
the Arabidopsis outgroup via `reroot_tree.py` (see the `-o`-is-cosmetic
-only gotcha in [Environment and installs](#environment-and-installs) for
why that script is needed for IQ-TREE's output too, not just FastTree's).

## Methods

Protein sequences for the FAD2 (omega-6 fatty acid desaturase) gene family
were retrieved from UniProtKB (accessions G8GUK7, Q5FBA1, Q19AK8; UniProt
Consortium) via the UniProt REST API for the core set, with an optional
expansion of 7 further legume species (UniProtKB: A0A1S2Y1S8, A0A4D6M8B4,
A0A151U3S7, E2JFD6, A0A9D4W617, A0A1S3V1H8; NCBI: XP_013467668.1), all
fetched by pinned accession. Each expansion accession was selected from
several candidate isoforms per species by matching canonical microsomal
FAD2 length (~370-400 aa) and correct omega-6/delta-12 desaturase
annotation, explicitly excluding longer chloroplast-targeted paralogs and
the omega-3/FAD3-family desaturases that also appear in these searches;
strict 1:1 orthology across species was not independently confirmed
beyond this (see the caveat in Results). One additional, non-legume
sequence was included in both tiers as an outgroup: AtFAD2 from
*Arabidopsis thaliana* (UniProtKB P46313; Okuley et al., 1994), chosen
because Arabidopsis diverged from every legume species here before they
diverged from each other, giving the otherwise-unrooted tree a
biologically meaningful root. Sequences were aligned with MAFFT v7.526
(Katoh & Standley, 2013) using the `--auto` strategy, which selected the
L-INS-i algorithm (iterative refinement with local pairwise alignment
information) appropriate for this sequence set.

Two phylogenies were inferred from each resulting alignment. (1) An
approximate maximum-likelihood tree with FastTree v2.2.0 (Price et al.,
2010) under the Jones-Taylor-Thornton (JTT) amino acid substitution model
with CAT-approximated among-site rate variation (20 categories); branch
lengths were further optimized under this model and support was estimated
with FastTree's SH-like local support test. (2) A maximum-likelihood tree
with IQ-TREE v3.1.3 (Nguyen et al., 2015), with substitution model
selected independently for each tier by ModelFinder (Kalyaanamoorthy et
al., 2017) under the Bayesian information criterion (best-fit model:
Q.PLANT+G4 for the 4-sequence core set, Q.PLANT+I+G4 for the 11-sequence
legume-expansion set -- both a plant-specific empirical amino acid model
with 4-category gamma-distributed rate variation, the expansion set's
model additionally estimating a proportion of invariable sites) and
branch support from 1000 ultrafast bootstrap replicates (Hoang et al.,
2018). IQ-TREE was run with a fixed random seed (`-seed 42`) and forced
to a single thread (`-nt 1`, not `-nt AUTO`), since IQ-TREE's own
documentation notes that a multi-threaded run is not guaranteed to
reproduce identical output even with the same seed, and `-nt AUTO` would
itself pick a different thread count on a different machine -- both
matter for exact reproducibility, which this repository is set up to
support (see `FOLLOWME.md`). Both the FastTree and IQ-TREE trees for each tier were rerooted on
the outgroup with a standard outgroup-rooting procedure (insert a new
root on the branch leading to the outgroup, splitting that branch's
length in half, implemented in `scripts/reroot_tree.py`) -- necessary for
FastTree, which has no built-in rooting option, and, less obviously, also
for IQ-TREE, whose `-o` flag only affects which taxon its own report
draws at the root rather than actually rerooting its `.treefile` output
(see the gotcha in [Environment and installs](#environment-and-installs)).

Pairwise percent sequence identity was computed directly from the MAFFT
alignment (identical residues / aligned, non-gap columns). Trees were
rendered with a minimal custom Newick parser and matplotlib (no
third-party tree-plotting library required), with internal nodes labeled
by their support value (FastTree SH-like support or IQ-TREE ultrafast
bootstrap percentage, as applicable).

## Results

**Core set (`results/`, 4 sequences: 3 legumes + Arabidopsis outgroup):**

Pairwise percent identity (`results/pairwise_identity.tsv`):

|  | AhFAD2 (peanut) | GmFAD2-1A | GmFAD2-1B | AtFAD2 (Arabidopsis) |
|---|---|---|---|---|
| AhFAD2 (peanut) | 100.0 | 83.9 | 83.9 | 71.5 |
| GmFAD2-1A | 83.9 | 100.0 | 93.8 | 69.8 |
| GmFAD2-1B | 83.9 | 93.8 | 100.0 | 69.3 |
| AtFAD2 (Arabidopsis) | 71.5 | 69.8 | 69.3 | 100.0 |

Arabidopsis is clearly the outgroup by identity alone -- the least similar
sequence to every legume here (69.3-71.5%), well below the legume-to-legume
range (83.9-93.8%). Both trees confirm this structurally: a genuine
bifurcation at the root (not the 3-branch trifurcation an unrooted tree
would show), with Arabidopsis alone on one side.

*FastTree topology, rooted* (`results/fad2_tree_rooted.nwk`, SH-like
support shown):

```
(Arabidopsis,
  (peanut, (soybean-1A, soybean-1B)) [0.998]
);
```

*IQ-TREE topology, rooted* (`results/fad2_tree_iqtree.nwk`, ultrafast
bootstrap % shown):

```
(Arabidopsis,
  (peanut, (soybean-1A, soybean-1B)) [100]
);
```

Both methods agree completely and with maximal confidence (SH-like 0.998,
bootstrap 100%): the two soybean paralogs are sister to each other, with
peanut as the next-closest relative, matching the ~94% vs ~84% identity
split above -- consistent with FAD2-1A/1B arising from a comparatively
recent within-soybean duplication, while the peanut/soybean split and the
Arabidopsis/legume split reflect progressively deeper divergence. IQ-TREE's
best-fit model for this tier was Q.PLANT+G4 (see
[Methods](#methods)). With only three ingroup taxa
plus one outgroup there are just three possible rooted topologies to
choose between, so this result -- while simple -- is a genuine,
well-supported statistical statement, not a foregone conclusion. See
`results/fad2_tree_rooted.png` and `results/fad2_tree_iqtree.png` for the
rendered trees.

**Legume expansion (`results/legume_expansion/`, 11 sequences: 10 legumes
+ Arabidopsis outgroup):**

Pairwise percent identity ranges from 67.7% (pigeon pea vs. mung bean) to
93.8% (the two soybean paralogs, as above) among the legumes; Arabidopsis
is again the least similar sequence to every legume (67.5-76.0%),
consistently below the legume-to-legume range -- full matrix in
`results/legume_expansion/pairwise_identity.tsv`.

*FastTree topology, rooted* (`results/legume_expansion/fad2_tree_rooted.nwk`,
SH-like support shown):

```
(Arabidopsis,
  (
    (cowpea, common bean) [0.993],
    (pigeon pea,
      (
        (peanut, (chickpea, (mung bean, (soybean-1A, soybean-1B) [0.949]) [0.969]) [0.911]) [1.000],
        (garden pea, Medicago) [0.996]
      ) [0.435]
    ) [0.522]
  )
);
```

*IQ-TREE topology, rooted* (`results/legume_expansion/fad2_tree_iqtree.nwk`,
ModelFinder-selected model Q.PLANT+I+G4, ultrafast bootstrap % shown):

```
(Arabidopsis,
  (
    (
      (peanut, (chickpea, (mung bean, (soybean-1A, soybean-1B) [89]) [96]) [87]) [100],
      (cowpea, common bean) [100]
    ) [34],
    (pigeon pea, (garden pea, Medicago) [100]) [55]
  )
);
```

Both trees agree, with high confidence, on three things: cowpea + common
bean are sisters (SH-like 0.993, bootstrap 100%); garden pea + Medicago
are sisters (0.996 / 100%); and peanut, chickpea, mung bean, and the
soybean pair form a single nested clade in the *same* branching order in
both trees (soybean pair innermost, peanut outermost) with strong support
throughout (0.911-1.000 SH-like; 87-100% bootstrap).

Where they disagree is how pigeon pea and the garden-pea/Medicago pair
attach to the rest of the tree. FastTree splits pigeon pea off on its own
(SH-like 0.522) before separately resolving the peanut-clade/garden-pea
-Medicago split (0.435); IQ-TREE instead groups pigeon pea with
garden-pea/Medicago as a single clade, sister to everything else -- but
at only **34% bootstrap support**, i.e. essentially no confidence, well
below the ~70% threshold often used as a rule-of-thumb lower bound for
trusting a node. This is a real result worth stating plainly rather than
smoothing over: each species' sequence was chosen independently (closest
length + correct omega-6 desaturase annotation to the canonical FAD2, per
species -- see [Methods](#methods)), without
confirming that all ten legume sequences are the same 1:1 orthologous
copy. Species with multiple FAD2-family paralogs (most of them, per the
UniProt searches run during data fetch) can easily place a given
comparison one duplication node off from the "true" ortholog set, which
is enough to weaken support at exactly the deepest, most ambiguous nodes
-- which is what the 34% bootstrap value above suggests is happening.
Chickpea and mung bean's placement relative to the soybean pair also
isn't a clean match to expected Fabaceae tribal relationships (mung bean,
cowpea, common bean, and pigeon pea are all Phaseoleae, but do not form a
single clade in either tree). A more rigorous version of this comparison
would build a small gene family tree per species first (all isoforms, not
one pick per species) and identify the orthologous clade before building
the cross-species comparison -- a natural next step if this pipeline were
extended further.

Rooting on Arabidopsis clarifies the direction of evolution (which end of
the tree is ancestral) but does not, and should not be expected to,
resolve this underlying ambiguity -- rooting and paralog-confirmation are
different problems, and only the second would fix this in a rebuild. See
`results/legume_expansion/fad2_tree_rooted.png` and
`results/legume_expansion/fad2_tree_iqtree.png` for the rendered trees,
and `results/legume_expansion/fad2_iqtree.iqtree` for IQ-TREE's full
report (model selection table, log-likelihood, etc.).

## References

- Katoh, K. & Standley, D.M. (2013). MAFFT multiple sequence alignment
  software version 7. *Mol Biol Evol* 30(4):772-780.
- Price, M.N., Dehal, P.S., Arkin, A.P. (2010). FastTree 2 -- approximately
  maximum-likelihood trees for large alignments. *PLoS ONE* 5(3):e9490.
- Nguyen, L.T., Schmidt, H.A., von Haeseler, A., Minh, B.Q. (2015).
  IQ-TREE: a fast and effective stochastic algorithm for estimating
  maximum-likelihood phylogenies. *Mol Biol Evol* 32(1):268-274.
- Kalyaanamoorthy, S., Minh, B.Q., Wong, T.K.F., von Haeseler, A.,
  Jermiin, L.S. (2017). ModelFinder: fast model selection for accurate
  phylogenetic estimates. *Nat Methods* 14:587-589.
- Hoang, D.T., Chernomor, O., von Haeseler, A., Minh, B.Q., Vinh, L.S.
  (2018). UFBoot2: improving the ultrafast bootstrap approximation.
  *Mol Biol Evol* 35(2):518-522.
- UniProt Consortium (2023). UniProt: the Universal Protein Knowledgebase.
  *Nucleic Acids Res* 51(D1):D523-D531.
- Okuley, J. et al. (1994). Arabidopsis FAD2 gene encodes the enzyme
  essential for polyunsaturated lipid synthesis. *Plant Cell* 6(1):147-158.
