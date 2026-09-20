# FOLLOWME: running this pipeline from scratch

A step-by-step walkthrough for setting up the environment, fetching data,
running the pipeline, and checking the results. See `README.md` for the
full background, methods, and results -- this file is just the "what do I
type" version.

## 1. Get the environment set up

From the repository root:

```bash
conda env create -f environment.yml   # first time only
conda activate fad2-legume-phylo
```

If you already have the environment and `environment.yml` has changed
since:

```bash
conda env update -f environment.yml
conda activate fad2-legume-phylo
```

`mamba env update -f environment.yml` is a drop-in faster solver if
`conda` is slow to resolve.

**Confirm the toolchain is actually on PATH and the right versions:**

```bash
mafft --version
fasttree -help | head -1
iqtree --version
python3 --version
```

Expected: MAFFT v7.526, FastTree v2.2.0 (Double precision), IQ-TREE
v3.1.3, Python 3.11.x. If `iqtree --version` says "command not found,"
double check you activated the `fad2-legume-phylo` environment -- there
is no `iqtree2` binary, only `iqtree` (see the gotcha in `README.md` if
this trips you up).

## 2. Fetch the sequence data

Everything is pulled by pinned accession over plain `curl` -- no live
search queries, so this step is exactly reproducible.

**Core set only** (4 sequences: peanut, 2 soybean paralogs, Arabidopsis
outgroup):

```bash
bash scripts/fetch_data.sh
```

**Core set + legume expansion** (11 sequences total -- adds 7 more
legume species):

```bash
bash scripts/fetch_data.sh --legumes
```

Either command writes into `data/` (gitignored, so you'll need to re-run
this after a fresh clone). Watch the output for any `WARNING: failed to
fetch <accession>` lines -- if you see one, check your network can reach
`rest.uniprot.org` and `eutils.ncbi.nlm.nih.gov`, then re-run.

**Sanity check** -- confirm the right number of sequences landed:

```bash
grep -c '^>' data/fad2_seed_core.fasta      # should print 4
grep -c '^>' data/fad2_seed_legumes.fasta   # should print 11 (only exists if you ran --legumes)
```

## 3. Run the pipeline

**Core set only:**

```bash
bash scripts/run_pipeline.sh
```

**Both tiers in one go** (requires having run `fetch_data.sh --legumes`
in step 2):

```bash
bash scripts/run_pipeline.sh --both
```

Each run prints a numbered `[1/6] ... [6/6]` progress log: MAFFT
alignment, FastTree, rerooting FastTree's tree on the Arabidopsis
outgroup, IQ-TREE (ModelFinder + 1000 ultrafast bootstrap, also
rerooted), pairwise identity, and figure rendering. The core set
finishes in seconds; the legume expansion takes a minute or two, mostly
in IQ-TREE's bootstrap step.

Outputs land in `results/` (core set) and, with `--both`, also
`results/legume_expansion/`.

## 4. Check the results

**The fast way -- one command, pass/fail:**

```bash
bash scripts/verify_results.sh results 4
bash scripts/verify_results.sh results/legume_expansion 11
```

Each prints `OK`/`FAIL` for every check (file presence, sequence count,
genuine outgroup rooting, ModelFinder ran, pairwise-identity matrix looks
sane) and exits non-zero if anything's wrong -- this is the version to
use if you're scripting this or handing it to another agent, rather than
reading output by eye.

**The manual way, if you want to look yourself:**

**Did it actually finish?** Each stage should have produced these files
(shown for `results/`; same names apply under `results/legume_expansion/`):

```bash
ls -1 results/
```

Expect to see: `fad2_aligned.fasta`, `fad2_tree.nwk`,
`fad2_tree_rooted.nwk` + `.png`, `fad2_tree_iqtree.nwk` + `.png`,
`fad2_iqtree.iqtree`, `pairwise_identity.tsv`, and the `*.log` files
(`mafft.log`, `fasttree.log`, `iqtree.log`).

**Look at the trees:**

```bash
# on a machine with a GUI / file browser, just open these:
results/fad2_tree_rooted.png       # FastTree, rooted
results/fad2_tree_iqtree.png       # IQ-TREE, rooted
```

**Confirm the rooting is real**, not cosmetic -- a genuinely rooted tree
has exactly one comma at the top level of its Newick string (two children
at the root), with the outgroup (`P46313`) as one of them:

```bash
python3 -c "
with open('results/fad2_tree_iqtree.nwk') as f:
    s = f.read().strip().rstrip(';')
depth = 0; commas = 0
for c in s:
    if c == '(': depth += 1
    elif c == ')': depth -= 1
    elif c == ',' and depth == 1: commas += 1
print('top-level commas:', commas, '(expect 1)')
print('outgroup present:', 'P46313' in s)
"
```

**Look at the numbers:**

```bash
cat results/pairwise_identity.tsv
grep -i "Best-fit model" results/fad2_iqtree.iqtree
```

## 5. Re-running after a change

If you edit `scripts/fetch_data.sh` (e.g. to add a species) or just want
a completely clean run, clear out the generated files first -- nothing
in `data/` or `results/` is hand-edited, so it's always safe to delete:

```bash
rm -rf data/*.fasta results/
bash scripts/fetch_data.sh --legumes
bash scripts/run_pipeline.sh --both
```

## 6. Committing changes

`data/` is gitignored (raw fetched sequences aren't committed -- re-run
`fetch_data.sh` to regenerate them). `results/` **is** committed as the
worked example. After a run you're happy with:

```bash
git add results/ README.md environment.yml scripts/
git status   # check nothing unexpected (e.g. stray data/ files) is staged
git commit -m "your message here"
git push
```

## Troubleshooting quick reference

| Symptom | Likely cause | Fix |
|---|---|---|
| `mafft: command not found` | `fad2-legume-phylo` env not activated | `conda activate fad2-legume-phylo` |
| `iqtree2: command not found` | Wrong binary name assumed | The binary is `iqtree`, not `iqtree2` -- the scripts here already call it correctly |
| `curl` fetch warnings in step 2 | Network can't reach UniProt/NCBI | Check connectivity, re-run `fetch_data.sh` |
| IQ-TREE errors about bootstrap and taxon count | Fewer than 4 sequences in the input | Shouldn't happen with the provided seed files (always 4+); if you're supplying a custom FASTA, add more sequences |
| A tree's `.png` shows the outgroup nested *inside* the legume clade, not split off alone | Rerooting didn't run or found the wrong leaf | Re-check `fad2_tree.nwk` and `fad2_tree_iqtree.nwk` (raw, pre-reroot) headers actually contain `P46313`; re-run `run_pipeline.sh` |

See `README.md` for the full story behind each of these (including why
IQ-TREE's `-o` flag alone isn't enough to root its output).
