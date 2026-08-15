# DevKidCC v0.5.1 — validation of the container and the two entry points

Run 2026-08-15. Everything below is reproducible from the scripts in this
directory; §6 gives the commands.

## 1. What is being checked

DevKidCC v0.5.1 can be reached four ways, and this release is the first in which
the input file decides which stack reads it:

|  | `.h5ad` (AnnData) | `.rds` (Seurat) |
|---|---|---|
| **command line** | `/opt/dkcc` → `/opt/run_dkcc.py` | `/opt/dkcc` → `/opt/run_dkcc.R` |
| **interactive** | `classify_h5ad.py` (Python API) | `classify_seurat.R` (R API) |

All four end in the same `DevKidCC::DKCC()` call — the Python side shells out to
R over a CSV handoff — so the question worth answering is whether they arrive at
the same answer, and whether the container arrives at the same answer as the
machine the package was developed on.

Three claims are tested:

1. **All four paths classify.** Not "one path works and the others exist".
2. **The two routes agree with each other**, on the same cells.
3. **The container agrees with the host**, on the same cells.

## 2. Test data

Two 600-cell subsets, one organoid and one fetal kidney, so the comparison is
not resting on a single tissue:

| Name | Source | Cells | Genes |
|---|---|---|---|
| `organoid_howden` | Howden 2019 kidney organoids (`Howden_2019_Organoids_qc.h5ad`, 5,365 cells) | 600 | 22,073 |
| `fetal_menon` | Menon 2018 human fetal kidney (`Menon_2018_HFK_qc.h5ad`) | 600 | 26,469 |

Drawn with `numpy.random.default_rng(42).choice(n, 600, replace=False)`, sorted —
`make_test_data.py` reproduces exactly these barcodes. Both are **raw counts**
(the `_qc` files, not the `_processed` ones, which are log-normalised and would
silently produce nonsense).

Each subset is written twice from one source: `.h5ad` for the Python route, and
a MatrixMarket trio that `build_rds_from_mtx.R` turns into a Seurat `.rds`. The
two files therefore hold the same cells and the same counts, which is what makes
the route comparison a comparison rather than two unrelated runs.

600 cells is small on purpose: it classifies in about a minute and fits in a
couple of GB, while still being large enough that a real disagreement between
the routes would show.

## 3. Environments

| | Host | Container |
|---|---|---|
| R | 4.5.3 | 4.4.3 |
| Seurat | 5.4.0 | 5.5.1 |
| SeuratObject | 5.3.0 | 5.4.0 |
| scPred | 1.9.2 | 1.9.2 |
| DevKidCC | 0.5.1 | 0.5.1 |
| Python | 3.12.13 | 3.12.13 |
| anndata | 0.12.10 | 0.13.2 |
| numpy | 2.2.6 | 2.4.6 |
| scipy | 1.16.3 | 1.18.0 |

Image: `ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1`, converted to a
3.9 GB SIF. The container runs were made under

```
singularity exec --contain --cleanenv --no-home ...
```

so no host `$HOME`, no host `/tmp`, and no host environment variables cross the
boundary — an external user's conditions, not a developer's. R and Python
versions differ across the two columns; that is the point of comparing them.

## 4. Results

### 4.1 All four paths run

Ten runs, in the container, no failures:

| Dataset | Path | Result | Time |
|---|---|---|---|
| organoid_howden | `.h5ad` CLI | PASS | 91 s |
| organoid_howden | `.h5ad` interactive | PASS | 90 s |
| organoid_howden | `.rds` CLI | PASS | 77 s |
| organoid_howden | `.rds` interactive | PASS | 62 s |
| organoid_howden | route agreement | PASS | 2 s |
| fetal_menon | `.h5ad` CLI | PASS | 60 s |
| fetal_menon | `.h5ad` interactive | PASS | 62 s |
| fetal_menon | `.rds` CLI | PASS | 52 s |
| fetal_menon | `.rds` interactive | PASS | 57 s |
| fetal_menon | route agreement | PASS | 2 s |

The same matrix was run on the host beforehand, with the same outcome.

A separate check: the Python route returns the caller's object at its **input
width** — 600 × 22,073 in, 600 × 22,073 out. The reference-gene projection
happens inside the handoff to R and does not reach back into the user's AnnData.
This was a real defect in the pre-release wrapper, where callers silently got
back roughly 10,000 genes.

### 4.2 The two routes agree

Python (`.h5ad`) against R (`.rds`), same 600 cells:

| Dataset | Column | Agreement | Cohen's κ |
|---|---|---|---|
| organoid_howden | LineageID | 99.00% | 0.9845 |
| organoid_howden | DKCC | 98.83% | 0.9850 |
| fetal_menon | LineageID | 96.33% | 0.9408 |
| fetal_menon | DKCC | 95.00% | 0.9421 |

**The residual disagreement is deterministic, not noise** — the host and the
container produced these figures to four decimal places independently, so a
re-run does not move them.

It is also structured. Every disagreement sits between adjacent states rather
than across the lineage tree: organoid is NPC ↔ NPC-like (3 cells) and
NPC-like ↔ Stroma (3); fetal is NPC-like ↔ Nephron (12 of the 22). Cells on the
boundary between a progenitor and what it becomes are exactly where a
probability threshold of 0.7 is doing the most work.

Where it does *not* come from, both checked in the run logs:

- **Not gene coverage.** The per-model reference overlap is identical between the
  routes — 6,047 / 6,027 / 5,735 / 5,910 / 5,880 / 5,817 features across the six
  scPred models on organoid, both ways. The Python side's reference-gene
  pre-filter (which exists to keep the CSV handoff from ballooning) is lossless
  with respect to what scPred sees.
- **Not normalisation depth.** The Python route carries full-matrix library sizes
  across in `dkcc_full_library_size` and normalises with those, so the pre-filter
  does not distort sequencing-depth correction.

What does differ is the object each route hands to the KNN smoothing step at the
end of `DKCC()`: 6,602 genes on the Python side against the full 22,073 on the R
side, and the R route additionally computes a UMAP that the Python route does
not. A neighbour graph built over different feature widths reassigns a handful of
borderline cells differently. That accounts for the size and the placement of the
disagreement, though it has not been isolated by ablation.

**For interpretation:** the routes are interchangeable at the level of composition
(§4.4) and for any cell not sitting on a decision boundary. Two runs of the same
sample should still go down the same route, and a figure comparing samples should
state which one.

### 4.3 The container agrees with the host

Same cells, same route (Python), host against container:

| Dataset | Column | Agreement | Cohen's κ |
|---|---|---|---|
| organoid_howden | LineageID | 100.00% | 1.0000 |
| organoid_howden | DKCC | 100.00% | 1.0000 |
| fetal_menon | LineageID | 100.00% | 1.0000 |
| fetal_menon | DKCC | 100.00% | 1.0000 |

Cell for cell identical, across R 4.4.3 vs 4.5.3 and Seurat 5.5.1 vs 5.4.0. The
summary CSVs are byte-identical too, and the R-route composition tables likewise
match. Near-identical was the pass condition; identical is what came back.

### 4.4 Composition

`organoid_howden`, 600 cells, 100% assigned:

| LineageID | n | % |
|---|---|---|
| Stroma | 265 | 44.17 |
| Nephron | 226 | 37.67 |
| NPC-like | 57 | 9.50 |
| NPC | 52 | 8.67 |

Top DKCC classes: CS 250 (41.67%), EN 68 (11.33%), NPC-like 57, EDT 55, EPT 54,
NPC 52, EPod 25, PEC 15.

`fetal_menon`, 600 cells, 599 assigned (99.83%):

| LineageID | n | % |
|---|---|---|
| Nephron | 340 | 56.67 |
| NPC-like | 115 | 19.17 |
| Stroma | 73 | 12.17 |
| UrEp | 46 | 7.67 |
| Endo | 25 | 4.17 |
| unassigned | 1 | 0.17 |

Top DKCC classes: EDT 155 (25.83%), NPC-like 115, EN 75, EPT 37, MS 35, EPod 27,
Endo 25, UOS 24, Pod 21.

The two profiles separate the way the biology says they should. The organoid is
stroma-dominant with no endothelium and no ureteric epithelium — both are
well-documented absences in this protocol. The fetal kidney is nephron-dominant
and carries both the UrEp (46 cells across UOS/UIS/UTip) and Endo (25) that the
organoid lacks, plus mature podocytes (21 Pod against the organoid's 3). This is
a sanity check on the classifier, not a finding.

## 5. Two traps, documented because both cost time

**Seurat v5 objects are not normalised for you.** `DKCC()` calls `NormalizeData`
only on its v4 branch; for an Assay5 object it joins layers and reads whatever is
in the `data` layer. Hand it an unnormalised v5 object and it classifies raw
counts without complaint. `classify_seurat.R` checks
`length(Layers(seu, search = "data")) > 0` and normalises when there is nothing
there.

**Log-normalised input fails silently.** The scPred models expect counts.
`classify_h5ad.py` refuses input whose maximum value looks log-scaled rather than
letting it through — the failure mode otherwise is plausible-looking labels that
are wrong.

## 6. Reproducing this

Inside the container, nothing else installed:

```bash
singularity pull dkcc.sif docker://ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1

# build the inputs (needs the source h5ads; skip if you have the test data)
python /opt/examples/make_test_data.py --out ./data
Rscript /opt/examples/build_rds_from_mtx.R \
    data/organoid_howden.mtx data/organoid_howden.genes.txt \
    data/organoid_howden.cells.txt data/organoid_howden.rds

# run all four paths on both datasets and compare the routes
singularity exec --contain --cleanenv --no-home \
    --bind $PWD/data:/data dkcc.sif \
    bash /opt/examples/run_examples.sh /data organoid_howden fetal_menon
```

`run_examples.sh` writes per-run outputs and summary CSVs to `<data>/results`
and prints the PASS/FAIL table in §4.1. `compare_routes.py` produces §4.2.

Any raw-count h5ad with HGNC gene symbols works in place of the two here.

## 7. Scope

- Both test sets are 600 cells. Timings do not extrapolate to a full sample;
  Howden's 5,365 cells take about 90 s on the host.
- The route comparison is one pair of samples per tissue. It establishes that the
  routes agree closely, not the exact rate at which they disagree in general.
- HPC (SLURM/Apptainer) is covered separately in `HPC_TESTING.md` and is not
  included above.
