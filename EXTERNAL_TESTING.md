# Testing the DevKidCC container — protocol for external testers

DevKidCC classifies cells in single-cell RNA-seq data of the developing human
kidney and of kidney organoids. This release ships as a container carrying both
an R and a Python stack, so it runs the same way on a laptop and on a cluster
without you installing anything.

**We are asking you to try to break it.** The container has been validated on the
machine it was developed on; what that cannot tell us is whether it works on
*your* machine, with your data, your scheduler and your filesystem. A failure
here is a useful result, not a wasted afternoon — please report what happens
either way.

Budget about **90 minutes**, most of it waiting on a 4 GB download.

<!-- ============================================================
     NOT READY TO SEND — read this first (delete before sending)

     The image currently on GHCR as :containerise-v0.5.1 was built from
     commit 40762f4, which is #11 of the 16 on this branch. Five commits
     are missing from it, and two of them break this protocol:

       453cfc9  numba/matplotlib cache dirs  -> `import devkidcc` fails
                under --contain --no-home. Worked around in S3 with two
                --env lines.
       1afc43f  COPY examples/ /opt/examples -> NOT IN THE IMAGE AT ALL.
                S5, S6 and Appendix A all invoke /opt/examples/... and
                will fail with "no such file or directory".

     There is no workaround for the second one worth shipping; the
     examples have to be in the image. Required before sending:

       1. push 197da46
       2. confirm CI rebuilds :containerise-v0.5.1 from the branch HEAD
       3. re-verify:
            singularity exec --contain --cleanenv --no-home dkcc.sif \
                bash -c "ls /opt/examples && python -c 'import devkidcc'"
       4. delete the two --env lines in S3 and the note explaining them
       5. delete this box
     ============================================================ -->

**Contents**

1. [What you need](#1-what-you-need)
2. [Get the image](#2-get-the-image)
3. [Set up the run command](#3-set-up-the-run-command)
4. [Self-test — no data required](#4-self-test--no-data-required)
5. [Classify the reference data](#5-classify-the-reference-data)
6. [Classify your own data](#6-classify-your-own-data)
7. [Run it on a cluster](#7-run-it-on-a-cluster)
8. [Things that are known to bite](#8-things-that-are-known-to-bite)
9. [What to send back](#9-what-to-send-back)

---

## 1. What you need

### Both environments

| | Requirement |
|---|---|
| Disk | ~10 GB free (the image is 4 GB, and building a SIF needs room for two copies) |
| RAM | 8 GB is enough for the test data; budget roughly 1 GB per 1,000 cells for your own |
| Network | Only to pull the image. Classification itself needs no network. |

No R, no Python, no DevKidCC installation. Everything is inside the container.

### On a laptop

Either **Docker Desktop** (macOS, Windows, Linux) or **Apptainer/Singularity**
(Linux). Docker is the easier route on macOS and Windows.

```bash
docker --version          # any recent version
# or
apptainer --version       # or: singularity --version
```

> **Windows:** run everything below inside **WSL2** or Git Bash, not PowerShell —
> the commands assume a POSIX shell. Docker Desktop with the WSL2 backend works
> well.
>
> **Apple Silicon (M1–M4):** the image is `linux/amd64` only. Docker will run it
> under emulation with a warning; expect it to be roughly 2–4× slower, but it
> should still produce identical labels. **Please tell us if it does not** — that
> is exactly the kind of thing we cannot test here.

### On a server or cluster

**Apptainer** or **Singularity**. You will not have Docker, and you do not need it.

```bash
module load apptainer     # or: module load singularity — site-dependent
apptainer --version
```

Also worth checking now, because §2 may need it:

```bash
which skopeo || module spider skopeo
```

---

## 2. Get the image

The image is public on the GitHub Container Registry. **No login, no token.**

> ### Use the explicit tag below — do not use `:latest`
>
> `:latest` currently points at an older, known-broken build. It will be corrected
> when this release merges, but until then it is a trap, and it is the tag most
> tools reach for by default. Every command in this document names the tag
> explicitly. If you find yourself typing `:latest`, stop.

<!-- TODO before sending: if v0.5.1 has been tagged and merged to main, replace
     :containerise-v0.5.1 with :v0.5.1 throughout and delete the box above. -->

### Laptop, with Docker

```bash
docker pull ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1
docker images ghcr.io/kidneyregeneration/dkcc        # expect ~4 GB
```

### Laptop or server, with Apptainer/Singularity

```bash
mkdir -p ~/dkcc-test && cd ~/dkcc-test

singularity pull dkcc-v051.sif \
    docker://ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1

ls -lh dkcc-v051.sif  # expect ~3.9 GB
```

#### If that pull dies partway through

You may see `stream ID 15; PROTOCOL_ERROR` or an unexpected EOF. This is an
HTTP/2 transport fault between Singularity and GHCR — not permissions, not disk.
**Retrying will not help**: Singularity restarts the whole download each time
rather than resuming, so a 4 GB image never reaches the end on a flaky link. We
hit this on every single attempt from one of our sites.

Use `skopeo`, which reconnects *inside* a layer and carries on from the byte it
reached:

```bash
export GODEBUG=http2client=0     # drop to HTTP/1.1

skopeo copy --retry-times 10 \
    docker://ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1 \
    oci:$PWD/dkcc-oci:v051

# Build the SIF from the local copy — no network involved.
# Note there is no :tag on the path here. Singularity treats everything after
# `oci:` as a filesystem path, so `oci:$PWD/dkcc-oci:v051` would look for a
# directory literally named "dkcc-oci:v051".
singularity build dkcc-v051.sif oci:$PWD/dkcc-oci

rm -rf dkcc-oci                  # the OCI directory is only an intermediate
```

Log lines like `Reading blob body ... failed (unexpected EOF), reconnecting after
1614807040 bytes` are the resume working. That is success, not an error.

**Please report which of these two routes you needed.** If plain `singularity
pull` works at your site, that is worth knowing too.

### If you already have the image locally

Some sites mirror images into an internal registry, or hand you a `.sif` on
shared storage rather than letting compute nodes reach the internet. If someone
has given you the file directly, skip this whole section — nothing above is
required.

The rest of this document never names the image again. It refers to it only as
`$SIF`, which you set once in §3, so pointing that at a local file is the only
change:

```bash
SIF=/path/to/the/file/you/were/given.sif
```

Everything from §3 onward then runs unchanged. If you do this, please say so in
your report and tell us where the file came from — we need to know which build
you actually exercised, and a `.sif` on shared storage carries no tag to check.
You can read the build it came from with:

```bash
singularity inspect $SIF | grep -E 'image.version|image.revision|image.created'
```

---

## 3. Set up the run command

Every command from here on runs inside the container. Rather than repeat a long
invocation, define it once.

**Apptainer/Singularity** — from `~/dkcc-test`:

```bash
cd ~/dkcc-test

SIF=$PWD/dkcc-v051.sif
DKCC="singularity exec --contain --cleanenv --no-home \
      --env NUMBA_CACHE_DIR=/tmp/numba-cache \
      --env MPLCONFIGDIR=/tmp/matplotlib"
```

**Docker** — the equivalent:

```bash
cd ~/dkcc-test

SIF=ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1
DKCC="docker run --rm \
      -e NUMBA_CACHE_DIR=/tmp/numba-cache \
      -e MPLCONFIGDIR=/tmp/matplotlib"
```

Then `$DKCC $SIF <command>` runs `<command>` inside the container, and where a
command needs your files visible inside it, a bind goes between the two:

```bash
$DKCC --bind $PWD/data:/data $SIF <command>      # Singularity
$DKCC -v "$PWD/data:/data"    $SIF <command>     # Docker
```

The order matters — the bind is an option to `singularity exec` / `docker run`, so
it must come **before** the image. Re-define both variables if you open a new shell.

> **Why `--contain --cleanenv --no-home`?** They stop your home directory, your
> environment variables and your host `/tmp` from leaking in. That makes this a
> test of the *image* rather than of your machine, and it is the configuration
> most likely to expose a packaging bug. Please keep them on.

> **Why the two `NUMBA_CACHE_DIR` / `MPLCONFIGDIR` variables?** Interim
> workaround. numba compiles some of scanpy's routines with caching on; in a
> read-only container with no home directory it has nowhere to put the cache and
> raises rather than skipping it, so `import devkidcc` fails outright. The fix is
> committed and will be baked into the image at the next rebuild — until then,
> these two lines supply it from outside. **If you have time, please also try
> once without them and tell us what you see**; that tells us whether the rebuild
> has landed.

<!-- TODO before sending: once CI has rebuilt the image from a commit that
     includes the NUMBA_CACHE_DIR/MPLCONFIGDIR ENV, delete the two --env lines
     from both blocks above and the note explaining them. Verify first with:
       singularity exec --contain --cleanenv --no-home dkcc-v051.sif \
           python -c "import devkidcc; print('ok')" -->

---

## 4. Self-test — no data required

The image carries its own test. It synthesises a small counts matrix and pushes
it through all three ways a user can reach the classifier, checking that each one
returns labels for every cell.

```bash
$DKCC $SIF python /opt/smoke_test.py ; echo "exit=$?"
```

Takes 3–6 minutes. It checks three things:

1. the Python API, `devkidcc.classify_kidney_cells()`
2. an `.h5ad` through `/opt/dkcc`, which must route to Python
3. an `.rds` through `/opt/dkcc`, which must route to R

The labels it produces are meaningless — the input is synthetic noise — and the
test deliberately does not check them. What it checks is that each path completes
and returns `DKCC` and `LineageID` for every cell. Seeing every cell come back
`unassigned` here is the correct result, not a failure.

**Expect**, as the last two lines:

```
PASS: the Python API, h5ad routing and .rds routing all classified the input.
exit=0
```

### Confirm the versions

Record these for your report:

```bash
$DKCC $SIF python -c "from importlib.metadata import version; import sys; \
print('python  ', sys.version.split()[0]); \
[print(f'{p:8}', version(p)) for p in ('anndata','numpy','scipy','scanpy')]"

$DKCC $SIF Rscript -e 'cat(sprintf("%-8s %s\n","R",getRversion())); \
for (p in c("Seurat","scPred","DevKidCC")) cat(sprintf("%-8s %s\n", p, packageVersion(p)))'
```

**Expect exactly:**

```
python   3.12.13
anndata  0.13.2
numpy    2.4.6
scipy    1.18.0
scanpy   1.12.3
R        4.4.3
Seurat   5.5.1
scPred   1.9.2
DevKidCC 0.5.1
```

Any difference means you have pulled a different image than we think you have.
Stop and tell us before going further.

### Check the routing

The headline change in this release is that **the input file decides which stack
reads it**, rather than the caller choosing:

| Input | Goes to |
|---|---|
| `.h5ad` | `/opt/run_dkcc.py` — Python, via anndata |
| `.rds` `.RData` `.h5seurat` `.h5` | `/opt/run_dkcc.R` — R, via Seurat |

Both end in the same `DevKidCC::DKCC()` call. Confirm the dispatcher agrees:

```bash
$DKCC $SIF /opt/dkcc --input nonexistent.h5ad --output /tmp/o.h5ad 2>&1 | head -1
$DKCC $SIF /opt/dkcc --input nonexistent.rds  --output /tmp/o.rds  2>&1 | head -1
$DKCC $SIF /opt/dkcc --input nonexistent.csv  --output /tmp/o.csv  2>&1 | head -2
```

**Expect:**

```
Routing to the Python entry point (.h5ad is AnnData's format)
Routing to the R entry point (.rds is read by Seurat)
ERROR: unsupported input extension '.csv'
       Supported: .h5ad (Python) .rds .RData .h5seurat .h5 (R)
```

The first two then fail on the missing file. That is fine — the routing line is
the whole point.

---

## 5. Classify the reference data

This is the strongest check available to you, because it has a known right answer.

Two 600-cell subsets, one kidney organoid and one human fetal kidney:

| Name | Source | Accession | Cells | Genes |
|---|---|---|---|---|
| `organoid_howden` | Howden et al. 2019, iPSC-derived kidney organoids | GSE118184 | 600 | 22,073 |
| `fetal_menon` | Menon et al. 2018, human fetal kidney | GSE109205 | 600 | 26,469 |

Each is provided twice — as `.h5ad` and as `.rds`, built from the *same cells* —
so the two entry points can be compared directly.

Download (~32 MB):

```bash
mkdir -p ~/dkcc-test/data && cd ~/dkcc-test/data
curl -L -O <<TEST_DATA_URL>>/organoid_howden.h5ad
curl -L -O <<TEST_DATA_URL>>/organoid_howden.rds
curl -L -O <<TEST_DATA_URL>>/fetal_menon.h5ad
curl -L -O <<TEST_DATA_URL>>/fetal_menon.rds
cd ~/dkcc-test
```

<!-- TODO before sending: replace <<TEST_DATA_URL>> with the real location.
     32 MB is too large to email comfortably; a GitHub release asset or a Zenodo
     record is the right home, and a release asset gives a stable curl-able URL. -->

### Run everything

One command runs the full matrix — both datasets, both formats, CLI and
interactive, plus a cell-by-cell comparison of the two routes.

**Apptainer/Singularity:**

```bash
$DKCC --bind $PWD/data:/data $SIF \
    bash /opt/examples/run_examples.sh /data organoid_howden fetal_menon
```

**Docker** — same shape, `-v` instead of `--bind`:

```bash
$DKCC -v "$PWD/data:/data" $SIF \
    bash /opt/examples/run_examples.sh /data organoid_howden fetal_menon
```

Substitute `-v "$PWD/…:/data"` for `--bind $PWD/…:/data` in every command from
here on if you are on Docker; nothing else changes.

**Expect:** ten runs, all `PASS`, in roughly 10–12 minutes. Outputs land in
`data/results/`. Individual runs took 52–91 s each on our reference machine; much
slower is worth noting but is not a failure.

```
  PASS  91s  organoid_howden  .h5ad  CLI (/opt/dkcc)
  PASS  90s  organoid_howden  .h5ad  interactive (classify_h5ad.py)
  PASS  77s  organoid_howden  .rds   CLI (/opt/dkcc)
  PASS  62s  organoid_howden  .rds   interactive (classify_seurat.R)
  PASS   2s  organoid_howden  route agreement (Python h5ad vs R rds)
  ... and the same five for fetal_menon

  10 passed, 0 failed
```

### Check the answers

The labels are deterministic, so yours should match ours.

```bash
cd ~/dkcc-test/data/results
column -s, -t organoid_howden_py_summary.csv
column -s, -t fetal_menon_py_summary.csv
```

**`organoid_howden`** — 600 cells, 100% assigned:

| LineageID | n | % |
|---|---|---|
| Stroma | 265 | 44.17 |
| Nephron | 226 | 37.67 |
| NPC-like | 57 | 9.5 |
| NPC | 52 | 8.67 |

Top `DKCC` classes: CS 250, EN 68, NPC-like 57, EDT 55, EPT 54, NPC 52, EPod 25,
PEC 15.

**`fetal_menon`** — 600 cells, 599 assigned:

| LineageID | n | % |
|---|---|---|
| Nephron | 340 | 56.67 |
| NPC-like | 115 | 19.17 |
| Stroma | 73 | 12.17 |
| UrEp | 46 | 7.67 |
| Endo | 25 | 4.17 |
| unassigned | 1 | 0.17 |

Top `DKCC` classes: EDT 155, NPC-like 115, EN 75, EPT 37, MS 35, EPod 27, Endo 25,
UOS 24, Pod 21.

These should match **exactly**. We have reproduced them across R 4.4.3 vs 4.5.3
and Seurat 5.5.1 vs 5.4.0 with cell-for-cell identity (Cohen's κ = 1.0000), so a
difference of even a few cells is a real finding — please send us the CSVs.

The two profiles separate the way the biology says they should, which is a sanity
check rather than a result: the organoid is stroma-dominant with no endothelium
and no ureteric epithelium, both well-documented absences in that protocol; the
fetal kidney is nephron-dominant and carries both.

### The two routes are close but not identical

```bash
column -s, -t organoid_howden_agreement.csv
column -s, -t fetal_menon_agreement.csv
```

**Expect:**

| Dataset | Column | Agreement | Cohen's κ |
|---|---|---|---|
| organoid_howden | LineageID | 0.99 | 0.9845 |
| organoid_howden | DKCC | 0.9883 | 0.985 |
| fetal_menon | LineageID | 0.9633 | 0.9408 |
| fetal_menon | DKCC | 0.95 | 0.9421 |

**This gap is expected, understood and deterministic** — it is not noise, and
re-running will not move it. The two routes hand different-width objects to the
KNN smoothing step at the end of classification (6,602 genes on the Python side
against the full 22,073 on the R side), and the R route additionally computes a
UMAP the Python route does not. The disagreements sit almost entirely between
adjacent progenitor states — NPC ↔ NPC-like, NPC-like ↔ Nephron — which is
exactly where a probability threshold is doing the most work.

**What this means in practice:** the routes are interchangeable for composition,
and for any cell not sitting on a decision boundary. But run all samples in a
comparison down the *same* route, and state which one in your methods.

---

## 6. Classify your own data

This is the part we most want tested, because it is the part we cannot simulate.

### Requirements for your input

- **Raw counts.** Not log-normalised, not scaled, not SCTransform residuals. The
  models were fit against counts and normalisation happens inside.
- **HGNC gene symbols** as feature names — `PAX2`, `SIX2`, `GATA3`. Ensembl IDs
  (`ENSG…`) will not match. Neither will Cellranger multi-genome names
  (`GRCh38_PAX2`); strip the prefix first.
- **Human kidney or kidney organoid.** It will happily label other tissues, and
  the labels will be meaningless.
- Any size, but start small — a few thousand cells — for a first run.

### Preferred: the worked example

```bash
mkdir -p ~/dkcc-test/mydata
cp /path/to/your_sample.h5ad ~/dkcc-test/mydata/
cd ~/dkcc-test

$DKCC --bind $PWD/mydata:/data $SIF \
    python /opt/examples/classify_h5ad.py \
        --input   /data/your_sample.h5ad \
        --output  /data/your_sample_DKCC.h5ad \
        --summary /data/your_sample_summary.csv
```

Use this rather than the bare CLI for a first run: it prints diagnostics as it
goes, **refuses log-normalised input rather than silently producing wrong
labels**, and writes a composition summary alongside the output.

A Seurat object goes through the R equivalent, which also normalises v5 objects
for you (see §8):

```bash
$DKCC --bind $PWD/mydata:/data $SIF \
    Rscript /opt/examples/classify_seurat.R \
        --input   /data/your_sample.rds \
        --output  /data/your_sample_DKCC.rds \
        --summary /data/your_sample_summary.csv
```

### The bare CLI

`/opt/dkcc` is the minimal machine-facing path — one command for both formats,
routed on the extension:

```bash
$DKCC --bind $PWD/mydata:/data $SIF \
    /opt/dkcc --input /data/your_sample.h5ad --output /data/your_sample_DKCC.h5ad

$DKCC --bind $PWD/mydata:/data $SIF \
    /opt/dkcc --input /data/your_sample.rds  --output /data/your_sample_DKCC.rds
```

Options beyond `--input`/`--output` pass straight through to whichever entry point
runs: `--threshold` (assignment probability, default 0.7), `--max-iter`,
`--knn-iter` (0 disables KNN smoothing).

### Check what came back

```bash
$DKCC --bind $PWD/mydata:/data $SIF python -c "
import anndata as ad
a = ad.read_h5ad('/data/your_sample_DKCC.h5ad')
print(a)
print(a.obs['LineageID'].value_counts())
print(a.obs['DKCC'].value_counts())
"
```

Three things to verify:

1. **`DKCC` and `LineageID` are present** in `.obs`, with a value for every cell.
2. **`n_vars` is unchanged** from your input. The classifier projects onto its
   reference genes internally, but that must not reach back into your object.
   Getting back ~10,000 genes when you supplied 30,000 is a bug we fixed in this
   release, and a regression we want caught.
3. **The labels are biologically plausible** for what you gave it. You know your
   data; we do not. This is the judgement we cannot make from here.

### Batch mode

Point it at a folder and it processes every input once, routing each file on its
own extension — a folder may hold both formats:

```bash
$DKCC --bind $PWD/mydata:/data $SIF bash /opt/run_dkcc_batch.sh -i /data
```

Worth watching the file list it prints before it starts: it should enumerate your
inputs and *not* pick up its own `_DKCC.` outputs as new work.

---

## 7. Run it on a cluster

Skip this if you are only testing a laptop.

The repository ships a launcher that submits to SLURM. Clone it next to your SIF:

```bash
cd ~/dkcc-test
git clone -b containerise-v0.5.1 \
    https://github.com/KidneyRegeneration/DevKidCC dkcc-repo
cd dkcc-repo
```

> **Pass `--container <path to your .sif>`.** Without it the launcher pulls
> `:latest`, which is the broken tag from §2. This default will be corrected; for
> now, always pass the file.

```bash
./run_dkcc.sh \
    --mode slurm \
    --container ~/dkcc-test/dkcc-v051.sif \
    --partition <your partition> \
    --time 02:00:00 \
    --mem 32G \
    --cpus 4 \
    --mounts /scratch \
    --logs ~/dkcc-test \
    ~/dkcc-test/data/organoid_howden.h5ad
```

`--mounts` takes a comma-separated list of host paths to bind in — set it to
whichever filesystems your data actually lives on. The defaults in the script
(`/group`, partition `prod_med`) are our institution's conventions and are almost
certainly wrong for you. **Please tell us what you had to change**, since making
this portable is one of the things we are testing.

Watch it, then check the output:

```bash
squeue -u $USER
tail -f ~/dkcc-test/dkcc_organoid_howden_*.log
ls -la ~/dkcc-test/data/organoid_howden_DKCC.h5ad
```

**Expect** the same labels as your local run in §5. 5,365 cells took ~90 s on our
reference machine, so a 600-cell job on a 2-hour partition is generous.

If your site uses PBS, LSF or SGE rather than SLURM, please say so — the launcher
only speaks SLURM today, and we would like to know what else to support. Calling
`singularity exec` from your own submission script works regardless of scheduler,
and is a fine substitute for this section.

---

## 8. Things that are known to bite

**Log-normalised input fails silently through `/opt/dkcc`.** The classifier will
run and return confident, plausible, wrong labels. `classify_h5ad.py` refuses such
input; the bare CLI path does not yet check. If your `X.max()` is under ~50 and
the values are not integers, you have log-normalised data — use the raw counts
file, not the processed one.

**Seurat v5 objects are not normalised for you.** `DKCC()` normalises on its v4
branch only; hand it an Assay5 object with nothing in the `data` layer and it
classifies raw counts without complaint. `classify_seurat.R` checks for this and
normalises when needed — another reason to prefer it over the bare CLI.

**`:latest` is the wrong image.** Covered in §2, repeated here because it is the
single most likely way to waste an afternoon.

**`import devkidcc` fails without a writable cache directory.** Covered in §3.
If you see `RuntimeError: cannot cache function 'agg_sum_csr-parallel': no locator
available`, you have dropped the two `--env` lines too early.

**Your home directory can leak in.** Without `--contain --no-home`, Singularity
bind-mounts `$HOME`, which can put your own `~/.local/lib/python3.x/site-packages`
ahead of the container's. The image sets `PYTHONNOUSERSITE=1` to defend against
this, but the flags are still the safer habit.

**HDF5 version warning.** You will see `h5py is running against HDF5 2.2.0 when it
was built against 2.1.0`. It is harmless and on our list to silence. Please ignore
it — but do mention any *other* warning you see.

---

## 9. What to send back

Please report what happened either way. A clean run is as useful to us as a
failure, because it tells us the image travels.

Send to **<<CONTACT_EMAIL>>**:

```
ENVIRONMENT
  Machine        laptop / server / cluster
  OS + version
  CPU            e.g. Apple M2, Intel Xeon
  RAM
  Container tool docker X.Y / apptainer X.Y / singularity X.Y
  Scheduler      none / SLURM / other

§2  IMAGE PULL
  Route that worked   singularity pull / skopeo / docker pull
  Time taken
  Errors seen (paste them, even if a retry fixed it)

§4  SELF-TEST
  smoke_test.py           PASS / FAIL
  Version block           matched exactly? (paste if not)
  Routing check           PASS / FAIL
  Did it work without the two --env lines?   yes / no / didn't try

§5  REFERENCE DATA
  run_examples.sh         N passed, N failed
  Composition matched?    yes / no  (attach the *_summary.csv files)
  Route agreement         (attach the *_agreement.csv files)
  Wall time

§6  YOUR OWN DATA
  Dataset       cells x genes, tissue, how it was generated
  Format        h5ad / rds
  Result        PASS / FAIL
  DKCC + LineageID present for every cell?   yes / no
  n_vars unchanged?                          yes / no
  Do the labels look biologically sensible to you?
  Wall time, and peak memory if you have it

§7  CLUSTER (if attempted)
  Scheduler
  What you had to change in the run_dkcc.sh command
  Result

ANYTHING ELSE
  Confusing instructions, missing steps, unclear output, warnings we did
  not warn you about. Blunt is welcome.
```

For any failure, the most useful things to attach are **the full console output**
(not just the last line) and the exact command you ran.

---

## Appendix A — rebuilding the test data from scratch

Not required for this protocol, included so the reference numbers in §5 are
reproducible rather than asserted.

Both source datasets are public: GSE118184 (Howden) and GSE109205 (Menon). The
600-cell subsets are drawn with a fixed seed, `numpy.random.default_rng(42)
.choice(n, 600, replace=False)`, sorted. The same script will build a test pair
from any raw-count h5ad of your own:

```bash
$DKCC --bind $PWD/data:/data $SIF \
    python /opt/examples/make_test_data.py \
        --out /data \
        --source mysample=/data/my_raw_counts.h5ad

$DKCC --bind $PWD/data:/data $SIF \
    Rscript /opt/examples/build_rds_from_mtx.R \
        /data/mysample.mtx /data/mysample.genes.txt \
        /data/mysample.cells.txt /data/mysample.rds
```

One caveat: the cell counts in §5 are **after our QC filtering**, so a subset
drawn from the raw GEO matrices will not be the same 600 cells and will not give
the same numbers. The figures in §5 are reproducible from the distributed files,
which is why we distribute them.

The full validation record, including the environment table and the analysis of
the route difference, is `examples/VALIDATION.md` — inside the image at
`/opt/examples/VALIDATION.md`, and in the repository.

---

## Appendix B — citing this

If the container ends up in a manuscript, `CITATION.cff` in the repository has the
machine-readable form. The paper is:

> Wilson SB, Howden SE, Vanslambrouck JM, Dorison A, Alquicira-Hernandez J,
> Powell JE, Little MH. DevKidCC allows for robust classification and direct
> comparisons of kidney organoid datasets. *Genome Medicine* 14, 19 (2022).
> doi:10.1186/s13073-022-01023-z

Cite an immutable tag rather than `:latest`, so the version you ran is the version
a reader can pull.

<!-- TODO before publishing: once v0.5.1 is tagged and a Zenodo DOI minted, add
     the DOI here and in CITATION.cff. -->
