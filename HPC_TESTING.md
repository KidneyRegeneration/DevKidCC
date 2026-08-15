# HPC test checklist — DevKidCC container v0.5.1

Run this on the MCRI cluster to confirm the `containerise-v0.5.1` image works
under Apptainer/SLURM. Everything below is copy-paste; §7 says what to send back.

**Prerequisite:** the branch is pushed and the `docker-publish.yml` run for
`containerise-v0.5.1` is green, so
`ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1` exists.

**What is being tested that wasn't before:** the image is now self-contained.
Previous runs bind-mounted fixed copies of `run_dkcc.R` / `run_dkcc_batch.sh`
over the ones inside the SIF. Those fixes are baked in now, and `build_binds()`
no longer mounts any script. If a step below fails in a way that would once have
been papered over by a bind-mount, that is the finding — don't re-add the mount,
send the log.

---

## 1. Set up (login node)

Compute nodes have no outbound network, so the SIF must be pulled here.

```bash
export DKCC_TEST=/group/$USER/dkcc_v051_test      # adjust to your group path
mkdir -p "$DKCC_TEST" && cd "$DKCC_TEST"

module load singularity 2>/dev/null || module load apptainer
```

## 2. Pull the image (login node)

`run_dkcc.sh --container remote` pulls `:latest`, which is *not* the branch build.
Pull the tag explicitly and pass the file:

```bash
singularity pull dkcc-v051.sif \
    docker://ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1

ls -lh dkcc-v051.sif        # expect roughly 4 GB
```

If this 401s, the package is private — make it public in the GHCR package
settings, or `singularity remote login -u <github-user> docker://ghcr.io` with a
PAT that has `read:packages`.

## 3. Sanity-check the image (login node, no SLURM)

```bash
singularity exec dkcc-v051.sif Rscript -e '
  cat("DevKidCC:", as.character(packageVersion("DevKidCC")), "\n")
  cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
  cat("knn.iter present:", "knn.iter" %in% names(formals(DevKidCC::DKCC)), "\n")
  cat("sceasy  :", requireNamespace("sceasy", quietly=TRUE), "\n")
  cat("reticulate:", requireNamespace("reticulate", quietly=TRUE), "\n")'
```

**Expected:**

```
DevKidCC: 0.5.1
Seurat  : 5.x.x
knn.iter present: TRUE
sceasy  : TRUE
reticulate: TRUE
```

`knn.iter present: FALSE` means the image built the R package from the wrong
branch — stop here and tell me.

Then run the image's own self-test, which drives a synthetic matrix through both
the R and the Python entry point. It needs no input file and no network, so it
works on a login node, and it is the same check CI runs before publishing:

```bash
singularity exec dkcc-v051.sif python /opt/smoke_test.py
```

**Expected:** ends with `PASS: both the R and Python entry points classified the
input.` The `LineageID` lines above it will read `{'unassigned': 200}` — the
input is noise, so that is the correct answer, and the point of the test is that
nothing crashed.

If this fails, stop and send me the output; nothing below it will work either.

And confirm the scripts inside the image are the fixed ones:

```bash
# reads h5ad via sceasy rather than the silently-failing SeuratDisk::Convert
singularity exec dkcc-v051.sif grep -c sceasy /opt/run_dkcc.R              # expect 4

# the runtime DKCC() rewrite is gone — the only remaining hit is the comment
# explaining that it used to be there, so no *code* line may match
singularity exec dkcc-v051.sif \
    grep -c '^[[:space:]]*[^#[:space:]].*assignInNamespace' /opt/run_dkcc.R  # expect 0

# batch mode enumerates up front instead of reprocessing its own outputs
singularity exec dkcc-v051.sif grep -c mapfile /opt/run_dkcc_batch.sh      # expect 1
```

## 4. Get a test file across

The reference input is Howden 2019 (5,365 cells, 166 MB), which has a known-good
answer from the host run:

```bash
# from the homeserver
rsync -avP \
  /data/homeserver/data/datasets/Howden_2019_Organoids/processed/Howden_2019_Organoids_qc.h5ad \
  <user>@<hpc-login>:/group/$USER/dkcc_v051_test/
```

Any h5ad of raw counts with HGNC symbols works if that one is awkward to move.

## 5. Single-file run via SLURM

```bash
cd "$DKCC_TEST"
git clone -b containerise-v0.5.1 https://github.com/KidneyRegeneration/DevKidCC dkcc-repo
cd dkcc-repo

./run_dkcc.sh \
    --mode slurm \
    --container "$DKCC_TEST/dkcc-v051.sif" \
    --partition prod_short \
    --time 02:00:00 \
    --mem 32G \
    --cpus 4 \
    --mounts /group \
    --logs "$DKCC_TEST" \
    "$DKCC_TEST/Howden_2019_Organoids_qc.h5ad"
```

5,365 cells took ~90 s on the homeserver, so `prod_short` with 32 G is generous.
Scale `--mem` up for larger inputs; the h5ad path holds the whole matrix in R.

Watch it:

```bash
squeue -u $USER
tail -f "$DKCC_TEST"/dkcc_Howden_2019_Organoids_qc_*.log
```

**Expected in the log:**

```
Job started: ...
Node      : ...
Reading input: /data/Howden_2019_Organoids_qc.h5ad
...
Running DKCC classification...
Saving output...
Job complete: ...
```

**Expected on disk:** `Howden_2019_Organoids_qc_DKCC.h5ad`, next to the input.

## 6. Check the output

```bash
singularity exec "$DKCC_TEST/dkcc-v051.sif" python -c "
import anndata as ad
a = ad.read_h5ad('$DKCC_TEST/Howden_2019_Organoids_qc_DKCC.h5ad')
print(a.shape)
print(a.obs['LineageID'].value_counts())
print(a.obs['DKCC'].value_counts().head(10))
"
```

**Expected** — the homeserver run of the same file, KNN on, for comparison:

| LineageID | n |
|---|---|
| Stroma | 2,733 |
| Nephron | 1,828 |
| NPC | 402 |
| NPC-like | 273 |
| unassigned | 122 |
| UrEp | 4 |
| Endo | 3 |

Near-identical is the pass condition, not bit-identical: the container pins
`r-base=4.4` against the host's R 4.5.3, and the container path reads the h5ad
through sceasy rather than the wrapper's CSV handoff. A few dozen cells moving is
fine; whole classes appearing or vanishing is not.

## 7. Batch mode (optional but worth it)

This exercises the de-race fix — the old loop reprocessed its own outputs.

```bash
mkdir -p "$DKCC_TEST/batch"
cp "$DKCC_TEST/Howden_2019_Organoids_qc.h5ad" "$DKCC_TEST/batch/sample_a.h5ad"
cp "$DKCC_TEST/Howden_2019_Organoids_qc.h5ad" "$DKCC_TEST/batch/sample_b.h5ad"

./run_dkcc.sh --mode slurm --container "$DKCC_TEST/dkcc-v051.sif" \
    --partition prod_short --time 02:00:00 --mem 32G --cpus 4 \
    --logs "$DKCC_TEST" "$DKCC_TEST/batch"

ls "$DKCC_TEST/batch"
```

**Expected:** exactly four files — `sample_a.h5ad`, `sample_b.h5ad`,
`sample_a_DKCC.h5ad`, `sample_b_DKCC.h5ad`. Anything named `*_DKCC_DKCC.h5ad`
means the fix did not take.

## 8. Send back

- The full SLURM log from §5
- The output of §3 (both blocks) and §6
- The `ls` from §7
- If anything failed: `sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed`

An `OUT_OF_MEMORY` state is a resource problem, not a code problem — re-run with
`--mem 128G --partition himem`. Anything else, send it over and I'll fix it.

---

## Known rough edges

- **`--container remote` pulls `:latest`.** There is no flag for an arbitrary
  tag, hence the manual pull in §2. Once this branch merges and `:latest` moves,
  `--container remote` is the right call again.
- **`--mounts` defaults to `/group`.** Add `--mounts /group,/scratch` if the data
  lives elsewhere; paths that don't exist are skipped with a warning rather than
  failing the run.
- **`R_PROFILE_USER=/dev/null`** is set at exec so a `~/.Rprofile` on the cluster
  (renv, in particular) cannot hijack the library paths inside the container.
