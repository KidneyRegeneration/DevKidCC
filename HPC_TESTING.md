# HPC test checklist — DevKidCC container v0.5.1

Run this on the MCRI cluster to confirm the `containerise-v0.5.1` image works
under Apptainer/SLURM. Everything below is copy-paste; §7 says what to send back.

**Prerequisite:** the branch is pushed and the `docker-publish.yml` run for
`containerise-v0.5.1` is green, so
`ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1` exists.

**What is being tested that wasn't before:** two things.

1. **The image is self-contained.** Previous runs bind-mounted fixed copies of
   `run_dkcc.R` / `run_dkcc_batch.sh` over the ones inside the SIF. Those fixes
   are baked in now, and `build_binds()` no longer mounts any script. If a step
   below fails in a way that would once have been papered over by a bind-mount,
   that is the finding — don't re-add the mount, send the log.
2. **The file now chooses the entry point.** `.h5ad` goes to Python
   (`/opt/run_dkcc.py`), R-native formats go to R (`/opt/run_dkcc.R`), and
   `/opt/dkcc` routes on the extension. Previously *everything* went to R, and an
   h5ad was read by R calling back into Python through reticulate — which is
   where every h5ad failure in this image came from.

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

**If it dies with `stream ID <n>; PROTOCOL_ERROR`**, that is an HTTP/2 transport
fault between Singularity and GHCR, not a permissions or disk problem. Retrying
does not help: Singularity restarts the whole pull each time rather than resuming,
so a 4 GB image never gets to the end. This happened on every attempt during local
testing. Use skopeo instead — it reconnects *inside* a blob and carries on from
the byte it reached, which is the difference that gets a 4 GB image home:

```bash
# Drop to HTTP/1.1 (both are Go binaries, so GODEBUG reaches their HTTP stack)
export GODEBUG=http2client=0

skopeo copy --retry-times 10 \
    docker://ghcr.io/kidneyregeneration/dkcc:containerise-v0.5.1 \
    oci:$PWD/dkcc-oci:v051

# Then build the SIF from the local copy -- no network involved.
# Note: no :tag here. Singularity treats everything after `oci:` as a path,
# so `oci:$PWD/dkcc-oci:v051` looks for a directory literally named "dkcc-oci:v051".
singularity build dkcc-v051.sif oci:$PWD/dkcc-oci
```

Expect log lines like `Reading blob body ... failed (unexpected EOF), reconnecting
after 1614807040 bytes` — that is the resume working, not an error.

`skopeo` is widely available on HPC systems; if it is not, `module spider skopeo`
or ask the helpdesk. The local OCI directory can be deleted once the SIF is built.

## 3. Sanity-check the image (login node, no SLURM)

```bash
singularity exec dkcc-v051.sif Rscript -e '
  cat("DevKidCC:", as.character(packageVersion("DevKidCC")), "\n")
  cat("Seurat  :", as.character(packageVersion("Seurat")), "\n")
  cat("knn.iter present:", "knn.iter" %in% names(formals(DevKidCC::DKCC)), "\n")'

singularity exec dkcc-v051.sif python -c '
import devkidcc, anndata
print("wrapper :", devkidcc.__file__)
print("anndata :", anndata.__version__)'
```

**Expected:**

```
DevKidCC: 0.5.1
Seurat  : 5.x.x
knn.iter present: TRUE
wrapper : /opt/micromamba/envs/devkid/lib/python3.12/site-packages/devkidcc/__init__.py
anndata : 0.x.x
```

`knn.iter present: FALSE` means the image built the R package from the wrong
branch — stop here and tell me.

Then run the image's own self-test, which drives a synthetic matrix through the
Python API, an h5ad through `/opt/dkcc` (asserting it routes to Python) and an
`.rds` through `/opt/dkcc` (asserting it routes to R). It needs no input file and
no network, so it works on a login node, and it is the same check CI runs before
publishing:

```bash
singularity exec dkcc-v051.sif python /opt/smoke_test.py
```

**Expected:** ends with `PASS: the Python API, h5ad routing and .rds routing all
classified the input.` The `LineageID` lines above it will read
`{'unassigned': 200}` — the input is noise, so that is the correct answer, and
the point of the test is that nothing crashed. Takes about three minutes.

If this fails, stop and send me the output; nothing below it will work either.

And confirm the scripts inside the image are the fixed ones:

```bash
# the R script no longer reads or writes h5ad at all: no sceasy, no reticulate,
# no scCustomize::as.anndata — that work belongs to Python now
singularity exec dkcc-v051.sif \
    grep -c 'library(sceasy)\|convertFormat\|as\.anndata(' /opt/run_dkcc.R  # expect 0

# and the routing exists, with a Python side to route to
singularity exec dkcc-v051.sif ls -l /opt/dkcc /opt/run_dkcc.py

# the runtime DKCC() rewrite is gone — the only remaining hit is the comment
# explaining that it used to be there, so no *code* line may match
singularity exec dkcc-v051.sif \
    grep -c '^[[:space:]]*[^#[:space:]].*assignInNamespace' /opt/run_dkcc.R  # expect 0

# batch mode enumerates up front instead of reprocessing its own outputs
singularity exec dkcc-v051.sif grep -c mapfile /opt/run_dkcc_batch.sh      # expect 1
```

Routing is worth one direct check too, since a regression here would still
classify and so would look like a pass everywhere else:

```bash
singularity exec dkcc-v051.sif /opt/dkcc --input /tmp/x.h5ad --output /tmp/y.h5ad 2>&1 | head -1
# expect: Routing to the Python entry point (.h5ad is AnnData's format)
# (it then fails on the missing file — that is fine, the routing line is the point)
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
Scale `--mem` up for larger inputs. An h5ad input goes through the Python
wrapper, which projects onto the reference genes before handing the counts to R,
so it holds rather less in R than the old all-R path did.

`--format` no longer takes a value that contradicts the input: h5ad in gives
h5ad out, an R-native file gives `.rds` out. Passing a mismatched `--format` is
now an error rather than a silent conversion.

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
`r-base=4.4` against the host's R 4.5.3. Both now reach R the same way — the h5ad
goes through the Python wrapper's CSV handoff on either machine — so the R
version is the only remaining source of drift. A few dozen cells moving is fine;
whole classes appearing or vanishing is not.

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
