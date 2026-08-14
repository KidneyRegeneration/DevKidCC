#!/usr/bin/env bash
#
# run_dkcc.sh — Launch DevKidCC cell classification inside a container
#
# Usage:
#   ./run_dkcc.sh [OPTIONS] <input_file_or_folder>
#
# Options:
#   -f, --format      Output format: h5ad | rds              [default: h5ad]
#   -c, --container   remote | <path/to/image.sif>           [default: remote]
#                       remote       = pull ghcr.io/kidneyregeneration/dkcc:latest
#                       <path.sif>   = use an existing local .sif file
#   -s, --sif         Where to save the pulled .sif file     [default: ./dkcc.sif]
#                     Ignored when --container points to an existing .sif
#   -m, --mode        local | slurm                          [default: local]
#                       local  = run immediately with Singularity
#                       slurm  = submit as a SLURM job using Apptainer
#   --mounts          Comma-separated host paths to bind into the container
#                     [default: /group]
#                     e.g. --mounts /group,/scratch
#
#   SLURM options (--mode slurm only):
#   -p, --partition   Partition name   [default: prod_med]
#        prod_short   max  2 h
#        prod_med     max 24 h  ← default
#        prod_long    max 14 d
#        himem        max  7 d  (high-memory nodes)
#   --mem             Memory per job   [default: 64G]
#   --cpus            CPUs per task    [default: 8]
#   --time            Time limit       [default: 23:00:00]
#   --job-name        Job name         [default: dkcc_<input_stem>]
#   --logs            Log directory    [default: same directory as input]

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

REMOTE_URI="ghcr.io/kidneyregeneration/dkcc:latest"

FORMAT="h5ad"
CONTAINER_OPT="remote"
MODE="local"
SIF="${SCRIPT_DIR}/dkcc.sif"

MOUNTS="/group"

PARTITION="prod_med"
MEM="64G"
CPUS="8"
TIME="23:00:00"
JOB_NAME=""
LOG_DIR=""

# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------

usage() {
    sed -n '3,32p' "$0" | sed 's/^# \?//'
    exit 1
}

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

INPUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--format)    FORMAT="$2";        shift 2 ;;
        -c|--container) CONTAINER_OPT="$2"; shift 2 ;;
        -s|--sif)       SIF="$2";           shift 2 ;;
        -m|--mode)      MODE="$2";          shift 2 ;;
        -p|--partition) PARTITION="$2";     shift 2 ;;
        --mem)          MEM="$2";           shift 2 ;;
        --cpus)         CPUS="$2";          shift 2 ;;
        --time)         TIME="$2";          shift 2 ;;
        --mounts)       MOUNTS="$2";        shift 2 ;;
        --job-name)     JOB_NAME="$2";      shift 2 ;;
        --logs)         LOG_DIR="$2";       shift 2 ;;
        -h|--help)      usage ;;
        -*) echo "ERROR: Unknown option: $1"; echo ""; usage ;;
        *)  INPUT="$1"; shift ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "ERROR: No input specified."
    echo ""
    usage
fi

INPUT=$(realpath "$INPUT")

# ---------------------------------------------------------------------------
# Validate
# ---------------------------------------------------------------------------

if [[ "$FORMAT" != "h5ad" && "$FORMAT" != "rds" ]]; then
    echo "ERROR: --format must be h5ad or rds (got: $FORMAT)"
    exit 1
fi

if [[ "$MODE" != "local" && "$MODE" != "slurm" ]]; then
    echo "ERROR: --mode must be local or slurm (got: $MODE)"
    exit 1
fi

# ---------------------------------------------------------------------------
# Resolve SIF path
# ---------------------------------------------------------------------------

if [[ "$CONTAINER_OPT" == "remote" ]]; then
    : # SIF stays as the pull target; resolved in ensure_sif()
else
    # User passed an explicit .sif path
    SIF="$CONTAINER_OPT"
    if [[ ! -f "$SIF" ]]; then
        echo "ERROR: .sif file not found: $SIF"
        exit 1
    fi
fi

ensure_sif() {
    if [[ -f "$SIF" ]]; then
        echo "Using existing SIF: $(realpath "$SIF")"
        return
    fi

    if ! command -v singularity &>/dev/null; then
        echo "ERROR: singularity not found. Cannot pull image."
        exit 1
    fi

    echo "Pulling container to ${SIF} ..."
    echo "  Source : docker://${REMOTE_URI}"
    echo "  Dest   : ${SIF}"
    echo ""
    singularity pull "$SIF" "docker://${REMOTE_URI}"
}

# ---------------------------------------------------------------------------
# Resolve inner container command from input path
# ---------------------------------------------------------------------------

if [[ -f "$INPUT" ]]; then
    DATA_DIR=$(dirname "$INPUT")
    BASE=$(basename "$INPUT")

    case "$BASE" in
        *.h5ad)     STEM="${BASE%.h5ad}" ;;
        *.h5seurat) STEM="${BASE%.h5seurat}" ;;
        *.RData)    STEM="${BASE%.RData}" ;;
        *.rdata)    STEM="${BASE%.rdata}" ;;
        *)          STEM="${BASE%.*}" ;;
    esac

    INNER_CMD="Rscript /opt/run_dkcc.R \
        --input  /data/${BASE} \
        --output /data/${STEM}_DKCC.${FORMAT} \
        --format ${FORMAT}"

    [[ -z "$JOB_NAME" ]] && JOB_NAME="dkcc_${STEM}"

elif [[ -d "$INPUT" ]]; then
    DATA_DIR="$INPUT"
    INNER_CMD="bash /opt/run_dkcc_batch.sh -i /data -o ${FORMAT}"

    [[ -z "$JOB_NAME" ]] && JOB_NAME="dkcc_batch_$(basename "$INPUT")"

else
    echo "ERROR: Input does not exist: $INPUT"
    exit 1
fi

[[ -z "$LOG_DIR" ]] && LOG_DIR="$DATA_DIR"

# ---------------------------------------------------------------------------
# Build --bind flags from DATA_DIR + MOUNTS
# ---------------------------------------------------------------------------

build_binds() {
    local flags="--bind ${DATA_DIR}:/data"

    # No script bind-mounts here: the h5ad-reading, Assay5, orig.ident/PAX2 and
    # batch-loop fixes are baked into the image, so the SIF is self-contained.

    IFS=',' read -ra mnt_list <<< "$MOUNTS"
    for mnt in "${mnt_list[@]}"; do
        mnt="${mnt## }"   # trim leading spaces
        mnt="${mnt%% }"   # trim trailing spaces
        if [[ -z "$mnt" ]]; then
            continue
        elif [[ -e "$mnt" ]]; then
            flags="$flags --bind ${mnt}:${mnt}"
        else
            echo "WARNING: skipping --mounts entry '$mnt' (does not exist on this host)" >&2
        fi
    done
    echo "$flags"
}

# ---------------------------------------------------------------------------
# Local mode — Singularity
# ---------------------------------------------------------------------------

run_local() {
    if ! command -v singularity &>/dev/null; then
        echo "ERROR: singularity not found."
        exit 1
    fi

    ensure_sif
    SIF=$(realpath "$SIF")
    BIND_FLAGS=$(build_binds)

    echo "Running DevKidCC locally via Singularity"
    echo "  Input  : $INPUT"
    echo "  Format : $FORMAT"
    echo "  Image  : $SIF"
    echo "  Mounts : /data (data dir), ${MOUNTS}"
    echo ""

    singularity exec \
        ${BIND_FLAGS} \
        --env R_PROFILE_USER=/dev/null \
        --env RETICULATE_PYTHON=/opt/micromamba/envs/devkid/bin/python \
        "$SIF" bash -c "$INNER_CMD"
}

# ---------------------------------------------------------------------------
# SLURM mode — Apptainer
# ---------------------------------------------------------------------------

run_slurm() {
    if ! command -v sbatch &>/dev/null; then
        echo "ERROR: sbatch not found — are you on the HPC login node?"
        exit 1
    fi

    # Pull the SIF on the login node now so compute nodes don't need internet
    ensure_sif
    SIF=$(realpath "$SIF")
    BIND_FLAGS=$(build_binds)

    local SBATCH_SCRIPT
    SBATCH_SCRIPT=$(mktemp /tmp/dkcc_XXXXXX.sh)

    cat > "$SBATCH_SCRIPT" <<SBATCH
#!/usr/bin/env bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --mem=${MEM}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --time=${TIME}
#SBATCH --output=${LOG_DIR}/${JOB_NAME}_%j.log
#SBATCH --error=${LOG_DIR}/${JOB_NAME}_%j.log

set -euo pipefail

module load apptainer

echo "Job started: \$(date)"
echo "Node      : \$(hostname)"
echo "Input     : ${INPUT}"
echo "Image     : ${SIF}"
echo ""

apptainer exec \
    ${BIND_FLAGS} \
    --env R_PROFILE_USER=/dev/null \
    --env RETICULATE_PYTHON=/opt/micromamba/envs/devkid/bin/python \
    "${SIF}" bash -c "${INNER_CMD}"

echo ""
echo "Job complete: \$(date)"
SBATCH

    echo "Submitting SLURM job: $JOB_NAME"
    echo "  Input     : $INPUT"
    echo "  Format    : $FORMAT"
    echo "  Image     : $SIF"
    echo "  Partition : $PARTITION"
    echo "  Memory    : $MEM"
    echo "  CPUs      : $CPUS"
    echo "  Time      : $TIME"
    echo "  Log       : ${LOG_DIR}/${JOB_NAME}_<jobid>.log"
    echo ""

    sbatch "$SBATCH_SCRIPT"
    rm -f "$SBATCH_SCRIPT"
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

case "$MODE" in
    local) run_local ;;
    slurm) run_slurm ;;
esac
