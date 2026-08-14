#!/usr/bin/env bash

# Run DKCC on a single file or recursively on a folder.
#
# Usage (single file):
#   bash /opt/run_dkcc_batch.sh -i sample.h5ad [-o h5ad|rds]
#
# Usage (folder):
#   bash /opt/run_dkcc_batch.sh -i /path/to/folder [-o h5ad|rds]
#
# Optional:
#   -r /path/to/run_dkcc.R   override the run script location

set -euo pipefail

OUTPUT_FORMAT="h5ad"
RSCRIPT_PATH="/opt/run_dkcc.R"

usage() {
    echo ""
    echo "Usage: bash run_dkcc_batch.sh -i <input_file_or_folder> [-o h5ad|rds] [-r run_dkcc.R]"
    echo ""
    exit 1
}

while getopts "i:o:r:" opt; do
    case ${opt} in
        i) INPUT_PATH="$OPTARG" ;;
        o) OUTPUT_FORMAT="$OPTARG" ;;
        r) RSCRIPT_PATH="$OPTARG" ;;
        *) usage ;;
    esac
done

if [[ -z "${INPUT_PATH:-}" ]]; then
    usage
fi

SUPPORTED_EXTENSIONS=("*.h5ad" "*.h5" "*.h5seurat" "*.rds" "*.Rds" "*.RData" "*.rdata")

process_file() {
    local infile="$1"
    local dir base stem outfile

    dir=$(dirname "$infile")
    base=$(basename "$infile")

    case "$base" in
        *.h5ad)    stem="${base%.h5ad}" ;;
        *.h5seurat) stem="${base%.h5seurat}" ;;
        *.RData)   stem="${base%.RData}" ;;
        *.rdata)   stem="${base%.rdata}" ;;
        *)         stem="${base%.*}" ;;
    esac

    outfile="${dir}/${stem}_DKCC.${OUTPUT_FORMAT}"

    echo "----------------------------------------"
    echo "Processing:"
    echo "  Input : $infile"
    echo "  Output: $outfile"
    echo "----------------------------------------"

    Rscript "$RSCRIPT_PATH" \
        --input  "$infile" \
        --output "$outfile" \
        --format "$OUTPUT_FORMAT"
}

if [[ -f "$INPUT_PATH" ]]; then
    process_file "$INPUT_PATH"
elif [[ -d "$INPUT_PATH" ]]; then
    # Enumerate the full input list up front, before any processing starts.
    # Streaming `find | while read` here is racy: process_file() writes its
    # "<stem>_DKCC.<ext>" output into the same directory being scanned, and
    # since outputs share an extension with the inputs (e.g. both are
    # .h5ad), a still-running `find` can pick up an output file that was
    # just written and reprocess it as if it were new input. Excluding the
    # "_DKCC." naming convention is a second layer of defense in case this
    # folder already contains previous outputs alongside fresh input files.
    mapfile -t files < <(
        for pattern in "${SUPPORTED_EXTENSIONS[@]}"; do
            find "$INPUT_PATH" -type f -name "$pattern"
        done | grep -v '_DKCC\.' | sort -u
    )
    for file in "${files[@]}"; do
        process_file "$file"
    done
else
    echo "ERROR: Input path does not exist: $INPUT_PATH"
    exit 1
fi

echo ""
echo "All processing complete."
