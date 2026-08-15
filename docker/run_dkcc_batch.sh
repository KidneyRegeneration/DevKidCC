#!/usr/bin/env bash

# Run DKCC on a single file or recursively on a folder.
#
# Usage (single file):
#   bash /opt/run_dkcc_batch.sh -i sample.h5ad
#
# Usage (folder):
#   bash /opt/run_dkcc_batch.sh -i /path/to/folder
#
# Each file is routed by extension: .h5ad goes to the Python entry point,
# R-native formats to the R one. The output format therefore follows the input
# rather than being chosen for the whole folder -- a folder may hold both.
#
# Optional:
#   -r /path/to/dkcc   override the dispatcher location

set -euo pipefail

DISPATCH_PATH="/opt/dkcc"

usage() {
    echo ""
    echo "Usage: bash run_dkcc_batch.sh -i <input_file_or_folder> [-r /opt/dkcc]"
    echo ""
    exit 1
}

while getopts "i:o:r:" opt; do
    case ${opt} in
        i) INPUT_PATH="$OPTARG" ;;
        o) echo "NOTE: -o is ignored; each file's output format follows its input." ;;
        r) DISPATCH_PATH="$OPTARG" ;;
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

    local out_format
    case "$base" in
        *.h5ad)     stem="${base%.h5ad}";     out_format="h5ad" ;;
        *.h5seurat) stem="${base%.h5seurat}"; out_format="rds"  ;;
        *.RData)    stem="${base%.RData}";    out_format="rds"  ;;
        *.rdata)    stem="${base%.rdata}";    out_format="rds"  ;;
        *)          stem="${base%.*}";        out_format="rds"  ;;
    esac

    outfile="${dir}/${stem}_DKCC.${out_format}"

    echo "----------------------------------------"
    echo "Processing:"
    echo "  Input : $infile"
    echo "  Output: $outfile"
    echo "----------------------------------------"

    "$DISPATCH_PATH" \
        --input  "$infile" \
        --output "$outfile"
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
