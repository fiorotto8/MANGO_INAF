#!/bin/bash
set -euo pipefail

INPUT_FILE="${1:-../build/output_files/Sr90_t0.root}"
PARTICLE="${2:-e-}"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file does not exist: $INPUT_FILE" >&2
    exit 1
fi

INPUT_BASENAME="$(basename "$INPUT_FILE")"
OUTPUT_FILE="elab_${INPUT_BASENAME}"

root -l -b -q "RecoTrack.C(\"${INPUT_FILE}\",\"${PARTICLE}\",true)"

echo "Executing angular study on ${OUTPUT_FILE}..."
./study "$OUTPUT_FILE"
