#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# EDIT CELL RANGER, REFERENCE, FEATURE REFERENCE, AND OUTPUT PATHS
# ============================================================
CELLRANGER=/path/to/cellranger-8.0.0/cellranger
REFERENCE=/path/to/refdata-gex-GRCh38-2024-A
FEATURE_REF=/path/to/feature_reference.csv
OUTDIR=/path/to/cellranger_results
CORES=16
MEMORY=128
# ============================================================

# ============================================================
# EDIT HIGH-MOI SAMPLE NAMES AND LIBRARIES CSV PATHS
# ============================================================
HIGH_SAMPLES=(
    high_reaction1
    high_reaction2
)

HIGH_LIBRARIES=(
    /path/to/high_reaction1_libraries.csv
    /path/to/high_reaction2_libraries.csv
)
# ============================================================

# ============================================================
# EDIT LOW-MOI SAMPLE NAMES AND LIBRARIES CSV PATHS
# ============================================================
LOW_SAMPLES=(
    low_reaction1
    low_reaction2
)

LOW_LIBRARIES=(
    /path/to/low_reaction1_libraries.csv
    /path/to/low_reaction2_libraries.csv
)
# ============================================================

SAMPLES=("${HIGH_SAMPLES[@]}" "${LOW_SAMPLES[@]}")
LIBRARIES=("${HIGH_LIBRARIES[@]}" "${LOW_LIBRARIES[@]}")

[[ ${#SAMPLES[@]} -eq ${#LIBRARIES[@]} ]] || {
    echo "The number of sample names and libraries CSV files must match." >&2
    exit 1
}

mkdir -p "$OUTDIR"
cd "$OUTDIR"

for i in "${!SAMPLES[@]}"; do
    "$CELLRANGER" count \
        --id="${SAMPLES[$i]}" \
        --transcriptome="$REFERENCE" \
        --libraries="${LIBRARIES[$i]}" \
        --feature-ref="$FEATURE_REF" \
        --localcores="$CORES" \
        --localmem="$MEMORY"
done

echo "sample_id,molecule_h5" > aggr.csv
> sceptre_matrix_dirs.txt

for sample in "${SAMPLES[@]}"; do
    echo "$sample,$OUTDIR/$sample/outs/molecule_info.h5" >> aggr.csv
    echo "$OUTDIR/$sample/outs/filtered_feature_bc_matrix" >> sceptre_matrix_dirs.txt
done

"$CELLRANGER" aggr \
    --id=combined_high_low_moi \
    --csv=aggr.csv \
    --normalize=none \
    --localcores="$CORES" \
    --localmem="$MEMORY"
