#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Fast glycan masking library from existing Rosetta models
# ============================================================
#
# Edit the USER SETTINGS section below, then run:
#
#   bash run_existing_models_library.sh
#
# MODE:
#   analyze   -> only compare pos_* total_score values with pos_0
#   empirical -> build a library directly from Rosetta glycans
#   expanded  -> empirical library + additional torsion sampling
#

# -------------------- USER SETTINGS --------------------------

MODE="empirical"

# Directory containing pos_0, pos_1, pos_2, ...
INPUT_DIR="/path/to/all_pos_folders"

# Output prefix. Parent directory is created automatically.
OUT_PREFIX="libraries/Man5_empirical"

# Position-selection criterion:
# mean(total_score, pos_i) - mean(total_score, pos_0) < REU_DELTA
REU_DELTA="0"

# Optional known native N-glycan ASN positions.
# Examples:
#   NATIVE_POSITIONS="200"
#   NATIVE_POSITIONS="141,200,305"
#   NATIVE_POSITIONS=""
NATIVE_POSITIONS=""

# Include glycans from pos_0 itself in the library?
# pos_0 is always used as the WT score baseline.
INCLUDE_WT_GLYCANS="false"

# Optional expansion settings. Used only when MODE="expanded".
SAMPLES_PER_GLYCAN="4"
ATTACHMENT_SIGMA_DEG="20"
GLYCOSIDIC_SIGMA_DEG="20"

# Reproducible random seed.
SEED="2026"

# Optional maximum final library size.
# Leave empty for no limit.
MAX_LIBRARY_SIZE=""

# Multi-model PDBs can become very large. NPZ is always written.
WRITE_PDB="false"

# If true, stop immediately when an invalid/mismatched glycan is found.
# If false, skip it and record the reason in *_skipped.tsv.
STRICT="false"

# Analyze-only report name.
ANALYZE_OUT="position_scores.csv"

# ------------------------------------------------------------

if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "[ERROR] INPUT_DIR does not exist: ${INPUT_DIR}" >&2
    exit 1
fi

mkdir -p "$(dirname "${OUT_PREFIX}")"

COMMON_ARGS=(
    --input-dir "${INPUT_DIR}"
    --reu-delta "${REU_DELTA}"
)

case "${MODE}" in
    analyze)
        echo "[INFO] Running position score analysis only..."
        python -m fast_glycan_masking_from_existing_models analyze \
            "${COMMON_ARGS[@]}" \
            --out "${ANALYZE_OUT}"
        ;;

    empirical|expanded)
        BUILD_ARGS=(
            "${COMMON_ARGS[@]}"
            --out-prefix "${OUT_PREFIX}"
            --seed "${SEED}"
        )

        if [[ -n "${NATIVE_POSITIONS}" ]]; then
            BUILD_ARGS+=(--native-positions "${NATIVE_POSITIONS}")
        fi

        if [[ "${INCLUDE_WT_GLYCANS}" == "true" ]]; then
            BUILD_ARGS+=(--include-wt-glycans)
        fi

        if [[ -n "${MAX_LIBRARY_SIZE}" ]]; then
            BUILD_ARGS+=(--max-library-size "${MAX_LIBRARY_SIZE}")
        fi

        if [[ "${WRITE_PDB}" == "true" ]]; then
            BUILD_ARGS+=(--write-pdb)
        fi

        if [[ "${STRICT}" == "true" ]]; then
            BUILD_ARGS+=(--strict)
        fi

        if [[ "${MODE}" == "expanded" ]]; then
            BUILD_ARGS+=(
                --samples-per-glycan "${SAMPLES_PER_GLYCAN}"
                --attachment-sigma-deg "${ATTACHMENT_SIGMA_DEG}"
                --glycosidic-sigma-deg "${GLYCOSIDIC_SIGMA_DEG}"
            )
        fi

        echo "[INFO] MODE=${MODE}"
        echo "[INFO] INPUT_DIR=${INPUT_DIR}"
        echo "[INFO] OUT_PREFIX=${OUT_PREFIX}"
        echo "[INFO] REU_DELTA=${REU_DELTA}"

        python -m fast_glycan_masking_from_existing_models build \
            "${BUILD_ARGS[@]}"
        ;;

    *)
        echo "[ERROR] MODE must be analyze, empirical, or expanded." >&2
        exit 1
        ;;
esac

echo "[DONE]"
