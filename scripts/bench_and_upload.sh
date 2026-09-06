#!/usr/bin/env bash
set -euo pipefail

# Script to run large-dataset Chazelle triangulation benchmarks and upload results to Bencher (bencher.dev).

BENCH_CMD="cargo bench --bench bench_triangulation"

# Check if bencher CLI is installed
if ! command -v bencher &> /dev/null; then
    echo "Error: 'bencher' CLI is not found in PATH."
    echo "Install it via: curl --proto '=https' --tlsv1.2 -sSfL https://bencher.dev/download/install-cli.sh | sh"
    exit 1
fi

echo "================================================================="
echo "  Chazelle Linear-Time Triangulation Benchmark Suite & Bencher   "
echo "================================================================="

# Check for Bencher credentials
if [[ -n "${BENCHER_API_TOKEN:-}" || -n "${BENCHER_API_KEY:-}" ]] && [[ -n "${BENCHER_PROJECT:-}" ]]; then
    echo "Bencher project: ${BENCHER_PROJECT}"
    echo "Running benchmarks and uploading results to bencher.dev..."
    
    BRANCH_ARGS=()
    if [[ -n "${BENCHER_BRANCH:-}" ]]; then
        BRANCH_ARGS+=(--branch "${BENCHER_BRANCH}")
    fi
    if [[ -n "${BENCHER_TESTBED:-}" ]]; then
        BRANCH_ARGS+=(--testbed "${BENCHER_TESTBED}")
    fi

    ${BENCH_CMD} | bencher run \
        --project "${BENCHER_PROJECT}" \
        --adapter json \
        "${BRANCH_ARGS[@]}" \
        "$@"
    
    echo "Upload complete! View your results on https://bencher.dev/perf/${BENCHER_PROJECT}"
else
    echo "NOTE: BENCHER_PROJECT and BENCHER_API_TOKEN / BENCHER_API_KEY are not set."
    echo "Running benchmark in dry-run mode (no data uploaded to bencher.dev)..."
    echo ""
    echo "To upload to Bencher:"
    echo "  export BENCHER_PROJECT=\"your-project-slug\""
    echo "  export BENCHER_API_TOKEN=\"your-token\" (or BENCHER_API_KEY=\"...\")"
    echo "  ./scripts/bench_and_upload.sh"
    echo ""

    ${BENCH_CMD} | bencher run --adapter json --dry-run
fi
