#!/usr/bin/env bash
# Create (or update) the conda environment for rnaseq_pipeline.sh
#   ./install.sh            -> environment named "rnaseq"
#   ./install.sh myenv      -> environment named "myenv"
#   ./install.sh -p /path   -> environment at a prefix (e.g. shared project space)
set -euo pipefail

cd "$(dirname "$0")"

if command -v mamba >/dev/null 2>&1; then
    CONDA=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA=conda
else
    echo "ERROR: conda/mamba not found. Install Miniforge first:" >&2
    echo "  https://github.com/conda-forge/miniforge" >&2
    exit 1
fi

if [[ ${1:-} == "-p" ]]; then
    [[ -n ${2:-} ]] || { echo "Usage: $0 -p /path/to/env" >&2; exit 1; }
    TARGET=(-p "$2")
    ACTIVATE=$2
    EXISTS=$([[ -d $2/conda-meta ]] && echo 1 || echo 0)
else
    NAME=${1:-rnaseq}
    TARGET=(-n "$NAME")
    ACTIVATE=$NAME
    EXISTS=$(conda env list | awk '{print $1}' | grep -qx "$NAME" && echo 1 || echo 0)
fi

if [[ $EXISTS == 1 ]]; then
    echo "Updating existing environment '$ACTIVATE' with $CONDA ..."
    "$CONDA" env update "${TARGET[@]}" -f environment.yml --prune
else
    echo "Creating environment '$ACTIVATE' with $CONDA ..."
    "$CONDA" env create "${TARGET[@]}" -f environment.yml
fi

chmod +x rnaseq_pipeline.sh

cat <<EOF

Done. Next steps:
  conda activate $ACTIVATE
  $(pwd)/rnaseq_pipeline.sh --check
EOF
