#!/usr/bin/env bash
# Create the `piezo1` conda environment with the full scientific stack.
# Usage:  bash scripts/create_env.sh   (then: conda activate piezo1)
set -euo pipefail
ENV_NAME="${1:-piezo1}"
CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
if command -v mamba >/dev/null 2>&1; then SOLVER=mamba; else SOLVER=conda; fi

echo "==> creating env '$ENV_NAME' (python 3.11) with $SOLVER"
$SOLVER create -y -n "$ENV_NAME" -c conda-forge python=3.11

echo "==> conda-forge scientific core"
$SOLVER install -y -n "$ENV_NAME" -c conda-forge \
    numpy scipy matplotlib pandas numba networkx \
    scikit-image scikit-learn biopython \
    mdanalysis mdtraj openmm pdbfixer biotite \
    requests pyyaml tqdm pytest pytest-qt ruff

echo "==> pip layer (GUI + GL + structural bioinformatics)"
conda run -n "$ENV_NAME" python -m pip install --upgrade pip
conda run -n "$ENV_NAME" python -m pip install \
    "PyQt6>=6.6" "PyQt6-Qt6>=6.6" moderngl PyOpenGL PyOpenGL-accelerate \
    pyqtgraph ProDy pydssp freesasa

echo "==> done. activate with:  conda activate $ENV_NAME"
