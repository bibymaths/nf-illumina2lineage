#!/usr/bin/env bash
set -euo pipefail

# 1) Create & activate Conda env (if not already)
ENV_NAME="nf-illumina2reads"
ENV_FILE="envs/environment.yaml"

echo "Checking Conda environment..."
if ! conda info --envs | grep -q "^${ENV_NAME}"; then
  echo "-> Creating environment ${ENV_NAME} from ${ENV_FILE}"
  conda env create -n "${ENV_NAME}" --file "${ENV_FILE}"
fi
echo "-> Activating ${ENV_NAME}"
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# 2) Install & configure Java via SDKMAN
if [ ! -d "$HOME/.sdkman" ]; then
  echo "💡 Installing SDKMAN…"
  curl -s "https://get.sdkman.io" | bash
fi
# shellcheck disable=SC1090
source "$HOME/.sdkman/bin/sdkman-init.sh"
JAVA_VER="17.0.8-tem"
if ! sdk current java | grep -q "$JAVA_VER"; then
  echo "-> Installing Java $JAVA_VER"
  sdk install java "$JAVA_VER"
  sdk default java "$JAVA_VER"
fi

# 3) Clear Conda's Java so Nextflow picks up SDKMAN's
unset JAVA_CMD JAVA_HOME

# 4) Run Nextflow
echo "Running Nextflow…"
nextflow run main.nf \
  --reads "${PWD}/data/*.fastq.gz" \
  --reference "${PWD}/data/reference.fasta" \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html
