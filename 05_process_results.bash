#!/usr/bin/env bash

set -e

set -a; source .env; set +a

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Python-based
PY_ENV="${THIS_FOLDER}/.venvs/py"
PY_SRC_PATH="${THIS_FOLDER}/scripts/py/drug-ora.py"
conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PY_SRC_PATH}

# R-based
R_ENV="${THIS_FOLDER}/.venvs/r"

fnames=(
    "04_read_filter_shap_drugbank.R"
    "05_heatmaps_balloon_tables.R"
    "06_RPhallmarks_analysis.R"
    "07_Drugs_ATC_statistics.R"
    "08_KDT_Kmeans_Clusters.R"
    "09_Drugs_SpiderPlot.R"
    "10_circle_chord_diagram.R"
)

for R_FNAME in ${fnames[@]}; do
    R_SRC_PATH="${THIS_FOLDER}/scripts/r/${R_FNAME}"
    echo "Running ${R_FNAME}"
    conda run --no-capture-output --live-stream -p ${R_ENV} Rscript ${R_SRC_PATH}
    wait
done
