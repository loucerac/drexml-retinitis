#!/usr/bin/env bash

set -e

set -a; source exp_design.env; set +a


THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

ATC_PATH="${THIS_FOLDER}/data/raw/WHO ATC-DDD 2021-12-03.csv"

# Python-based
PY_ENV="${THIS_FOLDER}/.venvs/py"
PY_SRC_PATH="${THIS_FOLDER}/scripts/py/drug-ora.py"
DB_NAME=${genes#"genes-"}
DB_PATH="${THIS_FOLDER}/data/final/${DB_NAME}"
echo "${DB_PATH}"
conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PY_SRC_PATH} "${DB_PATH}" "${ATC_PATH}"

# R-based
R_ENV="${THIS_FOLDER}/.venvs/r"

R_SRC_PATH="${THIS_FOLDER}/scripts/r/04_read_filter_shap_drugbank.R"
echo "Running ${R_FNAME}"
DB_FNAME="genes_drugbank-${DRUGBANK_VERSION}_mygene-${MYGENE_VERSION}.tsv"
conda run --no-capture-output --live-stream -p ${R_ENV} Rscript ${R_SRC_PATH} ${DB_FNAME}

fnames=(
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
