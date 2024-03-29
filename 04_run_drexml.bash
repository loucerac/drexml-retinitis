#!/usr/bin/env bash

set -e

set -a; source .env; set +a

if [[ -z "$USE_GPU" ]]; then
    n_gpus=0
else
    n_gpus=$USE_GPU
fi

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PY_ENV="${THIS_FOLDER}/.venvs/drexml"
EXP_PATH="${THIS_FOLDER}/exp_design.env"
conda run --no-capture-output --live-stream -p ${PY_ENV} drexml run --n-gpus $n_gpus --n-cpus -1 ${EXP_PATH}
ML_FOLDER="${THIS_FOLDER}/results/ml"
mkdir -p ${ML_FOLDER}
mv "${THIS_FOLDER}/results/shap"*".tsv" ${ML_FOLDER}/.
mv "${THIS_FOLDER}/results/stability"*".tsv" ${ML_FOLDER}/.
conda run --no-capture-output --live-stream -p ${PY_ENV} drexml rename ${ML_FOLDER}
conda run --no-capture-output --live-stream -p ${PY_ENV} drexml plot ${ML_FOLDER}/stability_results_symbol.tsv
rm -rf "${THIS_FOLDER}/results/tmp"
