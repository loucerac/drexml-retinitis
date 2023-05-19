#!/usr/bin/env bash

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PY_ENV="${THIS_FOLDER}/.venvs/drexml"
EXP_PATH="${THIS_FOLDER}/exp_design.env"
conda run -p ${PY_ENV} drexml run --n-gpus 0 --n-cpus -1 ${EXP_PATH}
ML_FOLDER="${THIS_FOLDER}/results/ml"
mkdir -p ${ML_FOLDER}
mv "${THIS_FOLDER}/results/shap*.tsv" ${ML_FOLDER}/.
mv "${THIS_FOLDER}/results/stability_results*.tsv" ${ML_FOLDER}/.
conda run -p ${PY_ENV} drexml rename ${ML_FOLDER}
conda run -p ${PY_ENV} drexml plot ${ML_FOLDER}/stability_results_symbol.tsv
rm -rf "${THIS_FOLDER}/results/tmp"
