#!/usr/bin/env bash

set -e

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

set -a; source .env; set +a

PY_ENV_FOLDER="$THIS_FOLDER/.venvs/py"
conda env remove -p ${PY_ENV_FOLDER} 
conda env create -p ${PY_ENV_FOLDER} -f "$THIS_FOLDER/environment_py.yml"

DREXML_ENV_FOLDER="$THIS_FOLDER/.venvs/drexml"
if [[ -z "$USE_GPU" ]]; then
    n_gpus=0
else
    n_gpus=$USE_GPU
fi

if [[ "$n_gpus" -gt 0 ]]; then
    if [ "$UPDATE" = false ] ; then
        env_file="environment_drexml_gpu.yml"
    else
        env_file="environment_update_drexml_gpu.yml"
    fi
else
    if [ "$UPDATE" = false ] ; then
        env_file="environment_drexml.yml"
    else
        env_file="environment_update_drexml.yml"
    fi
fi

echo "using ${env_file}"

conda env remove -p ${DREXML_ENV_FOLDER}
conda env create -p ${DREXML_ENV_FOLDER} -f "${THIS_FOLDER}/${env_file}"

R_ENV_FOLDER="$THIS_FOLDER/.venvs/r"
conda env remove -p ${R_ENV_FOLDER}
conda env create -p ${R_ENV_FOLDER} -f "$THIS_FOLDER/environment_r.yml"

# create data folders
mkdir -p "$THIS_FOLDER/data/raw"
mkdir -p "$THIS_FOLDER/data/interim"
mkdir -p "$THIS_FOLDER/data/final"
mkdir -p "$THIS_FOLDER/results/tables"
mkdir -p "$THIS_FOLDER/results/ml"
mkdir -p "$THIS_FOLDER/results/rds"
mkdir -p "$THIS_FOLDER/results/figures"
