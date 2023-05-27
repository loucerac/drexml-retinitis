#!/usr/bin/env bash

set -e

set -a; source .env; set +a

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PARSER_ENV_FOLDER="$THIS_FOLDER/.venvs/drugbank-parser"
conda env remove -p ${PARSER_ENV_FOLDER} 
conda env create -p ${PARSER_ENV_FOLDER} -f "$THIS_FOLDER/environment_drugbank-parser.yml"

DREXML_ENV_FOLDER="$THIS_FOLDER/.venvs/drexml"
if [[ -z "$USE_GPU" ]]; then
    env_file="environment_drexml.yml"
elif [[ "$USE_GPU" -gt 0 ]]; then
    env_file="environment_drexml_gpu.yml"
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
touch -a .env
