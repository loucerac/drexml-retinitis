#!/usr/bin/env bash

set -e

set -a; source .env; set +a

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CONDA_ENV="${THIS_FOLDER}/.venvs/drugbank-parser"
FPATH="${THIS_FOLDER}/modules/drugbank-parser/downloader.py"
conda run --no-capture-output --live-stream -p ${CONDA_ENV} python ${FPATH}
