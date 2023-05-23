#!/usr/bin/env bash

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CONDA_ENV="${THIS_FOLDER}/.venvs/drugbank-parser"
FPATH="${THIS_FOLDER}/modules/drugbank-parser/downloader.py"
conda run -p ${CONDA_ENV} python ${FPATH}
