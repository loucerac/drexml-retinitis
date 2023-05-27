#!/usr/bin/env bash

set -e

set -a; source .env; set +a

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PY_ENV_FOLDER="${THIS_FOLDER}/.venvs/py"
FPATH="${THIS_FOLDER}/scripts/py/downloader.py"
conda run --no-capture-output --live-stream -p ${PY_ENV_FOLDER} python ${FPATH}
