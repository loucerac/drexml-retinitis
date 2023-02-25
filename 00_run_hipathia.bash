#!/usr/bin/env bash

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
R_ENV="${THIS_FOLDER}/.venvs/r"
R_FNAME="00_GTEx_V8processing.R"
R_SRC_PATH="${THIS_FOLDER}/modules/r/${R_FNAME}"
conda run -p ${R_ENV} Rscript ${R_SRC_PATH}
