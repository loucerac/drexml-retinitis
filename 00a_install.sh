#!/usr/bin/env bash

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

PARSER_ENV_FOLDER="$THIS_FOLDER/.venvs/drugbank-parser"
conda env update --prune -p ${PARSER_ENV_FOLDER} -f "$THIS_FOLDER/environment_drugbank-parser.yml"

DREXML_ENV_FOLDER="$THIS_FOLDER/.venvs/drexml"
conda env update --prune -p ${DREXML_ENV_FOLDER} -f "$THIS_FOLDER/environment_drexml.yml"

R_ENV_FOLDER="$THIS_FOLDER/.venvs/r"
conda env update --prune -p ${R_ENV_FOLDER} -f "$THIS_FOLDER/environment_r.yml"
export PKG_CONFIG_PATH="${R_ENV_FOLDER}/lib/pkgconfig"
#export PKG_LIBS="-liconv"

conda run -p $R_ENV_FOLDER R --vanilla -e "source('$THIS_FOLDER/renv/activate.R'); renv::restore()"

# create data folders
mkdir -p "$THIS_FOLDER/data/raw"
mkdir -p "$THIS_FOLDER/data/interim"
mkdir -p "$THIS_FOLDER/data/final"
