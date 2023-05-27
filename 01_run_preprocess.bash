#!/usr/bin/env bash

set -e

set -a; source .env; set +a

if [ "$UPDATE" = true ] ; then
    MYGENE_VERSION=$( date +%Y%m%d )
    echo "Update Mygene version to ${MYGENE_VERSION}"
else
    MYGENE_VERSION="20230120"
    echo "Using Mygene version ${MYGENE_VERSION}"
fi

THIS_FOLDER=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
DATA_PATH="${THIS_FOLDER}/data"
RAW_FOLDER="${DATA_PATH}/raw"
mkdir -p $RAW_FOLDER
INTERIM_FOLDER="${DATA_PATH}/interim"
mkdir -p $INTERIM_FOLDER
FINAL_FOLDER="${DATA_PATH}/final"
mkdir -p $FINAL_FOLDER

# R
R_ENV="${THIS_FOLDER}/.venvs/r"
R_FNAME="00_GTEx_processing.R"
R_SRC_PATH="${THIS_FOLDER}/scripts/r/${R_FNAME}"

# Python
PY_ENV="${THIS_FOLDER}/.venvs/py"

# RAW
XML_PATH="${RAW_FOLDER}/drugbank-${DRUGBANK_VERSION}.xml.gz"
PARQUET_PATH="${FINAL_FOLDER}/expreset_Hinorm_gtex${GTEX_VERSION}.rds.feather"

# INTERIM
TSV_PATH="${INTERIM_FOLDER}/drugbank-${DRUGBANK_VERSION}.tsv"
DB_GENES_NAME="genes_drugbank-${DRUGBANK_VERSION}_mygene-${MYGENE_VERSION}.tsv"
DB_GENES_PATH="${INTERIM_FOLDER}/${DB_GENES_NAME}"
GTEX_GENES_NAME="genes_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"
GTEX_GENES_PATH="${INTERIM_FOLDER}/${GTEX_GENES_NAME}"

# FINAL
DB_FILTERED_PATH="${FINAL_FOLDER}/drugbank-${DRUGBANK_VERSION}_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"
GENES_FILTERED_PATH="${FINAL_FOLDER}/genes-drugbank-${DRUGBANK_VERSION}_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"

PARSER_FOLDER="${THIS_FOLDER}/scripts/py"

# Parsing
conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PARSER_FOLDER}/parser.py parse $XML_PATH $TSV_PATH

if test -f "$DB_GENES_PATH"; then
    echo "$DB_GENES_PATH exists."
else
    echo "$DB_GENES_PATH does no exist."
    if [ "$UPDATE" = true ] ; then
        echo "Updating"
        conda run --no-capture-output --live-stream -p ${R_ENV} Rscript ${R_SRC_PATH} ${GTEX_FNAME} ${GTEX_VERSION}
        conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PARSER_FOLDER}/parser.py translate $TSV_PATH $DB_GENES_PATH --kind drugbank 
        conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PARSER_FOLDER}/parser.py translate $PARQUET_PATH $GTEX_GENES_PATH --kind gtex 
    else
        cp "${RAW_FOLDER}/${DB_GENES_NAME}" $DB_GENES_PATH
        cp "${RAW_FOLDER}/${GTEX_GENES_NAME}" $GTEX_GENES_PATH

        if test -f "$DB_GENES_PATH"; then
            echo "Found $DB_GENES_PATH in assets."
        else
            echo "Either search for the corresponding Mygene $MYGENE_VERSION file or set UPDATE to true. Then, launch again the script."
        fi

    fi
fi

conda run --no-capture-output --live-stream -p ${PY_ENV} python ${PARSER_FOLDER}/parser.py filter \
    --drugbank-path ${TSV_PATH} \
    --drugbank-genes-path ${DB_GENES_PATH} \
    --gtex-genes-path ${GTEX_GENES_PATH} \
    ${DB_FILTERED_PATH} \
    ${GENES_FILTERED_PATH}


# recreate design environment after update
if [ "$UPDATE" = true ] ; then
    rm -f $THIS_FOLDER/exp_design.env
    echo "###### EXPERIMENT DESIGN  ######" >> $THIS_FOLDER/exp_design.env
    echo "data_path=${THIS_FOLDER}/data/final/" >> $THIS_FOLDER/exp_design.env
    echo "gene_exp=expreset_Hinorm_gtex${GTEX_VERSION}.rds.feather" >> $THIS_FOLDER/exp_design.env
    echo "pathvals=expreset_pathvals_gtex${GTEX_VERSION}.rds.feather" >> $THIS_FOLDER/exp_design.env 
    echo "circuits=circuits_RP.rds.feather" >> $THIS_FOLDER/exp_design.env
    echo "circuits_column=in_disease" >> $THIS_FOLDER/exp_design.env
    echo "genes=$GENES_FILTERED_PATH" >> $THIS_FOLDER/exp_design.env
    echo "genes_column=drugbank_approved_targets" >> $THIS_FOLDER/exp_design.env
else
    cp -f $THIS_FOLDER/exp_design_default.env $THIS_FOLDER/exp_design.env
fi
