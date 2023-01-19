#!/usr/bin/env bash
DATA_PATH="/home/carlos/github/drexml-retinitis/data"
RAW_FOLDER="${DATA_PATH}/raw"
INTERIM_FOLDER="${DATA_PATH}/interim"
FINAL_FOLDER="${DATA_PATH}/final"

GTEX_VERSION="V8"
DRUGBANK_VERSION="v050108"
MYGENE_VERSION=$(date +%Y%m%d)

# RAW
XML_PATH="${RAW_FOLDER}/drugbank-${DRUGBANK_VERSION}.xml.gz"
PARQUET_PATH="${FINAL_FOLDER}/expreset_Hinorm_gtex${GTEX_VERSION}.rds.feather"

# INTERIM
TSV_PATH="${INTERIM_FOLDER}/drugbank-${DRUGBANK_VERSION}.tsv"
DB_GENES_PATH="${INTERIM_FOLDER}/genes_drugbank-${DRUGBANK_VERSION}_mygene-${MYGENE_VERSION}.tsv"
GTEX_GENES_PATH="${INTERIM_FOLDER}/genes_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"

# FINAL
DB_FILTERED_PATH="${FINAL_FOLDER}/drugbank-${DRUGBANK_VERSION}_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"
DB_FILTERED_PATH="${FINAL_FOLDER}/genes-${DRUGBANK_VERSION}_gtex-${GTEX_VERSION}_mygene-${MYGENE_VERSION}.tsv"

PARSER_FOLDER="/home/carlos/github/drexml-retinitis/modules/drugbank-parser"

python ${PARSER_FOLDER}/parser.py parse $XML_PATH $TSV_PATH
python ${PARSER_FOLDER}/parser.py translate $TSV_PATH $DB_GENES_PATH --kind drugbank 
python ${PARSER_FOLDER}/parser.py translate $PARQUET_PATH $GTEX_GENES_PATH --kind gtex 
python ${PARSER_FOLDER}/parser.py filter \
 --drugbank-path $TSV_PATH \
 --drugbank-genes-path $DB_GENES_PATH \ 
 --gtex-path $PARQUET_PATH \
 --gtex-genes-path $GTEX_GENES_PATH \
 $DB_FILTERED_PATH

