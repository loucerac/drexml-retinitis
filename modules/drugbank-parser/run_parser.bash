#!/usr/bin/env bash
DATA_PATH="/home/carlos/github/drexml-retinitis/data"
RAW_FOLDER="${DATA_PATH}/raw"
INTERIM_FOLDER="${DATA_PATH}/interim"
FINAL_FOLDER="${DATA_PATH}/final"

XML_PATH="${RAW_FOLDER}/drugbank_v050108.xml.gz"
TSV_PATH="${INTERIM_FOLDER}/drugbank_v050108.tsv"
PARQUET_PATH="${FINAL_FOLDER}/expreset_Hinorm_gtexV8.rds.feather"

PARSER_FOLDER="/home/carlos/github/drexml-retinitis/modules/drugbank-parser"

python ${PARSER_FOLDER}/parser.py $INTERIM_FOLDER $FINAL_FOLDER parse $XML_PATH 
python ${PARSER_FOLDER}/parser.py $INTERIM_FOLDER $FINAL_FOLDER translate $TSV_PATH --kind drugbank 
python ${PARSER_FOLDER}/parser.py $INTERIM_FOLDER $FINAL_FOLDER translate $PARQUET_PATH --kind gtex 
