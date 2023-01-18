#!/usr/bin/env bash

XML_PATH="/home/carlos/github/drexml-retinitis/data/raw/drugbank_v050108.xml.gz"
INTERIM_FOLDER="/home/carlos/github/drexml-retinitis/data/interim"
FINAL_FOLDER="/home/carlos/github/drexml-retinitis/data/final"

PARSER_FOLDER="/home/carlos/github/drexml-retinitis/modules/drugbank-parser"

python ${PARSER_FOLDER}/parser.py $XML_PATH $INTERIM_FOLDER $FINAL_FOLDER
