PARSER_ENV_FOLDER="./.venvs/drugbank-parser"
conda env update --prune -p ${PARSER_ENV_FOLDER} -f environment_drugbank-parser.yml

DREXML_ENV_FOLDER="./.venvs/drexml"
conda env update --prune -p ${DREXML_ENV_FOLDER} -f environment_drexml.yml

R_ENV_FOLDER="./.venvs/r"
conda env update --prune -p ${R_ENV_FOLDER} -f environment_r.yml
export PKG_CONFIG_PATH="${R_ENV_FOLDER}/lib/pkgconfig"
export PKG_LIBS="-liconv"
conda run -p .venvs/r/ R --vanilla -e 'source("renv/activate.R"); renv::restore()'
