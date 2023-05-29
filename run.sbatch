#!/usr/bin/env bash

#SBATCH --job-name=drexml-rp
#SBATCH --cpus-per-task=50 
#SBATCH --nodes=1
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH --output=./logs/drexml-rp_%j.out.log
#SBATCH --error=./logs/drexml-rp_%j.err.log


if ! test -f .env; then
    echo "Copying default dotenv file from project."
    cp ./example.env .env
fi

set -a; source .env; set +a

if command -v conda >/dev/null 2>&1; then
    source $(conda info --base)/etc/profile.d/conda.sh
    echo "conda already available"
elif test -f ~/miniconda3/etc/profile.d/conda.sh; then
    . ~/miniconda3/etc/profile.d/conda.sh;
    echo "conda already available, sourcing installed conda"
else
    echo "Installing miniconda into the project's workspace"
    WORK_DIR=`mktemp -d`
    (cd $WORK_DIR & curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh)
    (cd $WORK_DIR & bash Miniconda3-latest-Linux-x86_64.sh -b -u -p ./.venvs/miniconda)
    source ./.venvs/miniconda/bin/activate
fi

echo "Using conda: $(which conda)"

start=`date +%s`
bash 00a_install.sh $USE_GPU
wait
end=`date +%s`
runtime=$((end-start))
echo "install duration: $runtime"

start=`date +%s`
bash 00b_download_data.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "download duration: $runtime"

start=`date +%s`
bash 01_run_preprocess.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "preprocess duration: $runtime"

start=`date +%s`
bash 02_run_annotator.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "annotator duration: $runtime"

start=`date +%s`
bash 03_run_map_builder.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "builder duration: $runtime"

start=`date +%s`
bash 04_run_drexml.bash $USE_GPU
wait
end=`date +%s`
runtime=$((end-start))
echo "drexml duration: $runtime"

start=`date +%s`
bash 05_process_results.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "analysis duration: $runtime"
