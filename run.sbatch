#!/usr/bin/env bash

#SBATCH --job-name=drexml-rp
#SBATCH --cpus-per-task=50 
#SBATCH --nodes=1
#SBATCH --mem=300g
#SBATCH --time=72:00:00
#SBATCH --output=./logs/drexml-rp_%j.out.log
#SBATCH --error=./logs/drexml-rp_%j.err.log

. ~/miniconda3/etc/profile.d/conda.sh;

start=`date +%s`
bash install.sh
wait
end=`date +%s`
runtime=$((end-start))
echo "install duration: $runtime"

start=`date +%s`
bash 00_run_hipathia.bash
echo "already run"
wait
end=`date +%s`
runtime=$((end-start))
echo "hipathia duration: $runtime"

start=`date +%s`
bash 01_run_parser.bash
wait
end=`date +%s`
runtime=$((end-start))
echo "parser duration: $runtime"

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
bash 04_run_drexml.bash
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
