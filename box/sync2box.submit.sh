#!/bin/bash
MY_DIR=$(dirname $0)
MY_NAME=$(basename $0 .submit.sh)
DIR_NAME=$(echo $PWD | cut -d / -f 5- | tr / .)

mkdir -p logs
sbatch --job-name=${MY_NAME}.${DIR_NAME} ${MY_DIR}/${MY_NAME}.job.sh
