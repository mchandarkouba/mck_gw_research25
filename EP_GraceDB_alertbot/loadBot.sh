#!/bin/bash

module purge
module load anaconda3/2020.07
source /opt/packages/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate ***env_name

cd ***/bot/dir/
# running the rest of our program
echo -e "crossmatch bot running\n"

nohup python -u ***bot_script.py > botoutput.txt 2> botoutput_err.txt < /dev/null &
echo $! > save_pid.txt
