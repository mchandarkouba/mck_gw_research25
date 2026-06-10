#!/bin/bash

module purge
module load anaconda3/2020.07
source /opt/packages/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda activate ***env_name

cd ***/bot/dir/

echo -e "env activated \n"
