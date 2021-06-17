#!/bin/bash
# Example: run.sh A2 eforie 2021-02-01 1 /data
source /root/.bashrc

# Clean input folder
rm /data/*.nc

# Download data
cd /opt/forcoast/preprocessing
python forcoast_download_yml.py $1 $2 $3 $4 $5
# Run model
cd /opt/forcoast/processing
python forcoast.py $1 $2 $3 $4 $5
# Run postprocessing
cd /opt/forcoast/postprocessing
# python SM-A2-Postprocess.py -y $2.yaml -d $5
