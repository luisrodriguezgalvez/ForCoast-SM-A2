#!/bin/bash
source /root/.bashrc

# Clean input folder
rm /data/*.nc

# Download data
cd /opt/forcoast/preprocessing
python forcoast_download.py $1 $2 $3
# Run model
cd /opt/forcoast/processing
python forcoast.py $1 $2 $3
# Run postprocessing
cd /opt/forcoast/postprocessing
python SM-A2-Postprocess.py -y SMA2-PostProcess-Eforie.yaml -d $3
