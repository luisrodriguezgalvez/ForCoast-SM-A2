#!/bin/bash
# Example: run.sh eforie 2021-09-20 10 /data [28.6803,44.2201] [28.6611,44.0535,28.6617,44.0553]
# Example: run.sh eforie 2021-09-20 10 /data [28.7303,44.2201] [28.6611,44.0535,28.6617,44.0553]

# Clean input folder
rm $4/*.nc

# Download data
cd /opt/forcoast/preprocessing
python forcoast_download_yml.py -T $2 -p $3 
# Run model
cd /opt/forcoast/processing
python forcoast.py $1 $2 $3 $4 $5 $6
# Run postprocessing
cd /opt/forcoast/postprocessing
python SM-A2-Postprocess.py -y $1.yaml -d $4 -s $5 -t $6
