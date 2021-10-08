#!/bin/bash
# Example: run.sh eforie 2021-09-20 10 [28.6803,44.2201] [28.6611,44.0535,28.6617,44.0553] 
# Example: run.sh eforie 2021-09-20 10 [28.7303,44.2201] [28.6611,44.0535,28.6617,44.0553]

source activate env

INITIAL_DIR="$(pwd)"

# Clean data folder
rm /usr/src/app/data/*.*

# Download data (using data dir from yml)
cd /usr/src/app/preprocessing
echo "python forcoast_download_yml.py -T $2 -p $3"
python forcoast_download_yml.py -T $2 -p $3 
# Run model
cd /usr/src/app/processing
echo "python forcoast.py -y $1 -T $2 -p $3 -s $4 -t $5"
python forcoast.py -y $1 -T $2 -p $3 -s $4 -t $5
# Run postprocessing (uses data directory from yaml file)
cd /usr/src/app/postprocessing
echo "python SM-A2-Postprocess.py -y $1.yaml -s $4 -t $5"
python SM-A2-Postprocess.py -y $1.yaml -s $4 -t $5

echo $INITIAL_DIR
cd ..
ls

cp *.nc $INITIAL_DIR
