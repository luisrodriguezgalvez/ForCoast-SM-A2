#!/bin/bash
# Example: run.sh eforie 2021-09-20 10 [28.6803,44.2201] [28.6611,44.0535,28.6617,44.0553] 
# Example: run.sh eforie 2021-09-20 10 [28.7303,44.2201] [28.6611,44.0535,28.6617,44.0553]

#source activate env

INITIAL_DIR="$(pwd)"

# Clean data folder
#rm /usr/src/app/data/*.*

DATA_DIR=/home/arthur/Desktop/DOCS/PROJECTS/FORECOAST/Pilot/SM/ForCoast-SM-A3/testdata/

mkdir ${DATA_DIR}

# Download data (using data dir from yml)
cd ./PreProcessing
#echo "python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}"
#python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}

# Run model
cd ../Processing
echo "python forcoast.py -y $1 -T $2 -p $3 -s $4 -t $5 -d ${DATA_DIR} "
python forcoast.py -y $1 -T $2 -p $3 -s $4 -t $5 -d ${DATA_DIR}

# Run postprocessing (uses data directory from yaml file)

cd ../PostProcessing
echo "python SM-A2-Postprocess.py -y $1.yaml -s $4 -t $5 -d ${DATA_DIR}"
python SM-A2-Postprocess.py -y $1.yaml -s $4 -t $5 -d ${DATA_DIR}

# Run bulletin spcript (uses data directory from yaml file)
cd ../BulletinScript
echo "bulletin_script.py -y $1.yaml -s $4 -t $5 -d ${DATA_DIR}"
python bulletin_script.py -y $1.yaml -s $4 -t $5 -d ${DATA_DIR}

echo $INITIAL_DIR
cd ..
ls

#cp /usr/src/app/data/*.png $INITIAL_DIR
