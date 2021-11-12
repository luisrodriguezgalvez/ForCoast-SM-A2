#!/bin/bash
# Example: run.sh eforie 2021-09-20 10 
# Example: run.sh eforie 2021-09-20 10 

#source activate env
INITIAL_DIR="$(pwd)"

# Clean data folder
#rm /usr/src/app/data/*.*
DATA_DIR=/home/arthur/Desktop/DOCS/PROJECTS/FORECOAST/Pilot/SM/ForCoast-SM-A3/testdata/

mkdir -p ${DATA_DIR}

## Download data (using data dir from yml)
cd ./PreProcessing
echo "python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}"
echo "-->Skipped"
#python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}

## Run model
## This is done once per source
cd ../Processing
sourcecount=0
while read -r s; do
    echo "python forcoast.py -y $1 -T $2 -p $3 -s $s -d ${DATA_DIR} "
    echo "-->Skipped"
    #python forcoast.py -y $1 -T $2 -p $3 -s $s -d ${DATA_DIR}
    #mv ../usr/$1/output/test.nc ../usr/$1/output/test_${sourcecount}.nc
    sourcecount=`expr $sourcecount + 1`
done < "../usr/$1/config/sources.txt"

## Run postprocessing
cd ../PostProcessing
sourcecount=0
targetcount=0
while read -r t; do # Loop on targets
    while read -r s; do # Loop on Sources
	echo "python SM-A2-Postprocess.py -y $1 -s $s -sc $sourcecount -t $t -tc $targetcount -d ${DATA_DIR}"
	python SM-A2-Postprocess.py -y $1 -s $s -c $sourcecount -t $t -k $targetcount -d ${DATA_DIR}
        sourcecount=`expr $sourcecount + 1`
    done < "../usr/$1/config/sources.txt"
    targetcount=`expr $targetcount + 1`
done < "../usr/$1/config/targets.txt"


# Run bulletin spcript 
cd ../BulletinScript
sourcecount=0
targetcount=0
while read -r t; do
    while read -r s; do
	echo "python bulletin_script.py -y $1 -s $s -sc $sourcecount -t $t -tc $targetcount -d ${DATA_DIR}"
	python bulletin_script.py -y $1 -s $s -c $sourcecount -t $t -k $targetcount -d ${DATA_DIR}
        sourcecount=`expr $sourcecount + 1`
    done < "../usr/$1/config/sources.txt"
    targetcount=`expr $targetcount + 1`
done < "../usr/$1/config/targets.txt"

# Remains TODO a multi-source bulletin per target.
## Can be done by collecting on the TS_Risk.png from each usr/${1}/output/source_X_target_Y/ directories

echo $INITIAL_DIR
cd ..
ls

#cp /usr/src/app/data/*.png $INITIAL_DIR
