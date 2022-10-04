#!/bin/bash

#(C) Copyright FORCOAST H2020 project under Grant No. 870465. All rights reserved.

# Example: ./run_py.sh eforie 2022-03-02 3

# Example: ./run_py.sh eforie 2022-03-02 3 github config source target datadir
#                      $1     $2         $3 $4     $5     $6     $7     $8

# Example: ./run_py.sh eforie 2022-03-02 3 github https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/config.yaml https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/sources.txt https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/targets.txt /usr/src/app/data/
# Example: ./run_py.sh eforie 2022-03-02 3 cli local [28.6803,44.2201,0.500] [[28.6611,44.0435],[28.67,44.0435],[28.67,44.0553],[28.6611,44.0553]] /usr/src/app/data/

# Example: ./run_py.sh eforie 2022-03-02 3 github config source target datadir token chat_id
#                      $1     $2         $3 $4     $5     $6     $7     $8     $9    $10

# Example: ./run_py.sh eforie 2022-03-02 3 github https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/config.yaml https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/sources.txt https://raw.githubusercontent.com/FORCOAST/ForCoast-A2-Settings/User_XYZ/targets.txt /usr/src/app/data/ 5267228188:AAGx60FtWgHkScBb3ISFL1dp6Oq_9z9z0rw -1001780197306

doPreProcessing=true
doProcessing=true
doPostProcessing=true
doBulletin=true

INITIAL_DIR="$(pwd)"

cd /usr/src/app

# Get source and target files from github based on links provided
if [[ "$4" = "github" ]] ; then
    echo 'Get files from github'
    wget -O ./usr/$1/config/config.yaml $5
    wget -O ./usr/$1/config/sources.txt $6
    echo sources.txt
    wget -O ./usr/$1/config/targets.txt $7
    echo targets.txt
fi

# Substitute values in source ($4) and target ($5) files
if [[ "$4" = "cli" ]] ; then
    echo 'Substitute values in sources.txt and targets.txt'
    > ./usr/$1/config/sources.txt
    echo $6 >> ./usr/$1/config/sources.txt
    > ./usr/$1/config/targets.txt
    echo $7 >> ./usr/$1/config/targets.txt
fi

# Set default data folder if not defined in input arguments
if [ "$#" -eq 3 ]; then
    DATA_DIR=/usr/src/app/data/
else
    DATA_DIR=$8
fi

mkdir -p ${DATA_DIR}

## Download data (using data dir from yml)
cd ./PreProcessing
echo "python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}"
if $doPreProcessing ; then
    python forcoast_download_yml.py -a $1 -T $2 -p $3 -d ${DATA_DIR}
else
    echo "-->Skipped"
fi

echo ''
echo '###########'
echo "Download done."
echo '###########'

## Run model
## This is done once per source
cd ../Processing
echo '###########'
echo "Starting Processing, according to ../usr/$1/config/sources.txt"
echo '###########'
sourcecount=0
while read -r s; do
    echo "python forcoast.py -y $1 -T $2 -p $3 -s $s -d ${DATA_DIR} "
    if $doProcessing; then
	python forcoast.py -y $1 -T $2 -p $3 -s $s -d ${DATA_DIR}
	if [ -e ../usr/$1/output/test.nc ]; then
	    mv ../usr/$1/output/test.nc ../usr/$1/output/test_${sourcecount}.nc
	else
	    echo "CRITICAL ERROR: PARCELS DID NOT GENERATE THE EXPECTED OUTPUT FILE"; exit
	fi
    else
	echo "-->Skipped"
    fi
    sourcecount=`expr $sourcecount + 1`
done < "../usr/$1/config/sources.txt"
echo '###########'
echo "Processing Done"
echo '###########'

## Run postprocessing
cd ../PostProcessing
echo '###########'
echo "Starting Post-Processing, according to ../usr/$1/config/sources.txt and ../usr/$1/config/targets.txt"
echo '###########'
targetcount=0
while read -r t; do # Loop on targets
    sourcecount=0
    while read -r s; do # Loop on Sources
	echo "python SM-A2-Postprocess.py -y $1 -s $s -sc $sourcecount -t $t -tc $targetcount -d ${DATA_DIR}"
	if $doPostProcessing; then
	    python SM-A2-Postprocess.py -y $1 -s $s -c $sourcecount -t $t -k $targetcount -d ${DATA_DIR}
	    # Generate Animation
	    convert -delay 50 -fuzz 2% -layers Optimize -loop 0 ../usr/$1/output/target_${targetcount}_source_${sourcecount}/AllTracks_Alarm*.png ../usr/$1/output/target_${targetcount}_source_${sourcecount}/AllTracks.gif
	else
            echo "-->Skipped"
	fi
	sourcecount=`expr $sourcecount + 1`    
    done < "../usr/$1/config/sources.txt"
    targetcount=`expr $targetcount + 1`
done < "../usr/$1/config/targets.txt"
echo '###########'
echo "Post - Processing Done"
echo '###########'


# Run bulletin spcript 
cd ../BulletinScript
echo '###########'
echo "Starting Bulletin Generation, according to ../usr/$1/config/sources.txt and ../usr/$1/config/targets.txt"
echo '###########'
targetcount=0
while read -r t; do
    sourcecount=0
    while read -r s; do
    	echo "python bulletin_script.py -y $1 -T $2 -s $s -sc $sourcecount -t $t -tc $targetcount -d ${DATA_DIR}"
    	if $doBulleting ; then
	        python bulletin_script.py -y $1 -T $2 -s $s -c $sourcecount -t $t -k $targetcount -d ${DATA_DIR}
            # Generate Animation
    	    convert -delay 50 -fuzz 2% -layers Optimize -loop 0 ../usr/$1/output/target_${targetcount}_source_${sourcecount}/bulletin*.png ../usr/$1/output/target_${targetcount}_source_${sourcecount}/bulletin.gif
	    else
	        echo "-->Skipped"
        fi
        cp ../usr/$1/output/target_${targetcount}_source_${sourcecount}/bulletin.mp4 ../usr/$1/output/bulletin_target_${targetcount}_source_${sourcecount}.mp4
        cp ../usr/$1/output/target_${targetcount}_source_${sourcecount}/bulletin.webm ../usr/$1/output/bulletin_target_${targetcount}_source_${sourcecount}.webm
        sourcecount=`expr $sourcecount + 1`
    done < "../usr/$1/config/sources.txt"
    targetcount=`expr $targetcount + 1`
done < "../usr/$1/config/targets.txt"
echo '###########'
echo "Bulletin Done"
echo '###########'

echo '###########'
echo "Share bulletin through Telegram"
echo '###########'

if [ "$#" -eq 10 ]; then

    cd ../Telegram
    python send_bulletin.py -T $9 -C ${10} -B /usr/src/app/usr/$1/output/target_0_source_0/bulletin.mp4 -M video

fi

# Remains TODO a multi-source bulletin per target.
## Can be done by collecting on the TS_Risk.png from each usr/${1}/output/source_X_target_Y/ directories

echo $INITIAL_DIR
cd ..
ls

cp /usr/src/app/usr/$1/output/bulletin*mp4 $INITIAL_DIR/bulletin.mp4
cp /usr/src/app/usr/$1/output/bulletin*mp4 ${DATA_DIR}/bulletin.mp4

cp /usr/src/app/usr/$1/output/bulletin*webm $INITIAL_DIR/bulletin.webm
cp /usr/src/app/usr/$1/output/bulletin*webm ${DATA_DIR}/bulletin.webm

rm /usr/src/app/data/*.nc