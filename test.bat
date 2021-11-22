:: Download input (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PreProcessing
python forcoast_download_yml.py -a eforie -T 2021-08-15 -p 5 -d c:\data
cd ..

:: Run model (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd Processing
python forcoast.py -y eforie -T 2021-08-15 -p 5 -s [28.6803,44.2201,0.5] -d c:\data\
cd ..

:: TO BE MODIFIED FOR FLEXIBILITY Change name of test.nc to test_1.nc (1 is because the sourcecount -c in this case in the Postprocessing will be 1)
cd usr\eforie\output
ren test.nc test_1.nc
cd ..\..\..

:: Run postprocessing (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PostProcessing
python SM-A2-Postprocess.py -y eforie -s [28.6803,44.2201] -c 1 -t [[28.6611,44.0435],[28.67,44.0435],[28.67,44.0553],[28.6611,44.0553]] -k 0 -d c:\data 
cd ..

:: Run bulletin script
cd BulletinScript
python bulletin_script.py -y eforie -s [28.6803,44.2201] -c 1 -t [[28.6611,44.0435],[28.67,44.0435],[28.67,44.0553],[28.6611,44.0553]] -T 2021-08-15 -k 0 -d c:\data