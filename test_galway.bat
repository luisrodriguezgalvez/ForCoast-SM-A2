:: Download input (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PreProcessing
python forcoast_download_yml.py -a galway -T 2021-11-01 -p 5 -d c:\data
cd ..

:: Run model (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd Processing
python forcoast.py -y galway -T 2021-11-01 -p 5 -s [-8.927046,53.142875,-1.0] -d c:\data\
cd ..

:: TO BE MODIFIED FOR FLEXIBILITY Change name of test.nc to test_1.nc (1 is because the sourcecount -c in this case in the Postprocessing will be 1)
cd usr\galway\output
ren test.nc test_0.nc
cd ..\..\..

:: Run postprocessing (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PostProcessing
python SM-A2-Postprocess.py -y galway -s [-8.927046,53.142875,-1.0] -c 0 -t [[-9.0175,53.16816958],[-9.041166667,53.16817019],[-9.042,53.18817033],[-9.018833333,53.18816981]] -k 0 -d c:\data 
cd ..

:: Run bulletin script
cd BulletinScript
python bulletin_script.py -y galway -s [-8.927046,53.142875,-1.0] -c 0 -t [[-9.0175,53.16816958],[-9.041166667,53.16817019],[-9.042,53.18817033],[-9.018833333,53.18816981]] -T 2021-09-20 -k 0 -d c:\data
cd ..
