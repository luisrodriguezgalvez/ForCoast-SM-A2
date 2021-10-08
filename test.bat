# Download input (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PreProcessing
python forcoast_download_yml.py -a eforie -T 2021-09-10 -p 20 -d c:\data
cd ..
# Run model (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd  Processing
python forcoast.py -y eforie -T 2021-09-20 -p 10 -s [28.6803,44.2201] -t [28.6611,44.0535,28.6617,44.0553] -d c:\data 
cd ..
# Run postprocessing (uses data directory from command line, such that we don't need to change the yaml file for testing)
cd PostProcessing
python SM-A2-Postprocess.py -y eforie.yaml -s [28.6803,44.2201] -t [28.6611,44.0535,28.6617,44.0553] -d c:\data 
cd ..