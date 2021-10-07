cd PreProcessing
# python forcoast_download_yml.py -T 2021-09-10 -p 20
cd ..
cd  Processing
python forcoast.py eforie 2021-09-20 10 c:\data [28.6803,44.2201] [28.6611,44.0535,28.6617,44.0553]
cd ..
cd PostProcessing
python SM-A2-Postprocess.py -y eforie.yaml -d c:\data -s [28.6803,44.2201] -t [28.6611,44.0535,28.6617,44.0553]
cd ..