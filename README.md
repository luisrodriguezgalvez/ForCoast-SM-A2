# ForCoast-SM-A2
ForCoast Service Module A2

Contains all required to : 
  * acquire data
  * compute tracks
  * derives stats
  * combine bulletin
   
   
For a user `<username>`, there should be `./usr/<username>/config/` with : 
  * config.yaml : provding instructions for both the processing and post-processing steps.
  * sources.txt : one row for each sources coordinates. Names for these sources should be given in the config.yaml.
  * targets.txt : one row for each targets coordinates. Names for these targets should be given in the config.yaml.
  * coastline/line.shp, containing coastline information. 

Final products are to be found in `./usr/<username>/output/`
  
USAGE : `run_py.sh <username> <start date> <duration>`
