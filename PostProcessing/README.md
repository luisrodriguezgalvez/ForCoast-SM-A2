Use SM-A2-Postprocess.py to generate output figures. 
USAGE : 

    python3 SM-A2-Postprocess.py -y USER_YAML_FILE

Two example of yaml file are given :

* SMA3-PostProcess-Eforie.yaml 
* SMA3-PostProcess-Galway.yaml

Required input files includes netcdf for particle tracks, and detailled coastlines. Those cannot be included in github.

The codes generates output files for display. 

* TS1.png : gives the mean age of particle reaching the sites and proportion of release reaching the site.
* TS2.png : gives the proportion of release reaching the site, for different age fractions.
* TS2_Alarm.png : highlights Alarms (on the basis of user-defined thresholds).
* AllTracks_AlarmX.png : A series of figure for each particular alarm : 

In addtion, each day of alarm is listed in a file _AlarmsList.csv_ , which list the dates for alarm and the name of the corresponding figure. 
