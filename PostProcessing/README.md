# ForCoast SM-A2 Post-processing routines
 
## Description
This is a set of routines to prepare the ouptputs to users of SM-A2 ForCoast service module. 
It should be executed after the preparation of Lagrangian tracks (Parcels), once per pollutant sources and per target.
 
Basically, the script acquire the particle tracks computed in the Processing steps. 
For certain output time intervals (default each 12h), it computes the number of particles being inside the target. 
This [number of particles per time] is divided is divided by the release rate to provide the fraction of the release reaching the farm.
Note that this number may be larger than 100% in case of recirculation (yet, each particle is counted only once).

In addition, particles are seprated for different class of Age, ie. their travel time between the source and the target.
This used used to compute a relative Risk factor, based on a user input specifications of risk factors for each class of ages.
The total Risk for time interval ti, is obtained as the weighted sum of fraction of release for the differant age classes. 

So, $R(t) = \sum_{ a \in A} (p_a(t).r_a)$.
--> If all Age classe have a risk factor of 100%, and are characterized by a reaching fraction of 50%, the total risk will be 50.
--> If a single class age has a risk of 70, and has a reaching rate of 100%, the total risk will be 70.

## Input Requirements
 
The scripts expects command line arguments and uses info in ./usr/USERNAME/config/config.yaml (same as the processing script).
 
USAGE : SM-A2-Postprocess.py -y <username> -s <source> -c <sourcecount> -t <target> -k <targetcount> -d <datadir>

   <username>      : user ID, eg. "eforie" or "galway", should match a ./usr/<username>/ directory.
   <source>        : source point, as an array [lon, lat], those are listed in ./usr/<username>/config/sources.txt, when the script is called from run_py.sh
   <sourcecount>   : interger index of current source, used in the frame of multiple source recursive calls.
   <target>        : array of coordinates for the target, [[lon1,lat1],...,[lonN, latN]], those are listed in ./usr/<username>/config/targets.txt, when the script is called from run_py.sh
   <sourcecount>   : interger index of current source, used in the frame of multiple source recursive calls.
   <datadir>       : NOT USED currently, instead processing results are stored in ./usr/<username>/output/


## Outputs
 
All production is directed towards : ./usr/<username>/output/target_<k>_source_<c>/.  This includes : 
    TS_violin.png            : a time series violin plot with age distributions and relative fraction
    TS_Risk.png              : a time serie of the total risk factor
    TS_Risk_chart.png        : a reminder of the age-specific risk factor used in the computation of the total risk
    AllTracks_Alarm<ti>.png  : Situation maps for each time index <ti>
    Risk.csv                 : a csv file providing the total risk factor for each date. Usefull for further selection for dissemination.

## Module load and Specific definition

    ...
