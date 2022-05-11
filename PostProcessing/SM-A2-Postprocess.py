#!/usr/bin/env python3
# coding: utf-8

# # ForCoast SM-A2 Post-processing routines
# 
# This is a set of routines to prepare the ouptputs to users of SM-A2 ForCoast service module. 
# It should be executed after the preparation of Lagrangian tracks (Parcels), once per pollutant sources and per target.
# 
# Basically, the script acquire the particle tracks computed in the Processing steps. 
# For certain output time intervals (default each 12h), it computes the number of particles being inside the target. 
# This [number of particles per time] is divided is divided by the release rate to provide the fraction of the release reaching the farm.
# Note that this number may be larger than 100% in case of recirculation (yet, each particle is counted only once).
#
# In addition, particles are seprated for different class of Age, ie. their travel time between the source and the target.
# This used used to compute a relative Risk factor, based on a user input specifications of risk factors for each class of ages.
# The total Risk for time interval ti, is obtained as the weighted sum of fraction of release for the differant age classes. 
#
# So, $R(t) = \sum_{ a \in A} (p_a(t).r_a).
# --> If all Age classe have a risk factor of 100%, and are characterized by a reaching fraction of 50%, the total risk will be 50.
# --> If a single class age has a risk of 70, and has a reaching rate of 100%, the total risk will be 70.
#
# # Input Requirements
#   ------------------
# 
# The scripts expects command line arguments and uses info in ./usr/USERNAME/config/config.yaml (same as the processing script).
# 
# USAGE : SM-A2-Postprocess.py -y <username> -s <source> -c <sourcecount> -t <target> -k <targetcount> -d <datadir>
#
#   <username>      : user ID, eg. "eforie" or "galway", should match a ./usr/<username>/ directory.
#   <source>        : source point, as an array [lon, lat], those are listed in ./usr/<username>/config/sources.txt, when the script is called from run_py.sh
#   <sourcecount>   : interger index of current source, used in the frame of multiple source recursive calls.
#   <target>        : array of coordinates for the target, [[lon1,lat1],...,[lonN, latN]], those are listed in ./usr/<username>/config/targets.txt, when the script is called from run_py.sh
#   <sourcecount>   : interger index of current source, used in the frame of multiple source recursive calls.
#   <datadir>       : NOT USED currently, instead processing results are stored in ./usr/<username>/output/
#
# # Outputs
#   -------
# 
#  All production is directed towards : ./usr/<username>/output/target_<k>_source_<c>/.  This includes : 
#     TS_violin.png            : a time series violin plot with age distributions and relative fraction
#     TS_Risk.png              : a time serie of the total risk factor
#     TS_Risk_chart.png        : a reminder of the age-specific risk factor used in the computation of the total risk
#     AllTracks_Alarm<ti>.png  : Situation maps for each time index <ti>
#     Risk.csv                 : a csv file providing the total risk factor for each date. Usefull for further selection for dissemination.
# 
# ## Module load and Specific definition
#
#    ...

###############
# DEBUG FLAGS #
debug   = False  # print outputs
fast    = False  # skip map tiles, which takes time and requires a connection.
skipmap = False  # skip maps
###############


import netCDF4
import matplotlib.pyplot as plt
import matplotlib.cm as cm 
from matplotlib import cm, colors, colorbar
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import xarray as xr
import numpy as np
import numpy.ma as ma
from datetime import timedelta as delta
from datetime import datetime

from cartopy.io import shapereader
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
request = cimgt.OSM()

import os
import yaml
import csv
import sys, getopt
import math

from shapely.geometry.polygon import Polygon
   
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hy:s:c:t:k:d:",["yamlfile="])
except getopt.GetoptError:
    print(opts)
    print ('USAGE : SM-A2-Postprocess.py -y <username> -s <source> -c <sourcecount> -t <target> -k <targetcount> -d <datadir>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('USAGE : SM-A2-Postprocess.py -y <username>')
        sys.exit()
    elif opt in ("-y", "--yamlfile"):
        USER_YAML_FILE = arg
    elif opt in ("-s"):
        # AC 10112021 - Not sure this is a safe parsing method .. Same below.
        exec('source = '+arg)
    elif opt in ("-c"):
        sourcecount = int(arg)
    elif opt in ("-t"):
        exec('target = '+arg)
    elif opt in ("-k"):
        targetcount = int(arg)
    elif opt in ("-d"): # AC 10112021 shouldn't be needed ? 
        datadir = arg
USER_YAML_FILE = '../usr/' + USER_YAML_FILE + '/config/config.yaml'
print ('yaml file is', USER_YAML_FILE)

with open(USER_YAML_FILE) as f:
    data   = yaml.load(f, Loader=yaml.FullLoader)
    fname         = '../usr/'+data['username'] +'/output/test_'+str(sourcecount)+'.nc'
    figdir        = '../usr/'+data['username'] +'/output/target_' + str(targetcount)+'_source_' + str(sourcecount)+'/'
    coastlinefile = '../usr/'+data['username'] +'/config/coastlines/lines.shp'

## ## ## Upon service subsription, such files should be downloaded and stored  ## ## ##
# To Get high res coastlines.
# https://clouds.eos.ubc.ca/~phil/courses/atsc301/coursebuild/html/hires_map.html
# To generate new shapefile for other domains : 
# ogr2ogr -skipfailures -f "ESRI Shapefile"  -clipsrc -10.0 53.0 -8.0 54.0 galway_coastlines coastlines-split-4326

# Particles release per hour (! assuming one file per source at this stage !)
    ReleaseRate    = data['npart']*60/data['repeatdt'] 

# Age after which we don't care.
    thres_Age      = data['thres_Age'] # hours

# Time bins for post-processing
    outtimestep    = data['outtimestep'] # hours

# Make animations or Not
    BuildAnim      = data['BuildAnim']

# Thresolds for alarm (should be a curve)
    Age_Alarm      = data['Age_Alarm']
    Fraction_Alarm = data['Fraction_Alarm']

    poly1           = target 
    sourcepoint     = source
    
# To build map domains
    domainpercentiles   = data['domainpercentiles']
    expansionfactors_x  = data['expansionfactors_x']
    expansionfactors_y  = data['expansionfactors_y']
    extformap           = data['extformap']

# Age Classes for risk output
    agesc = data['agesc']
    riskc = data['riskc']

#To adapt when multiple farms are considered
    sourcenames = data['sourcenames']
    targetnames = data['targetnames']

    polys=[poly1]

if not os.path.exists(figdir):
    os.makedirs(figdir)

print('Start Post-Processing for :')
print('user : ' + data['username']) 
print(' for source %s at lat: %s and lon: %s '%(sourcecount, sourcepoint[1], sourcepoint[0]))      
print(' for farm   %s which corners are at :'%(targetcount))
for i,p in enumerate(poly1):
    print('%s - Lat : %s , Lon : %s'%(i,p[1],p[0]))
          
# ## Functions
# Function to prepare axis for map plots.
def start_axes(title, extent, fig=None, sp=None, fast=fast):
    if fig is None:
        fig = plt.figure(figsize=(13, 5))        
    if sp is None:
        ax = fig.add_axes([0.01, 0.01, 0.85, 0.9],projection=ccrs.PlateCarree())
    else:
        ax = fig.add_subplot(sp,projection=ccrs.PlateCarree())
            
    ax.set_extent(extent)
    #ax.gridlines()
    if not fast:
        ax.add_image(request, 10, interpolation='spline36', regrid_shape=2000)
    
    if not fast and coastlinefile is not None:
        shape_project=ccrs.PlateCarree()
        shp = shapereader.Reader(coastlinefile)
        for record, geometry in zip(shp.records(), shp.geometries()):
            ax.add_geometries([geometry], shape_project,facecolor="none",edgecolor='black',lw=1)
    else:
        ax.coastlines(resolution='10m')
    #ax.set_xlim(-6, 36.5), ax.set_ylim(30, 46)
    #ax.set_aspect("equal")
    ax.set_title(title)
    #ax.gridlines(xlocs=range(25,42,1), ylocs=range(40,48,1),draw_labels=True)#
    return ax

# Function to identify if particles are in a polygon
import matplotlib.path as mpltPath

def mplpoly(xin,yin,poly):
    """
    Check if points are inside a polygon.

    Parameters
    ----------
    xin,yin : arrays of floats
        points coordinates
    poly : array 
        aray of coordinates for the polygon corners

    Returns
    -------
    A vector of boolean, of same size as xin and zin

    """
    if (isinstance(xin,np.ndarray))&(len(xin.shape)>1):
        fromarray=True
        x = xin.flatten()
        y = yin.flatten()
    else:
        fromarray=False
        x=xin
        y=yin

    path = mpltPath.Path(poly)
    points=np.array([x, y]).T
    inside = path.contains_points(points)
    
    if fromarray:
        inside= np.reshape(inside, xin.shape)
    return(inside)

# Functions to automatically set map extents
## The rationale is to consider the distance between sources and target, and to use it a lentgh scale.

def distance(origin, destination):
    """
    Calculate the Haversine distance.

    Parameters
    ----------
    origin : tuple of float
        (lat, long)
    destination : tuple of float
        (lat, long)

    Returns
    -------
    distance_in_km : float
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6378.1  # km

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) * math.sin(dlat / 2) +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) * math.sin(dlon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius * c
    return d

def destination(origin,brng,d):
    """
    Calculate destination given orign, bearing and distance

    Parameters
    ----------
    origin : tuple of float
        (lat, long)
    brng : float  
        Bearing in degrees, clockwise from north
    d : float
        ditsance, in km
    Returns
    -------
    destination : tuple of float
        (lat, long)
    """
    R = 6378.1 #Radius of the Earth

    brng = math.radians(brng)
    lat1 = math.radians(origin[0]) #Current lat point converted to radians
    lon1 = math.radians(origin[1]) #Current long point converted to radians

    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) + math.cos(lat1)*math.sin(d/R)*math.cos(brng))
    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)

    return((lat2,lon2))

def update_extent_sourcetarget(sourcepoint, poly, extfactor=1.5):
    """
    Compute the map extent based on a point (source location) and a polygon (target)

    A typical length scale is computed as the distance between the source point and the center of the polygon.
    The map extension is obtained by adding expfactor*lengthscale to the boundaries of the smaller rectangle including the source and the target.

    Parameters
    ----------
    sourcepoint : array 
       [lon,lat]
    poly : array
       [[lon1,lat1], ...]
    extfactor : 
       expansion factor for the length scale

    Returns
    -------
    domain extension : 
       [ext_lonmin, ext_lonmax, ext_latmin, ext_latmax]
    """
    #get the corners of the smallest rectangle surrounding source an target
    lonmax = np.array( [sourcepoint[0], np.array([ p[0] for p in poly]).max()]).max()
    lonmin = np.array( [sourcepoint[0], np.array([ p[0] for p in poly]).min()]).min()
    latmax = np.array( [sourcepoint[1], np.array([ p[1] for p in poly]).max()]).max()
    latmin = np.array( [sourcepoint[1], np.array([ p[1] for p in poly]).min()]).min()

    # get distance between targets (mean value of all points) and source
    targetmeanlon = np.array([ p[0] for p in poly]).mean()
    targetmeanlat = np.array([ p[1] for p in poly]).mean()
    d             = distance( (targetmeanlat,targetmeanlon) , (sourcepoint[1],sourcepoint[0]))

    ext_lonmin = destination(  (latmin,lonmin), 270.0, d*extfactor)[1]
    ext_lonmax = destination(  (latmax,lonmax), 90.0, d*extfactor)[1]

    ext_latmin = destination(  (latmin,lonmin), 180.0, d*extfactor)[0]
    ext_latmax = destination(  (latmax,lonmax), 0.0, d*extfactor)[0]

    return([ext_lonmin, ext_lonmax, ext_latmin, ext_latmax])

def scale_bar(ax, length=None, location=(0.5, 0.05), linewidth=3):
    """
    Add a kilometric scale bar to cartopy plots. 

    ax(axes)         : is the axes to draw the scalebar on.
    length(float)    : is the length of the scalebar in km.
    location(tuple)  : is center of the scalebar in axis coordinates (ie. 0.5 is the middle of the plot)
    linewidth(float) : is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom')


#######################################
### The Post-Processing starts here ###
#######################################

# ## Data Load
data_xarray = xr.open_dataset(fname)
np.set_printoptions(linewidth=160)

ns_per_hour = np.timedelta64(1, 'h') # nanoseconds in an hour
outputdt = delta(hours=outtimestep)

timerange = np.arange(np.nanmin(data_xarray['time'].values),
                      np.nanmax(data_xarray['time'].values)+np.timedelta64(outputdt), 
                      outputdt) # timerange in nanoseconds


###
##  Time Series
# Compute the number of particles entering the polygon with an age below `thres_Age`.
###

labels=[]; labelsm=[]; ages=[]; ninside=[]

if debug:
    print('timerange')
    print(timerange)

ninside =np.zeros((len(timerange)))

for i,t in enumerate(timerange):
    tinds       = (data_xarray['time'] > t)  & (data_xarray['time'] <= t+np.timedelta64(outputdt))
    points_x    = data_xarray['lon'].values[tinds]
    points_y    = data_xarray['lat'].values[tinds]
    areinside   = mplpoly(points_x,points_y, poly1)
    ageofinside = data_xarray['age'].values[tinds][areinside]/86400
    # Each individual particle should be counted only once !
    bid=tinds.values.copy()
    bid[bid==True]=areinside
    ninside[i]=bid.any(axis=1).sum()
    ages.append(ageofinside[~np.isnan(ageofinside)])
    # AC 11.11.2021 - this is boring date parsing. A very handy way out would be to use pandas bu that seems overkill to inlcude it as a project requirement.
    dstring = np.datetime_as_string(t, unit='m')
    y=dstring[2:4] ; M=dstring[5:7] ; d=dstring[8:10] ; h=dstring[11:13] ; m=dstring[14:16] 
    labels.append(d+'/'+M+'\n'+h+':'+m) # used as labels
    labelsm.append(d+'/'+M+' '+h+':'+m) # used on maps
    #labels.append(t.astype(datetime).strftime('%d%b -%hh'))

# Fraction of particles released during the same time interval
relninside = ninside/(ReleaseRate*outtimestep)*100

if debug:
    print('len(ages)')
    print(len(ages))

############
## Violin plot
############

# dealing with empty cases
nans = [float('nan'), float('nan')]
agesplot       = [ a if a.any() else nans for a in ages ] #if len(a)>0 ]

relninsideplot = relninside# [ r for r,a in zip(relninside, ages) if len(a)>0 ]

fig = plt.figure(figsize=(10,4))
ax=fig.add_axes([0.05,0.1,.85,.8], title = 'Influence of '+sourcenames[sourcecount] )

vpl = ax.violinplot(agesplot,showmedians=False, showextrema=False)

cmap = cm.get_cmap('RdPu')
ccolors = [cmap(r/100) for r in relninsideplot]

if debug:
    for c, r in zip(ccolors, relninsideplot):
        print(' r:%s, r/100:%s gives c: %s'%(r,r/100,c))

# Color to the violin elements according to fraction
for patch, color in zip(vpl['bodies'], ccolors):
    patch.set_color(color)
    patch.set_alpha(.8)

ax.xaxis.set_tick_params(direction='out')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(np.arange(1, len( agesplot ) + 1))
ax.set_xticklabels(labels, fontsize=8)
ax.set_xlim(0.25, len(agesplot) + 0.75)
#plt.xticks(rotation=45)
ax.set_ylabel('Time to reach ' + targetnames[targetcount] +' [days]', labelpad = 2)

# colorbar
cax=fig.add_axes([0.91,0.1,.03,.8])
norm = colors.Normalize(vmin=0, vmax=100)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
cax.set_ylabel('Fraction reaching ' + targetnames[targetcount] +' [%]', labelpad = 2)

fig.savefig( figdir+'TS_violin.png', dpi=200)

## Violin Plots for all time frames.

for ti,t in enumerate(timerange): 
    plt.close()
    fig = plt.figure(figsize=(10,4))
    ax=fig.add_axes([0.05,0.1,.85,.8], title = 'Influence of '+sourcenames[sourcecount] )
    vpl = ax.violinplot(agesplot,showmedians=False, showextrema=False)
    cmap = cm.get_cmap('RdPu')
    #    ccolors = [cmap(r/100) for r in relninsideplot]
    for pi, (patch, color) in enumerate(zip(vpl['bodies'], ccolors)):
        patch.set_color(color)
        patch.set_alpha(.8)
        if pi==ti:
            patch.set_edgecolor('black')

    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len( agesplot ) + 1))
    ax.set_xticklabels(labels, fontsize=8)

    ax.get_xticklabels()[ti].set_fontsize(10)
    ax.get_xticklabels()[ti].set_weight("bold")


    ax.set_xlim(0.25, len(agesplot) + 0.75)
    #plt.xticks(rotation=45)
    ax.set_ylabel('Time to reach ' + targetnames[targetcount] +' [days]', labelpad = 2)

    # colorbar
    cax=fig.add_axes([0.91,0.1,.03,.8])
    norm = colors.Normalize(vmin=0, vmax=100)
    cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='vertical')
    cax.set_ylabel('Fraction reaching ' + targetnames[targetcount] +' [%]', labelpad = 2)

    fig.savefig( figdir+'TS_violin_'+'%03d'%(ti)+'.png', dpi=200)

############
## Maps
# For now making maps for all time frames. 
############

if not skipmap : 

    for ti,t in enumerate(timerange): 
        plt.close()
        fig = plt.figure(figsize=(4,4))
        ax1   = start_axes(title='Map1', extent=update_extent_sourcetarget(sourcepoint, poly1,.8), fig=fig)
        tinds = (data_xarray['time'] > t)  & (data_xarray['time'] <= t+np.timedelta64(outputdt))
        sc1= ax1.scatter(data_xarray['lon'].values[tinds],
                        data_xarray['lat'].values[tinds], 5,
                        data_xarray['age'].values[tinds]/86400, vmin=0,vmax=5, alpha=.5, cmap='viridis_r')
        if polys[0] is not None:
            pol = Polygon(poly1)
            ax1.add_geometries([pol], facecolor='orange', edgecolor='red', alpha=0.8, crs=ccrs.PlateCarree())

        ax1.text( poly1[0][0],poly1[0][1],  targetnames[targetcount])

        ax1.scatter(sourcepoint[0],sourcepoint[1], 25,  'red')
        ax1.text( sourcepoint[0],sourcepoint[1],  sourcenames[sourcecount])

        cax = fig.add_axes([0.88,0.1,0.04,0.8])
        norm = colors.Normalize(vmin=0, vmax=5)
        cb = colorbar.ColorbarBase(cax, cmap=cm.get_cmap('viridis_r'), norm=norm, orientation='vertical')
        cax.set_title('Age [d]')

        ax1.set_title( labelsm[ti])
        scale_bar(ax1, 1)
        plt.savefig(figdir+'AllTracks_Alarm_'+'%03d'%(ti)+'.png', dpi=200)

############
# Risk
############ 

ninsideperclass =np.zeros((len(timerange), len(agesc)-1))

for i,t in enumerate(timerange):
    tinds       = (data_xarray['time'] > t)  & (data_xarray['time'] <= t+np.timedelta64(outputdt))
    points_x    = data_xarray['lon'].values[tinds]
    points_y    = data_xarray['lat'].values[tinds]
    points_age  = data_xarray['age'].values[tinds]/86400
    
    for j in range(len(agesc)-1):
        age_id = (points_age >= agesc[j]) & (points_age < agesc[j+1])
        points_xl  = points_x[age_id]
        points_yl  = points_y[age_id]
        areinside = mplpoly(points_xl,points_yl, poly1)

        # Each individual particle should be counted only once !
        bid=tinds.values.copy()
        bid2 = age_id.copy()
        bid2[bid2==True]=areinside             
        bid[bid==True]=bid2
        ninsideperclass[i,j]=bid.any(axis=1).sum()

relninsideperclass = ninsideperclass/(ReleaseRate*outtimestep)*100

# FIXME !! sum of relative fraction per classes isn't stricly equal to relative fractions  ????
#    Shall we renormalize as a quick patch ? 

'''
print('relninsideperclass')
print(relninsideperclass)

print('relninsideperclass.sum(axis=1)')
#print(relninsideperclass.shape)
print(relninsideperclass.sum(axis=1))
'''

## Patching the above ##
relninsideperclass_ref = relninsideperclass.copy()
for i,t in enumerate(timerange):
    for j in range(len(agesc)-1):
        relninsideperclass_ref[i,j]=relninsideperclass[i,j] * relninside[i]/relninsideperclass[i,:].sum() if relninsideperclass[i,:].sum()!=0.0 else 0.0

relninsideperclass = relninsideperclass_ref.copy()
##

'''
print('relninside')
#print(relninside.shape)
print(relninside)

print('relninsideperclass.sum(axis=1)')
#print(relninsideperclass.shape)
print(relninsideperclass.sum(axis=1))

print('relninsideperclass')
print(relninsideperclass)
'''

risk = np.zeros((len(timerange)))

for ti,t in enumerate(timerange):
    risk[ti] = (relninsideperclass[ti]*riskc).sum()/100

# Risk Figure
plt.close()
fig = plt.figure(figsize=(10,2))
ax = fig.add_axes([0.05,0.3,.85,.5])

ax.imshow( np.tile(risk, [2,1] ), cmap='RdYlGn_r', vmin=0, vmax=100)

ax.xaxis.set_tick_params(direction='out')
ax.xaxis.set_ticks_position('bottom')
ax.set_xticks(np.arange(0, len( agesplot )))
ax.set_xticklabels(labels, fontsize=8)#[ l for i,l in enumerate(labels) if len(ages[i])>0])
ax.set_yticks([])
ax.set_yticklabels([])
ax.set_ylabel('Risk')
ax.spines['top'].set_color('white') 
ax.spines['left'].set_color('white') 
ax.spines['right'].set_color('white') 
ax.spines['bottom'].set_color('white') 

fig.savefig( figdir+'TS_Risk.png', dpi=200)

# Risk Figure, per time frame
for ti,t in enumerate(timerange):
    plt.close()
    fig = plt.figure(figsize=(10,2))
    ax = fig.add_axes([0.05,0.3,.85,.5])

    ax.imshow( np.tile(risk, [2,1] ), cmap='RdYlGn_r', vmin=0, vmax=100)

    ax.xaxis.set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(0, len( agesplot )))
    ax.set_xticklabels(labels, fontsize=8)#[ l for i,l in enumerate(labels) if len(ages[i])>0])
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_ylabel('Risk')
    ax.spines['top'].set_color('white') 
    ax.spines['left'].set_color('white') 
    ax.spines['right'].set_color('white') 
    ax.spines['bottom'].set_color('white') 

    ax.get_xticklabels()[ti].set_fontsize(10)
    ax.get_xticklabels()[ti].set_weight("bold")

    fig.savefig( figdir+'TS_Risk_'+'%03d'%(ti)+'.png', dpi=200)

# Risk Panel figure
plt.close()
fig = plt.figure(figsize=(4,2))
ax = fig.add_axes([0.14,0.3,.85,.5])

agesc_mid = (np.array(agesc)[1:]+np.array(agesc)[:-1] )/2
agesc_int = (np.array(agesc)[1:]-np.array(agesc)[:-1] )

cmap = cm.get_cmap('RdYlGn_r')
ccolors = [cmap(r/100) for r in riskc]

ax.bar( agesc_mid, np.array(riskc), width=agesc_int , color =  ccolors )

ax.set_xticks(agesc)
ax.set_ylim([0,100])
ax.set_yticks([0,50,100])
ax.set_yticklabels(['Null','Mid','Strong'])
ax.set_title('Risk scale')
ax.set_xlabel('Time to reach - [days]')

fig.savefig( figdir+'TS_Risk_chart.png', dpi=200)

# Writing Risk values for further use
with open(figdir+'Risk.csv', 'w', newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerows(zip(timerange,risk))

############
# For animation I'd rather suggest to use imagemagick.. producing all png would make it easier for further diffusion, in any case (gif, or display with scrollable timebar..)
#
# The idea would be to complete the violin and risk time series with a 'current' pointer. Then processing all images, all buileins and collecting those in such a slide show.
#
'''
if not BuildAnim: 
    raise SystemExit("You choosed to skip the animation.")
else:
    if polys[0] is not None:
        apoly1=np.array(poly1)
        ml=np.min(apoly1[:,1])
        Ml=np.max(apoly1[:,1])
        mL=np.min(apoly1[:,0])
        ML=np.max(apoly1[:,0])
        extp=[mL-(ML-mL), ML+(ML-mL), ml-(Ml-ml), Ml+(Ml-ml)]


    from matplotlib.animation import FuncAnimation

    fig = plt.figure(figsize=(15,5))
    nexts=len(exts)

    axs = []
    scats=[]
    scatpolys=[]

    time_id = np.where(data_xarray['time'] == timerange[0]) # Indices of the data where time = 0

    for i,ext in enumerate(exts):     
        axs.append(start_axes('Zoom', fig=fig,sp=int('1'+str(nexts)+str(i+1)), extent=ext))
        

        scats.append( axs[i].scatter(data_xarray['lon'].values[time_id],
                                data_xarray['lat'].values[time_id],10,
                                data_xarray['age'].values[time_id]/86400, vmin=0,vmax=5))
        if polys[0] is not None:
            scatpolys.append(axs[i].scatter(np.array(poly1)[:,0],np.array(poly1)[:,1] ))

    clb3 =plt.colorbar(scats[i])
    clb3.ax.set_xlabel('Age [d]')

    t = np.datetime_as_string(timerange[0], unit='m')
    title = axs[0].set_title('Particles at t = '+t)

    def animate(i):
        t = np.datetime_as_string(timerange[i], unit='m')
        title.set_text('Particles at t = '+t)
        
        time_id = np.where((data_xarray['time'] >= timerange[i]) & (data_xarray['time'] < timerange[i+1]))
        
        
        for i,ext in enumerate(exts):  
            scats[i].set_offsets(np.c_[data_xarray['lon'].values[time_id], data_xarray['lat'].values[time_id]])
            scats[i].set_array(data_xarray['age'][time_id]/86400)
        
    anim = FuncAnimation(fig, animate, frames = len(timerange)-1, interval=500)

    from IPython.display import HTML
    HTML(anim.to_jshtml())
    anim.save('GAL1.mp4', fps=5, extra_args=['-vcodec', 'libx264'])
'''


