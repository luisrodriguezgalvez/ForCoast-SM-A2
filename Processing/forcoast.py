#!/data/anaconda3/envs/py3_parcels/bin/python

# ForCoast ServiceModule A2 -- Ocean Parcels Lagrangian tracking

from parcels import FieldSet, Field, NestedField, VectorField
from parcels import ParticleSet, ParticleFile, Variable, ErrorCode
from parcels import JITParticle #, ScipyParticle
from parcels import plotTrajectoriesFile
from parcels.kernel import Kernel
from glob import glob
import numpy as np
import sys
from datetime import timedelta as delta
import time
from forcoast_loaders import *
from forcoast_kernels import *
import yaml
import os
from pathlib import Path
import sys, getopt

argv = sys.argv[1:]

try:
    opts, args = getopt.getopt(argv,"hy:T:p:d:s:",["yamlfile="])
except getopt.GetoptError:
    print ('forcoast.py -y <yamlfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('forcoast.py -y <yamlfile>')
        sys.exit()
    elif opt in ("-y", "--yamlfile"):
        USER_YAML_FILE = arg
    elif opt in ("-T"):
        startdate = arg
    elif opt in ("-p"):
        period = arg 
    elif opt in ("-d"):
        datadir = arg
    elif opt in ("-s"):
        source = arg
#    elif opt in ("-t"):
#        target = arg

config_file = '../usr/' + USER_YAML_FILE + '/config/config.yaml'

# OPTIONS
with open(config_file) as f:
  options=yaml.load(f,Loader=yaml.Loader)
npart=options['npart']
if options['run3D']:
    print('running in 3D')
else:
    print('running in 2D')

# Replace options from yml file with options passed from command line
if not argv:
    print("Use input from yml file")
else:
    print("Replace selected input with input from command line")

    if "startdate" in locals():
        options['sdate'] = startdate
    if "period" in locals():
        options['simlength'] = int(period)
    # options['PHY_path'] = argv[5]
    if "source" in locals():
        # Parse pollutant release coordinate from command line
        bbox_input_temp = str(source)[1:-2]
        bbox_input_temp = bbox_input_temp.split(',')
        bbox_input_x0 = [float(bbox_input_temp[0])]
        bbox_input_y0 = [float(bbox_input_temp[1])]
        bbox_input_z0 = [float(bbox_input_temp[2])]
        options['x0'] = bbox_input_x0
        options['y0'] = bbox_input_y0
        options['z0'] = bbox_input_z0
#    if "target" in locals():
#        # Parse target area coordinates from command line
#        bbox_output_temp = str(target)[1:-1]
#        bbox_output_temp = bbox_output_temp.split(',')
#        bbox_output_left = float(bbox_output_temp[0])
#        bbox_output_bottom = float(bbox_output_temp[1])
#        bbox_output_right = float(bbox_output_temp[2])
#        bbox_output_top = float(bbox_output_temp[3])
#        bbox_output_targetx = [bbox_output_left,bbox_output_right]
#        bbox_output_targety = [bbox_output_bottom,bbox_output_top]
#        options['targetx'] = bbox_output_targetx
#        options['targety'] = bbox_output_targety
    if "datadir" in locals():
        options['PHY_path'] = datadir

print("Pollutant release coordinates: x0=" + str(options['x0']) + ', y0=' + str(options['y0']))
#print("Target area coordinates: targetx=" + str(options['targetx']) + ', targety=' + str(options['targety']))
print("Start date: " + options['sdate'])
print("Simulation length: " + str(options['simlength']))
print("Data directory: " + options['PHY_path'])

# LOADING EULERIAN VELOCITY FIELD
if options['PHY_type']=='ROMS':
  romsfiles=sorted(glob(options['PHY_path']+options['files']))
  print('romsfiles :' + romsfiles[0])
  fieldset=get_roms_fields(romsfiles,run3D=options['run3D'],chunksize=False,vdiffusion=options['vdiffusion'],beaching=options['beaching'])
elif options['PHY_type']=='NEMO':
  #print('PHY PATH:', options['PHY_path'])
  nbgrids=len(options['mfiles'])
  myfieldset=[]
  for g in range(nbgrids):
    print('PHY PATH:', options['PHY_path']+options['wfiles'][g])
    ufiles = sorted(glob(options['PHY_path']+options['ufiles'][g]))
    vfiles = sorted(glob(options['PHY_path']+options['vfiles'][g]))
    wfiles = sorted(glob(options['PHY_path']+options['wfiles'][g]))
    mesh_mask = options['PHY_path'] + options['mfiles'][g]
    indices = {'lon': range(options['range_i1'][g],options['range_i2'][g]), 'lat': range(options['range_j1'][g],options['range_j2'][g])} # NEMO puts zero along the (ghost) boundaries
    if options['range_i2'][g]>0:
      myfieldset.append(get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,run3D=options['run3D'],indices=indices,vdiffusion=options['vdiffusion'],beaching=options['beaching']))
    else:
      myfieldset.append(get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,run3D=options['run3D']                ,vdiffusion=options['vdiffusion'],beaching=options['beaching']))
      
  if options['nesting']==True:
    print('using ',nbgrids,' nested grids')

    # NestedFields
    ufields=[];vfields=[];wfields=[];Kzfields=[];KzEVDfields=[];tmaskfields=[];
    for g in range(nbgrids):
      ufields.append(myfieldset[g].U); vfields.append(myfieldset[g].V)
      if options['run3D']:
        wfields.append(myfieldset[g].W)
        if options['vdiffusion']:
          Kzfields.append(myfieldset[g].Kz) ; KzEVDfields.append(myfieldset[g].Kz_EVD);
      if options['beaching']>0:
        tmaskfields.append(myfieldset[g].tmask)
    U=NestedField('U', ufields)
    V=NestedField('V', vfields)
    if options['run3D']:
        W=NestedField('W', wfields)
        if options['vdiffusion']:
            Kz=NestedField('Kz', Kzfields)
            Kz_EVD=NestedField('Kz_EVD', KzEVDfields)
    if options['beaching']>0:
        tmask=NestedField('tmask', tmaskfields)
    
    if options['run3D']:
        if options['vdiffusion']:
            if options['beaching']>0:
                fieldset = FieldSet(U, V, {'W':W, 'Kz':Kz, 'Kz_EVD':Kz_EVD, 'tmask':tmask})
            else:
                fieldset = FieldSet(U, V, {'W':W, 'Kz':Kz, 'Kz_EVD':Kz_EVD})
            fieldset.add_constant('dres', 0.01)
        else:
            if options['beaching']>0:
                fieldset = FieldSet(U, V, {'W':W, 'tmask':tmask})
            else:
                fieldset = FieldSet(U, V, {'W':W})
    else:
        if options['beaching']>0:
            fieldset = FieldSet(U, V, {'tmask':tmask})
        else:
            fieldset = FieldSet(U, V)

  else:
    print('No nesting, using just 1 grid')
    fieldset=myfieldset[0]
    del myfieldset
else:
  print('unrecognized physical model, aborting...')
  sys.exit(1)

# LOADING WAVE FIELD
if options['stokes']==True:
    print('using currents and Stokes drift')
    mesh_mask  = options['WAV_path'] + options['wavemesh']
    get_wav_fields(fieldset,options['WAV_path'],options['run3D'],mesh_mask,chunksize=False)
else:
    print('using just currents, no Stokes drift')
fieldset.add_constant('maxage', options['maxage'])

# CREATE PARTICLE SET
if options['hdiffusion']:
    fieldset.add_constant('Kx',options['Kx'])
    fieldset.add_constant('Ky',options['Ky'])

class ForCoastParticle(JITParticle):
    age =  Variable('age', dtype=np.int32,   initial=0)
    beached = Variable('beached', dtype=np.int32,   initial=0.)   # 0: sea,     1: keep frozen on land    2: push beached particles back in sea
    if options['beaching']>0:                                        # prev* is used for pushing back in case (beaching==2)  or (beaching>=1 and crossdike detection)
        prevlon = Variable('prevlon', dtype=np.float64, initial=0., to_write=False)
        prevlat = Variable('prevlat', dtype=np.float64, initial=0., to_write=False)
        prevdep = Variable('prevdep', dtype=np.float64, initial=0., to_write=False)
        
if options['hdiffusion']==False and options['vdiffusion']==False:
    npart=1
lon0=np.zeros(npart*len(options['x0']));lat0=np.zeros(npart*len(options['x0']));depth0=np.zeros(npart*len(options['x0']));
for i in range(0,len(options['x0'])):
    lon0[i*npart:(i+1)*npart]   = options['x0'][i]
    lat0[i*npart:(i+1)*npart]   = options['y0'][i]
    depth0[i*npart:(i+1)*npart] = options['z0'][i]
    
repeatdt=delta(minutes=options['repeatdt']) # release particles every...
pset = ParticleSet.from_list(fieldset=fieldset, pclass=ForCoastParticle,
                             lon=lon0, lat=lat0, depth=depth0,
                             time=np.datetime64(options['sdate']),
                             repeatdt=repeatdt,
                             lonlatdepth_dtype=np.float64 )

# CHOOSE KERNELS
if options['stokes']:
    kernels = pset.Kernel(Dbl_AdvectionRK4_3D_clumsy) if options['run3D'] else pset.Kernel(Dbl_AdvectionRK4)
else:
    kernels = pset.Kernel(AdvectionRK4_3D) if options['run3D'] else pset.Kernel(AdvectionRK4, delete_cfiles=False)

if options['hdiffusion']:
    print('Using horizontal diffusion')
    kernels += DiffusionUniformKh

if options['run3D'] and options['vdiffusion']:
    print('Using vertical diffusion')
    kernels += DiffusionZ

# prevent beaching
if options['beaching']==1:
    print('Freezing beached particles')
    kernels += Frozenbeach
    if options['experiment']=="Eforie":
       c_includefile = path.join('parcels/c_kernels/crossdike1.h') # path.dirname(__file__)
elif options['beaching']==2:
    print('Un-beaching beached particles')
    if (options['run3D']):
       if options['PHY_type']=="NEMO":
          kernels += Unbeaching3D_nemo
       elif options['PHY_type']=='ROMS':
          kernels += Unbeaching3D_roms
    kernels += Unbeaching2D
    if options['experiment']=="Eforie":
       c_includefile = path.join('parcels/c_kernels/crossdike2.h') # path.dirname(__file__),
# prevent cross-dike
if options['experiment']=="Eforie" and options['beaching']>0:
    with open(c_includefile,'r') as f:
        c_include=f.read()
    kernels += pset.Kernel(ckernel_dikes, c_include=c_include, delete_cfiles=False)
if options['beaching']>0:
    kernels += save_position
else:
    print('Beaching: no action')
    

# E.Coli decay
kernels += decay
    
# RUN SIMULATION    
output_file = pset.ParticleFile(name=options['out_filename'], outputdt=delta(minutes=options['savedt']))
t0=time.time()
if options['through_surface']:
  pset.execute(kernels, runtime=delta(days=options['simlength']), dt=delta(minutes=options['dt']), output_file=output_file, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle1, ErrorCode.ErrorThroughSurface: DeleteParticle2})
else:
  pset.execute(kernels, runtime=delta(days=options['simlength']), dt=delta(minutes=options['dt']), output_file=output_file, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle,  ErrorCode.ErrorThroughSurface: DeleteParticle})  
              # verbose_progress=False)
t1=time.time() ; totaltime=t1-t0 ; print("computing time :",totaltime)
output_file.close()


# POST-TREATMENT DETECTION : this is de-activated, move to Arthur's post-processing code
'''
if options['detection']:
    import xarray as xr
    print("Particle detection in target area")
    events=0
    out = xr.open_dataset(options['out_filename'])
    timecount=np.size(out.lon.values,0)
    partcount=np.size(out.lon.values,1)
    for p in range(0,partcount-1):
        for t in range(0,timecount-1):
            if options['targetx'][0]<=out.lon.values[t,p]<=options['targetx'][1] and options['targety'][0]<=out.lat.values[t,p]<=options['targety'][1]:
                print('target attained ',t,p)


                
## PLOT PARCELS RESULTS as simple diagnostic
import xarray as xr
import matplotlib.pyplot as plt
data_xarray = xr.open_dataset(output_file.name)  # trajectories, time, lon, lat, z, prevlon, prevlat, prevdep, beached [obs=timestep, traj=particle]
if options['experiment']=='Eforie':
  x=data_xarray['lon'].values ; y=data_xarray['lat'].values ;
  fig,ax=plt.subplots(1,3)
  for sp in range(3):
    dx=fieldset.gridset.grids[sp].lon[1]-fieldset.gridset.grids[sp].lon[0] ; dy=fieldset.gridset.grids[sp].lat[1]-fieldset.gridset.grids[sp].lat[0]
    ax[sp].pcolormesh(fieldset.gridset.grids[sp].lon-(dx/2), fieldset.gridset.grids[sp].lat-(dy/2), fieldset.tmask[sp].data[0,0,:,:], shading='auto')
    ax[sp].plot(x.T,y.T)
    #plot dikes
    ax[sp].plot(np.array([28.671237945556641,28.673706054687500,28.673706054687500,28.676176071166992,28.681114196777344,28.681114196777344,28.700866699218750,28.700866699218750,28.703336715698242,28.703336715698242,28.705804824829102,28.705804824829102]),np.array([44.150741577148438,44.148887634277344,44.141479492187500,44.141479492187500,44.137775421142578,44.135925292968750,44.121109008789062,44.117404937744141,44.115554809570312,44.113704681396484,44.111850738525391,44.110000610351562]),linewidth=3,color='green')
    ax[sp].plot(np.array([28.688522338867188,28.690990447998047,28.690990447998047,28.693460464477539]),np.array([44.328517913818359,44.326667785644531,44.324813842773438,44.322963714599609]),linewidth=3,color='green')
    ax[sp].plot(np.array([28.644077301025391,28.646545410156250,28.649015426635742]),np.array([44.219257354736328,44.221111297607422,44.219257354736328]),linewidth=3,color='green')
    ax[sp].plot(np.array([28.663829803466797,28.668767929077148]),np.array([44.182220458984375,44.178516387939453]),linewidth=3,color='green')
    ax[sp].axis('equal')
    #ax[sp]._viewLim(Bbox([[np.min(fieldset.gridset.grids[sp].lon),np.max(fieldset.gridset.grids[sp].lon)],[np.min(fieldset.gridset.grids[sp].lat),np.max(fieldset.gridset.grids[sp].lat)]]))
  plt.show()
elif options['experiment']=='Galway':
  fig, (ax1, ax2) = plt.subplots(1, 2)
  x=data_xarray['lon'].values ; y=data_xarray['lat'].values ;
  #dx=fieldset.gridset.grids[0].lon[1]-fieldset.gridset.grids[0].lon[0] ; dy=fieldset.gridset.grids[0].lat[1]-fieldset.gridset.grids[0].lat[0]
  if options['beaching']>0:
    ax1.pcolormesh(fieldset.tmask.lon, fieldset.tmask.lat, fieldset.tmask.data[0,:,:], shading='auto')
  ax1.plot(x.T,y.T)   # transpose
  ax1.axis('equal')
  #ax[0]._viewLim(Bbox([[np.min(fieldset.gridset.grids[0].lon),np.max(fieldset.gridset.grids[0].lon)],[np.min(fieldset.gridset.grids[0].lat),np.max(fieldset.gridset.grids[0].lat)]]))
  z=data_xarray['z'].values;
  ax2.plot(x.T,z.T)
  plt.show()
'''

