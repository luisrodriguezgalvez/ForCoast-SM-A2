#!/data/anaconda3/envs/py3_parcels/bin/python

# ForCoast Pilot 7 -- Ocean Parcels Lagrangian tracking

from parcels import FieldSet, Field, NestedField, VectorField
from parcels import ParticleSet, ParticleFile, Variable, ErrorCode
from parcels import JITParticle #, ScipyParticle
from parcels import plotTrajectoriesFile
from parcels.kernel import Kernel
from glob import glob
import numpy as np
from datetime import timedelta as delta
import time
from forcoast_loaders import *
from forcoast_kernels import *
import yaml
import os
from pathlib import Path
import sys

argv = sys.argv[1:]

# OPTIONS
with open('pilot7.yaml') as f:
  options=yaml.load(f,Loader=yaml.Loader)
npart=options['npart']

# Replace options from yml file with options passed from command line (example; code could be improved) 
# Example: 
# - Windows: python forcoast.py 2021-02-01 4 c:\data
# - Linux: python forcoast.py 2021-02-01 4 /data
if not argv:
    print("Use input from yml file")
else:
    print("Replace selected input with input from command line")
    options['sdate'] = argv[0]
    options['simlength'] = int(argv[1])
    options['PHY_path'] = argv[2]
    options['out_filename'] = options['PHY_path'] + '/' + "EforieParticles.nc"

# LOAD DATA
if options['run3D']:
    print('running in 3D')
else:
    print('running in 2D')

# EFORIE
ufiles = sorted(glob(str(options['PHY_path'] + '/' + '2_EFORIE_1h_*_grid_U_*.nc')))
vfiles = sorted(glob(str(options['PHY_path'] + '/' + '2_EFORIE_1h_*_grid_V_*.nc')))
wfiles = sorted(glob(str(options['PHY_path'] + '/' + '2_EFORIE_1h_*_grid_W_*.nc')))

mesh_mask = str(Path(options['PHY_path']  + '/' + '2_mesh_mask.nc'))

indices = {'lon': range(1,107), 'lat': range(2,211)} # NEMO puts zero along the (ghost) boundaries
EFORIE=get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,run3D=options['run3D'],indices=indices,chunksize=False,vdiffusion=options['vdiffusion'],beaching=options['beaching'])

if options['nesting']==True:
    print('using 3 nested grids')

    # NWS
    ufiles = sorted(glob(options['PHY_path']+'*/*/1_NWS_1h_*_grid_U_*.nc*'))
    vfiles = sorted(glob(options['PHY_path']+'*/*/1_NWS_1h_*_grid_V_*.nc*'))
    wfiles = sorted(glob(options['PHY_path']+'*/*/1_NWS_1h_*_grid_W_*.nc*'))
    mesh_mask = options['PHY_path'] + '1_mesh_mask.nc'
    indices = {'lon': range(1,425), 'lat': range(3,346)} # NEMO puts zero along the (ghost) boundaries
    NWS=get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,run3D=options['run3D'],indices=indices,chunksize=False,vdiffusion=options['vdiffusion'],beaching=options['beaching'])

    # BSFS
    ufiles = sorted(glob(options['PHY_path']+'*/*/BS_1h_*_grid_U_*.nc*'))
    vfiles = sorted(glob(options['PHY_path']+'*/*/BS_1h_*_grid_V_*.nc*'))
    wfiles = sorted(glob(options['PHY_path']+'*/*/BS_1h_*_grid_W_*.nc*'))
    mesh_mask = options['PHY_path'] + 'mesh_mask_59levels.nc'
    BSFS=get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,run3D=options['run3D'],chunksize=False,vdiffusion=options['vdiffusion'],beaching=options['beaching'])

    # NestedFields
    U=NestedField('U', [EFORIE.U, NWS.U, BSFS.U])
    V=NestedField('V', [EFORIE.V, NWS.V, BSFS.V])
    if options['run3D']:
        W=NestedField('W', [EFORIE.W, NWS.W, BSFS.W])
        if options['vdiffusion']:
            Kz=NestedField('Kz', [EFORIE.Kz, NWS.Kz, BSFS.Kz])
            Kz_EVD=NestedField('Kz_EVD', [EFORIE.Kz_EVD, NWS.Kz_EVD, BSFS.Kz_EVD])
    if options['beaching']>0:
        tmask=NestedField('tmask', [EFORIE.tmask, NWS.tmask, BSFS.tmask])
    del EFORIE ; del NWS ; del BSFS
    
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
    fieldset=EFORIE
    del EFORIE

if options['stokes']==True:
    print('using currents and Stokes drift')
    mesh_mask  = options['WAV_path'] + 'BS-MFC_007_003_coordinates.nc'
    get_wav_fields(fieldset,options['WAV_path'],options['run3D'],mesh_mask,chunksize=False)
else:
    print('using just currents, no Stokes drift')
fieldset.add_constant('maxage', options['maxage'])

# CREATE PARTICLE SET
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
    c_includefile = path.join('c_kernels/crossdike1.h') # path.dirname(__file__),
elif options['beaching']==2:
    print('Un-beaching beached particles')
    if (options['run3D']):
        kernels += Unbeaching3D_nemo
    kernels += Unbeaching2D
    c_includefile = path.join('c_kernels/crossdike2.h') # path.dirname(__file__),
# prevent cross-dike
if options['beaching']>0:
    with open(c_includefile,'r') as f:
        c_include=f.read()
    kernels += pset.Kernel(ckernel_dikes, c_include=c_include, delete_cfiles=False)
    kernels += save_position
else:
    print('Beaching: no action')
    

# RUN SIMULATION    
output_file = pset.ParticleFile(name=options['out_filename'], outputdt=delta(minutes=options['savedt']))
t0=time.time()
pset.execute(kernels, runtime=delta(days=options['simlength']), dt=delta(minutes=options['dt']), output_file=output_file, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle, ErrorCode.ErrorThroughSurface: DeleteParticle})
              # verbose_progress=False)
t1=time.time() ; totaltime=t1-t0 ; print("computing time :",totaltime)
output_file.close()


# POST-TREATMENT DETECTION
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


## PLOT PARCELS RESULTS
#plt=plotTrajectoriesFile('EforieParticles.nc')

## PCOLOR TRAJECTORIES OVER TMASK
# DT 20-5-2021: Shouldn't this be in postprocessing? Gave an error, so commented out instead of fixing for now.
# import xarray as xr
# import matplotlib.pyplot as plt
# data_xarray = xr.open_dataset('EforieParticles.nc')  # trajectories, time, lon, lat, z, prevlon, prevlat, prevdep, beached [obs=timestep, traj=particle]
# x=data_xarray['lon'].values ; y=data_xarray['lat'].values ;
# fig,ax=plt.subplots(1,3)
# for sp in range(3):
#     dx=fieldset.gridset.grids[sp].lon[1]-fieldset.gridset.grids[sp].lon[0] ; dy=fieldset.gridset.grids[sp].lat[1]-fieldset.gridset.grids[sp].lat[0]
#     ax[sp].pcolormesh(fieldset.gridset.grids[sp].lon-(dx/2), fieldset.gridset.grids[sp].lat-(dy/2), fieldset.tmask[sp].data[0,0,:,:], shading='auto')
#     ax[sp].plot(x.T,y.T)
#     #plot dikes
#     ax[sp].plot(np.array([28.671237945556641,28.673706054687500,28.673706054687500,28.676176071166992,28.681114196777344,28.681114196777344,28.700866699218750,28.700866699218750,28.703336715698242,28.703336715698242,28.705804824829102,28.705804824829102]),np.array([44.150741577148438,44.148887634277344,44.141479492187500,44.141479492187500,44.137775421142578,44.135925292968750,44.121109008789062,44.117404937744141,44.115554809570312,44.113704681396484,44.111850738525391,44.110000610351562]),linewidth=3,color='green')
#     ax[sp].plot(np.array([28.688522338867188,28.690990447998047,28.690990447998047,28.693460464477539]),np.array([44.328517913818359,44.326667785644531,44.324813842773438,44.322963714599609]),linewidth=3,color='green')
#     ax[sp].plot(np.array([28.644077301025391,28.646545410156250,28.649015426635742]),np.array([44.219257354736328,44.221111297607422,44.219257354736328]),linewidth=3,color='green')
#     ax[sp].plot(np.array([28.663829803466797,28.668767929077148]),np.array([44.182220458984375,44.178516387939453]),linewidth=3,color='green')
#     ax[sp].axis('equal')
#     #ax[sp]._viewLim(Bbox([[np.min(fieldset.gridset.grids[sp].lon),np.max(fieldset.gridset.grids[sp].lon)],[np.min(fieldset.gridset.grids[sp].lat),np.max(fieldset.gridset.grids[sp].lat)]]))
# plt.show()


