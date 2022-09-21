#(C) Copyright FORCOAST H2020 project under Grant No. 870465. All rights reserved.

# ForCoast loader : loads netcdf files from NEMO, from CMEMS WAV, ...

from parcels import FieldSet, Field, NestedField, VectorField, ParticleFile, ParticleSet, JITParticle, Variable
from parcels import ErrorCode
from parcels.tools.converters import Geographic,GeographicPolar
import numpy as np
from datetime import timedelta as delta
from glob import glob
from pprint import pprint

# create a new fieldset from NEMO currents (1 grid)
# if 3D run, read also vertical currents
# if required, read also vertical diffusion coefficient
# if required, read also tmask (for unbeaching)
def get_nemo_fields(ufiles,vfiles,wfiles,mesh_mask,**kwargs):
    filenames =  {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                  'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles}}
    variables =  {'U': 'vozocrtx', 'V': 'vomecrty'}
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw'  , 'time': 'time_counter'},
                  'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw'  , 'time': 'time_counter'}}

    run3D=kwargs.get('run3D',False)
    if run3D:
        filenames['W']  = {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}
        variables['W']  = 'vovecrtz'
        dimensions['W'] = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

    vdiffusion=kwargs.get('vdiffusion',False)
    if vdiffusion:
        filenames['Kz']  = {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}
        variables['Kz']  = 'votkeavs'
        dimensions['Kz'] = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw', 'time': 'time_counter'}
        filenames['Kz_EVD']  = {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}
        variables['Kz_EVD']  = 'voevdavt'
        dimensions['Kz_EVD'] = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'}

    beaching=kwargs.get('beaching',0)
    if beaching>0:
        filenames['tmask'] = {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': mesh_mask}
        variables['tmask'] = 'tmask'
        dimensions['tmask'] = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw'}
        
    indices=kwargs.get('indices',None)
    fieldset=FieldSet.from_nemo(filenames, variables, dimensions, indices=indices)
        # this includes: fieldset.W.set_scaling_factor(-1.)
    
    if run3D:
        def compute(fieldset):
            fieldset.W.data[:, 0, :, :] = 0.
        fieldset.compute_on_defer = compute

    # at this stage, fieldset contains: U,V,compute_on_defer,gridset,time_origin
    return fieldset


# get_mohid_fields : create a new fieldset from MOHID currents (all variables co-located)
def get_mohid_fields(files,**kwargs):
    filenames = {'U': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                 'V': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files}}
    variables =  { 'U': 'u', 'V': 'v'}
    dimensions = {'U':       {'lon': 'lon', 'lat': 'lat', 'depth': 'depth' , 'time': 'time'},
                  'V':       {'lon': 'lon', 'lat': 'lat', 'depth': 'depth' , 'time': 'time'}}
    run3D=kwargs.get('run3D',True)
    if run3D:
        filenames['W'] = {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files}
        variables['W'] = 'velocity_W'
        dimensions['W']= {'lon': 'lon', 'lat': 'lat', 'depth': 'depth' , 'time': 'time'}
    beaching=kwargs.get('beaching',0)
    if beaching>0:
        filenames['tmask'] = {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files[0]}
        variables['tmask'] = 'mask'
        dimensions['tmask']= {'lon': 'lon', 'lat': 'lat', 'depth': 'depth'}

    indices=kwargs.get('indices',None)
    fieldset=FieldSet.from_netcdf(filenames, variables, dimensions, indices=indices, vmax=1.0e36) #allow_time_extrapolation=True

    fieldset.U.vmin=-1.0e14
    fieldset.V.vmin=-1.0e14
    
    if run3D:
        fieldset.W.vmin=-1e+14
        def compute(fieldset):
            fieldset.W.data[:, 0, :, :] = 0.
        fieldset.compute_on_defer = compute

        
    return fieldset

# get_MIT_fields: create a new fieldset from MIGgcm currents (all variables co-located)
def get_MIT_fields(files,maskfile,**kwargs):
    filenames = {'U': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                 'V': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files}}
    variables =  { 'U': 'uo', 'V': 'vo'}
    dimensions = {'U':       {'lon': 'longitude', 'lat': 'latitude', 'depth': 'depth' , 'time': 'time'},
                  'V':       {'lon': 'longitude', 'lat': 'latitude', 'depth': 'depth' , 'time': 'time'}}
    run3D=kwargs.get('run3D',True)
    if run3D:
        filenames['W'] = {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files}
        variables['W'] = 'wo'
        dimensions['W']= {'lon': 'longitude', 'latitude': 'lat', 'depth': 'depth' , 'time': 'time'}
    beaching=kwargs.get('beaching',0)
    if beaching>0:
        filenames['tmask'] = {'lon': maskfile, 'lat': maskfile, 'depth': maskfile, 'data': maskfile}
        variables['tmask'] = 'mask'
        dimensions['tmask']= {'lon': 'longitude', 'lat': 'latitude', 'depth': 'depth'}

    indices=kwargs.get('indices',None)
    fieldset=FieldSet.from_netcdf(filenames, variables, dimensions, indices=indices, vmax=1.0e36) #allow_time_extrapolation=True
    
    fieldset.U.vmin=-1.0e14
    fieldset.V.vmin=-1.0e14

    if run3D:
        fieldset.W.vmin=-1e+14
        def compute(fieldset):
            fieldset.W.data[:, 0, :, :] = 0.
        fieldset.compute_on_defer = compute

    return fieldset


# get_roms_fields : create a new fieldset from ROMS currents
#     grids are time-varying
#     no hypothesis on grids such as C-grid etc: each variable is defined on its own grid
#     masks are 2D fields
def get_roms_fields(files,**kwargs):
    run3D=kwargs.get('run3D',False)
    if run3D:
        filenames = {'depth_u':{'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                     'depth_v':{'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                     'depth_w':{'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                     'U': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                     'V': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                     'W': {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files},
                      'zeta':   {'lon': files[0], 'lat': files[0],                    'data': files}}
        variables =  { 'depth_u': 'z_u', 'depth_v': 'z_v', 'depth_w': 'z_w',
                       'U': 'u', 'V': 'v', 'W': 'w',
                       'zeta': 'zeta'}
        dimensions = {'depth_u': {'lon': 'lon_u',   'lat': 'lat_u',   'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'depth_v': {'lon': 'lon_v',   'lat': 'lat_v',   'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'depth_w': {'lon': 'lon_rho', 'lat': 'lat_rho', 'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'U':       {'lon': 'lon_u',   'lat': 'lat_u',   'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'V':       {'lon': 'lon_v',   'lat': 'lat_v',   'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'W':       {'lon': 'lon_rho', 'lat': 'lat_rho', 'depth': 'not_yet_set'  , 'time': 'ocean_time'},
                      'zeta':    {'lon': 'lon_rho', 'lat': 'lat_rho',                           'time': 'ocean_time'}}
    else:
        filenames = {'U': {'lon': files[0], 'lat': files[0], 'data': files},
                     'V': {'lon': files[0], 'lat': files[0], 'data': files}}
        variables =  { 'U': 'ubar',
                       'V': 'vbar'}
        dimensions = {'U':       {'lon': 'lon_u', 'lat': 'lat_u', 'time': 'ocean_time'},
                      'V':       {'lon': 'lon_v', 'lat': 'lat_v', 'time': 'ocean_time'}}

    vdiffusion=kwargs.get('vdiffusion',False)
    if vdiffusion:
        filenames['Kz']  = {'lon': files[0], 'lat': files[0], 'depth': files[0], 'data': files}
        variables['Kz']  = 'AKt'
        dimensions['Kz'] = {'lon': 'lon_rho', 'lat': 'lat_rho', 'depth': 'not_yet_set', 'time': 'ocean_time'}

    beaching=kwargs.get('beaching',0)
    if beaching>0:
        filenames['tmask']     = {'lon': files[0], 'lat': files[0], 'data': files[0]} # using s-layers, the mask is a 2D field (no depth dimension)
        variables['tmask']     = 'mask_rho'
        dimensions['tmask']    = {'lon': 'lon_rho', 'lat': 'lat_rho'}
        
    indices=kwargs.get('indices',None)
    fieldset=FieldSet.from_netcdf(filenames, variables, dimensions, indices=indices, vmax=1.0e36) #allow_time_extrapolation=True
    if run3D:
        fieldset.U.set_depth_from_field(fieldset.depth_u)
        fieldset.V.set_depth_from_field(fieldset.depth_v)
    fieldset.U.vmax=1.0e36
    fieldset.V.vmax=1.0e36
    if run3D:
        fieldset.W.set_depth_from_field(fieldset.depth_w)
        fieldset.W.vmax=1.0e36

    #fieldset.depth_u.set_scaling_factor(-1)
    #fieldset.depth_v.set_scaling_factor(-1)
    #if run3D:
    #    fieldset.depth_w.set_scaling_factor(-1)
    #    fieldset.W.set_scaling_factor(-1)

    return fieldset



# add to an existing fieldset, a vectorfield containing CMEMS WAV stokes drift currents
def get_wav_fields(fieldset,path,run3D,mesh_mask,**kwargs):
    basepath = path + '*/*/' + '/*WAVES*.nc'
    fnames=[];
    fnames += sorted(glob(str(basepath)))

    filenames  = {'Us': {'lon': mesh_mask, 'lat': mesh_mask, 'data': fnames},   'Vs': {'lon': mesh_mask, 'lat': mesh_mask, 'data': fnames}}
    dimensions = {'Us': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'},'Vs': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}}
    variables  = {'Us': 'VSDX', 'Vs': 'VSDY'}
    fs=FieldSet.from_netcdf(filenames,variables,dimensions)
    fieldset.add_field(fs.Us) ; fieldset.Us.units = GeographicPolar()
    fieldset.add_field(fs.Vs) ; fieldset.Vs.units = Geographic()
    UVs = VectorField('UVs', fs.Us, fs.Vs)
    fieldset.add_vector_field(UVs)

    if run3D:
        filenames  = {'WMP': {'lon': mesh_mask, 'lat': mesh_mask, 'data': fnames},   'SWH': {'lon': mesh_mask, 'lat': mesh_mask, 'data': fnames}}
        dimensions = {'WMP': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'},'SWH': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}}
        variables  = {'WMP': 'VTM10', 'SWH': 'VHM0'}
        W3Dfs = FieldSet.from_netcdf(filenames,variables,dimensions)
        fieldset.add_field(W3Dfs.WMP)
        fieldset.add_field(W3Dfs.SWH)
    
    return fieldset
