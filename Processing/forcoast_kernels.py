# Custom kernels for ForCoast Service Module
import math
import parcels.rng as ParcelsRandom

from parcels import FieldSet, Field, NestedField, VectorField, ParticleFile, ParticleSet, JITParticle, ScipyParticle, Variable
from parcels import ErrorCode
from os import path
import numpy as np
import math

# Kernel functions:
# * AdvectionRK4                  2D advection (RK4 scheme) using just ocean currents
# * AdvectionRK3_3D               3D advection (RK4 scheme) using just ocean currents
# * Dbl_AdvectionRK4 :            2D advection (RK4 scheme)   using 2D currents and surface Stokes drift
# * Dbl_AdvectionEE_3D :          3D advection (Euler scheme) using 3D currents and horizontal Stokes drift at different depths
# * Dbl_AdvectionRK4_3D_clumsy :  3D advection (RK4 scheme)   using 3D currents and horizontal Stokes drift at different depths
# * stokes_profile :              function from OpenDrift to propagate (surface) Stokes drift to depth. The actual code is c/p to other kernels (no subroutine calls)
# * DiffusionUniformKh :          horizontal diffusion kernel with uniform diffusion coefficient Kh (i.e. only a random term, no cross term between advection and Wiener)
# * DiffusionZ :                  vertical advection (Milstein 1st order) using Kz from Nemo
# * DeleteParticle :              kernel that deletes a particle (e.g. after an out-of-box exception)
# * Frozenbeach :                 kernel that sets the particle beached attribute to true if it is beached (preventing further movement)
# * Unbeaching?D :                2D and 3D kernels that pushes a beached particle back into the sea
# * ckernel_dikes :               C-kernel to freeze or push back a particle horizontally crossing one segment out of a hard-coded set of segments

# Functions will to be converted to Kernel object before execution

def AdvectionRK4(particle, fieldset, time):
    #Advection of particles using fourth-order Runge-Kutta integration.
    if (particle.beached==0):
        (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt


def AdvectionRK4_3D(particle, fieldset, time):
    #Advection of particles using fourth-order Runge-Kutta integration including vertical velocity.
    #print("before rk4 depth=%g",particle.depth)
    if (particle.beached==0):
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt
    #print("after: [%g %g %g] (%g %g %g %g)",particle.lon,particle.lat,particle.depth,w1,w2,w3,w4)


def Dbl_AdvectionRK4(particle, fieldset, time):
    # Double advection of particles using fourth-order Runge-Kutta integration
    # using velocities from both current and Stokes drift
    if (particle.beached==0):
        (u1c, v1c) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
        (u1s, v1s) = fieldset.UVs[time, particle.depth, particle.lat, particle.lon]
        u1 = u1c + u1s
        v1 = v1c + v1s
        lon1, lat1 = (particle.lon + u1*.5*particle.dt, particle.lat + v1*.5*particle.dt)
        (u2c, v2c) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
        (u2s, v2s) = fieldset.UVs[time + .5 * particle.dt, particle.depth, lat1, lon1]
        u2 = u2c + u2s
        v2 = v2c + v2s
        lon2, lat2 = (particle.lon + u2*.5*particle.dt, particle.lat + v2*.5*particle.dt)
        (u3c, v3c) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
        (u3s, v3s) = fieldset.UVs[time + .5 * particle.dt, particle.depth, lat2, lon2]
        u3 = u3c + u3s
        v3 = v3c + v3s
        lon3, lat3 = (particle.lon + u3*particle.dt, particle.lat + v3*particle.dt)
        (u4c, v4c) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
        (u4s, v4s) = fieldset.UVs[time + particle.dt, particle.depth, lat3, lon3]
        u4 = u4c + u4s
        v4 = v4c + v4s
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt

    
def Dbl_AdvectionEE_3D(particle, fieldset, time):
    # Double advection of particles using Forward Euler integration including vertical velocity.
    if (particle.beached==0):
        (u1c, v1c, w1c) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        (u1s, v1s) = fieldset.UVs[time, particle.depth, particle.lat, particle.lon]
        wmp = fieldset.WMP[time, particle.depth, particle.lat, particle.lon]
        swh = fieldset.SWH[time, particle.depth, particle.lat, particle.lon]
        stokes_surface_speed = math.sqrt(u1s*u1s + v1s*v1s)
        fm02 = fm02 = 1. / wmp
        total_transport = (2.*math.pi/16.)*fm02*swh*swh
        k = (stokes_surface_speed/(2*total_transport))
        stokes_speed = stokes_surface_speed*math.exp(2*k*particle.depth)
        zeromask = stokes_surface_speed == 0
        u1sz = stokes_speed*u1s/stokes_surface_speed
        v1sz = stokes_speed*v1s/stokes_surface_speed
        if zeromask:
            u1sz = 0.0
            v1sz = 0.0
        u1 = u1c + u1sz
        v1 = v1c + v1sz
        particle.lon += u1 * particle.dt
        particle.lat += v1 * particle.dt
        particle.depth += w1c * particle.dt



def Dbl_AdvectionRK4_3D_clumsy(particle, fieldset, time):
    # Double advection of particles using RK4 integration including vertical velocity.
    if (particle.beached==0):
        (u1c, v1c, w1c) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        (u1s, v1s) = fieldset.UVs[time, particle.depth, particle.lat, particle.lon]
        wmp = fieldset.WMP[time, particle.depth, particle.lat, particle.lon]
        swh = fieldset.SWH[time, particle.depth, particle.lat, particle.lon]
        stokes_surface_speed = math.sqrt(u1s*u1s + v1s*v1s)
        fm02 = 1. / wmp
        total_transport = (2.*math.pi/16.)*fm02*swh*swh
        k = (stokes_surface_speed/(2*total_transport))
        stokes_speed = stokes_surface_speed*math.exp(2*k*particle.depth)
        u1sz = stokes_speed*u1s/stokes_surface_speed    ;   v1sz = stokes_speed*v1s/stokes_surface_speed
        if stokes_surface_speed == 0:
            u1sz = 0.0 ; v1sz = 0.0
        u1 = u1c + u1sz   ;   v1 = v1c + v1sz
        lon1,lat1,dep1 = (particle.lon + u1*.5*particle.dt , particle.lat + v1*.5*particle.dt , particle.depth + w1c*.5*particle.dt)

        t1=time + .5 * particle.dt
        (u2c, v2c, w2c) = fieldset.UVW[t1, dep1, lat1, lon1]
        (u2s, v2s) = fieldset.UVs[t1, dep1, lat1, lon1]
        wmp = fieldset.WMP[t1, dep1, lat1, lon1]
        swh = fieldset.SWH[t1, dep1, lat1, lon1]
        stokes_surface_speed = math.sqrt(u2s*u2s + v2s*v2s)
        fm02 = 1. / wmp
        total_transport = (2.*math.pi/16.)*fm02*swh*swh
        k = (stokes_surface_speed/(2*total_transport))
        stokes_speed = stokes_surface_speed*math.exp(2*k*dep1          )
        u2sz = stokes_speed*u2s/stokes_surface_speed    ;   v2sz = stokes_speed*v2s/stokes_surface_speed
        if stokes_surface_speed == 0:
            u2sz = 0.0 ; v2sz = 0.0
        u2 = u2c + u2sz ;     v2 = v2c + v2sz
        lon2,lat2,dep2 = (particle.lon + u2*.5*particle.dt , particle.lat + v2*.5*particle.dt , particle.depth + w2c*.5*particle.dt)

        t2=t1
        (u3c, v3c, w3c) = fieldset.UVW[t2, dep2, lat2, lon2]
        (u3s, v3s) = fieldset.UVs[t2, dep2, lat2, lon2]
        wmp = fieldset.WMP[t2, dep2, lat2, lon2]
        swh = fieldset.SWH[t2, dep2, lat2, lon2]
        stokes_surface_speed = math.sqrt(u3s*u3s + v3s*v3s)
        fm02 = 1. / wmp
        total_transport = (2.*math.pi/16.)*fm02*swh*swh
        k = (stokes_surface_speed/(2*total_transport))
        stokes_speed = stokes_surface_speed*math.exp(2*k*dep1          )
        u3sz = stokes_speed*u3s/stokes_surface_speed    ;   v3sz = stokes_speed*v3s/stokes_surface_speed
        if stokes_surface_speed == 0:
            u3sz = 0.0 ; v3sz = 0.0
        u3 = u3c + u3sz ;     v3 = v3c + v3sz
        lon3,lat3,dep3 = (particle.lon + u3*particle.dt , particle.lat + v3*particle.dt , particle.depth + w3c*particle.dt)

        t3=time + particle.dt
        (u4c, v4c, w4c) = fieldset.UVW[t3, dep3, lat3, lon3]
        (u4s, v4s) = fieldset.UVs[t3, dep3, lat3, lon3]
        wmp = fieldset.WMP[t3, dep3, lat3, lon3]
        swh = fieldset.SWH[t3, dep3, lat3, lon3]
        stokes_surface_speed = math.sqrt(u4s*u4s + v4s*v4s)
        fm02 = 1. / wmp
        total_transport = (2.*math.pi/16.)*fm02*swh*swh
        k = (stokes_surface_speed/(2*total_transport))
        stokes_speed = stokes_surface_speed*math.exp(2*k*dep1          )
        u4sz = stokes_speed*u4s/stokes_surface_speed    ;   v4sz = stokes_speed*v4s/stokes_surface_speed
        if stokes_surface_speed == 0:
            u4sz = 0.0 ; v4sz = 0.0
        u4 = u4c + u4sz ;     v4 = v4c + v4sz

        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.depth += (w1c + 2*w2c + 2*w3c + w4c) / 6. * particle.dt



def stokes_profile(stokes_u_surface, stokes_v_surface,significant_wave_height,mean_wave_period,z):
    """
    Code copied from OpenDrift

    Calculate vertical Stokes drift profile from Breivik et al. 2016,
    A Stokes drift approximation based on the Phillips spectrum, Ocean Mod. 100

    Parcels C compiler cannot use subroutines, so this code is copy/pasted 4 times in the 3D RK4 kernels above
    """

    stokes_surface_speed = np.sqrt(stokes_u_surface**2 + stokes_v_surface**2)
    fm02 = fm02 = 1. / mean_wave_period
    total_transport = (2.*np.pi/16.)*fm02*np.power(significant_wave_height, 2)
    k = (stokes_surface_speed/(2*total_transport))
    stokes_speed = stokes_surface_speed*np.exp(2*k*z)
    zeromask = stokes_surface_speed == 0
    stokes_u = stokes_speed*stokes_u_surface/stokes_surface_speed
    stokes_v = stokes_speed*stokes_v_surface/stokes_surface_speed
    stokes_u[zeromask] = 0
    stokes_v[zeromask] = 0
    return stokes_u, stokes_v




def DiffusionUniformKh(particle, fieldset, time):
    """Kernel for simple 2D diffusion where diffusivity (Kh) is assumed uniform.
    Assumes that fieldset has fields Kh_zonal and Kh_meridional.

    This kernel neglects gradients in the diffusivity field and is
    therefore more efficient in cases with uniform diffusivity.
    Since the perturbation due to diffusion is in this case spatially
    independent, this kernel contains no advection and can be used in
    combination with a seperate advection kernel.

    The Wiener increment `dW` should be normally distributed with zero
    mean and a standard deviation of sqrt(dt). Instead, here a uniform
    distribution with the same mean and std is used for efficiency and
    to keep random increments bounded. This substitution is valid for
    the convergence of particle distributions. If convergence of
    individual particle paths is required, use normally distributed
    random increments instead. See Gräwe et al (2012)
    doi.org/10.1007/s10236-012-0523-y for more information.

    at (lon=30,lat=42.5), 1m is 1.2153752066978996e-5° longitude
                          1m is 8.993512963684225e-6° latitude
    so Kx = 40 m²/s = (6.32m)²/s = (6.32*1.2153752066978996e-5)² [°²/s] = 5.908547572223849e-9 °²/s
       Ky = 40 m²/s =                                                   = 3.235331017118248e-9 °²/s

    """
    if (particle.beached==0):
        # Wiener increment with zero mean and std of sqrt(dt)
        dWx = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        dWy = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        
        # bx = sqrt( 2 * Kh_x) , idem for by
        bx = math.sqrt(2 * fieldset.Kx)  # [°²/s]  instead of : fieldset.Kh_zonal[time, particle.depth, particle.lat, particle.lon])  , roughly corresponds to 0.2 m2/s
        by = math.sqrt(2 * fieldset.Ky)  # [°²/s]  instead of : fieldset.Kh_meridional[time, particle.depth, particle.lat, particle.lon])

        particle.lon += bx * dWx
        particle.lat += by * dWy


def DiffusionZ(particule, fieldset, time):
    """ Kernel for vertical diffusion using Milstein scheme at first order
    The approximation done here is that diffusion is separated from advection instead of applying the increment all at once
    
    """
    if (particle.beached==0):
        dWz = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
        
        dz = min(0.9*particle.depth, fieldset.dres)
        Kzp1 = fieldset.Kz[time, particle.depth + dz, particle.lat, particle.lon] - fieldset.Kz_EVD[time, particle.depth + dz, particle.lat, particle.lon]
        Kzm1 = fieldset.Kz[time, particle.depth - dz, particle.lat, particle.lon] - fieldset.Kz_EVD[time, particle.depth - dz, particle.lat, particle.lon]
        dKdz = (Kzp1 - Kzm1) / (2 * dz)
        bz = math.sqrt(2 * (fieldset.Kz[time, particle.depth, particle.lat, particle.lon]-fieldset.Kz_EVD[time, particle.depth, particle.lat, particle.lon]))

        # prevent particles to cross sea-surface boundary
        z = particle.depth + 0.5 * dKdz * (dWz**2 + particle.dt) + bz * dWz
        particle.depth = max(z,0.05)

        # prevent particles to cross more than 1 layer / timestep, to cross bottom ???

def decay(particle, fieldset, time):
    particle.age += math.fabs(particle.dt)
    if particle.age>=fieldset.maxage:
        particle.delete()

    
def DeleteParticle(particle, fieldset, time):
    print("Particle [%d] lost !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()
def DeleteParticle1(particle, fieldset, time):
    print("Particle [%d] out-of-bounds !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.delete()
def DeleteParticle2(particle, fieldset, time):
    print("Particle [%d] through surface !! (%g %g %g %g)" % (particle.id, particle.lon, particle.lat, particle.depth, particle.time))
    particle.depth=-0.1


def Frozenbeach(particle,fieldset,time):
    # northsea: if u<1e-14 and v<1e-14 --> update part.lon and lat : advect using pre-computed unbeaching velocity
    # forcoast: if tmask==0 --> freeze particle (set particle.beached = 1)
    if particle.beached == 0:
        if fieldset.tmask[time, particle.depth, particle.lat, particle.lon]<0.5:
            particle.beached=1
            #print("particle %d beached!" % particle.id)
            
def Unbeaching2D(particle,fieldset,time):
    # northsea: if u<1e-14 and v<1e-14 --> update part.lon and lat : advect using pre-computed unbeaching velocity
    # forcoast: if tmask==0 --> advect using computed unbeaching velocity (by continuing or by going back)
    # galway/roms: tmask is a 2D mask ; a particle could still be advected BELOW the sea (but still have mask == 1)
    if fieldset.tmask[time, particle.depth, particle.lat, particle.lon]<0.5:
        in_sea=0
        """
        lon_increment=particle.lon-particle.prevlon
        lat_increment=particle.lat-particle.prevlat
        if abs(lon_increment)>=abs(lat_increment):   # try to push even further in longitude direction
            test_tmask=fieldset.tmask[time, particle.depth, particle.lat, particle.lon+lon_increment]
            if test_tmask>0.5:
                particle.lon += lon_increment
                in_sea=1
                #print("particle %d pushed further!" % particle.id)
        if (in_sea==0 and abs(lon_increment)<abs(lat_increment)):  # try in latitude direction
            test_tmask=fieldset.tmask[time, particle.depth, particle.lat+lat_increment, particle.lon]
            if test_tmask>0.5:
                particle.lat += lat_increment
                in_sea=1
                #print("particle %d pushed further!" % particle.id)
        """
        if (in_sea==0):                         # if that did not work, push back to where it came from
            particle.lon = particle.prevlon
            particle.lat = particle.prevlat
            #print("particle %d pushed back!" % particle.id)

def Unbeaching3D_nemo(particle,fieldset,time):
    # northsea: if u<1e-14 and v<1e-14 --> update part.lon and lat : advect using pre-computed unbeaching velocity
    # forcoast: if tmask==0 --> check if pushing the particle back up brings it back in the sea, if not use 2D unbeaching
    if fieldset.tmask[time, particle.depth, particle.lat, particle.lon]<0.5:
        if particle.prevdep<particle.depth: # depth is zero at surface and increases going downward (nemo convention)
            if (fieldset.tmask[time, particle.prevdep, particle.lat, particle.lon]>0.5):
                particle.depth=particle.prevdep
                #print("particle %d pushed up!" % particle.id)
    # if that did not work: abandon for now, particle will be unbeached in the horizontal

def Unbeaching3D_roms(particle,fieldset,time):
    '''
    print("entering unbeaching3D_roms")
    
    print("particle depth %g %g",particle.depth,fieldset.depth_w[time,-1,particle.lon,particle.lat])
    if particle.depth>fieldset.depth_w[time,-1,particle.lon,particle.lat]:
        particle.depth=fieldset.depth_w[time,-1,particle.lon,particle.lat]-0.05

    if particle.depth<fieldset.depth_w[time,0,particle.lon,particle.lat]:
        particle.depth=fieldset.depth_w[time,0,particle.lon,particle.lat]+0.05

    if fieldset.U[time,particle.depth,particle.lat,particle.lon]>20.0:
        if particle.prevdep<particle.depth:
            if fieldset.U[time, particle.prevdep, particle.lat, particle.lon]<20.0:
                particle.depth=particle.prevdep
    '''

def save_position(particle,fieldset,time):
    # this should be the last kernel called !
    particle.prevlon=particle.lon
    particle.prevlat=particle.lat
    particle.prevdep=particle.depth



######################################### C KERNELS ###################################################################
# define C kernel function to not cross dikes
def ckernel_dikes(particle, fieldset, time):
    func('parcels_customed_Cfunc_pointer_args',    particle.lon,particle.lat,particle.depth,   particle.beached,  particle.prevlon,particle.prevlat,particle.prevdep)
             # fieldset.dikes['points'],fieldset.dikes['lon'],fieldset.dikes['lat'],fieldset.dikes['restarts'],fieldset.dikes['restart']
