#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import importlib
import shutil
import platform

OSX = platform.system()

from scipy.interpolate import make_interp_spline

import sys
if OSX == 'Darwin':
    sys.path.append('/Users/psummers8/Documents/MITgcm/MITgcm/elizaScripts/main_scripts')
else:
    sys.path.append('/storage/home/hcoda1/2/psummers8/MITgcmSandbox/elizaScripts/main_scripts')
import build_domain_funcs as build_domain 
import run_config_funcs as rcf # import helpter functions

setupNotes = open("setupReport.txt", "w") 

# ## Main run configuration
email = 'psummers8@gatech.edu'
# set high level run configurations

run_config = {}
run_config['ncpus_xy'] = [7, 1] # cpu distribution in the x and y directions
run_config['run_name'] = 'bergPlume_200'
run_config['ndays'] = 20 # simulaton time (days)
run_config['test'] = False # if True, run_config['nyrs'] will be shortened to a few time steps

run_config['horiz_res_m'] = 200 # horizontal grid spacing (m)
run_config['Lx_m'] = 35000 # domain size in x (m)
run_config['Ly_m'] = 5400 # domain size in y (m)
# NOTE: the number of grid points in x and y should be multiples of the number of cpus.

#run_config['evolve_salt'] = False
run_config['use_GMRedi'] = False # should be set to false for eddy permitting resolutions
run_config['periodic_forcing'] = False # note: the code is not yet set up to handle time-dependent forcing

MITgcm_release = 'MITgcm-checkpoint68z' #Sept 2024 release
#MITgcm_code_dir = os.path.join(group_home_dir, 'shared/mitgcm_releases', MITgcm_release)

# you probably don't need to touch this
run_config['use_MPI'] = True # for multi-processing
run_config['lf'] = '\r\n' # linebreak characters 
if OSX == 'Darwin':
    run_config['exps_dir'] = os.path.join('/Users/psummers8/Documents/MITgcm/MITgcm/experiments') 
else:
    run_config['exps_dir'] = os.path.join('/storage/home/hcoda1/2/psummers8/MITgcmSandbox/experiments') 
run_config['run_dir'] = os.path.join(run_config['exps_dir'], run_config['run_name'])
print('run_config is', run_config)
setupNotes.write('run_config is ' + str(run_config)+'\n')

# ## Generate new experiment directory and copy over defaults

# In[3]:


# create experimentary directory on SCRATCH and copy over default configuration
# NOTE: this step does not overwrite existing directories. 
run_subdir_list = ['build', 'code', 'input', 'results']
for subdir in run_subdir_list:
    run_config['%s_dir'% subdir] = os.path.join(run_config['run_dir'], subdir)
    os.makedirs(run_config['%s_dir'% subdir], exist_ok=True)
     
# copy over defaults
if OSX == 'Darwin':
    default_dirs = os.listdir('/Users/psummers8/Documents/MITgcm/MITgcm/DEFAULT_Berg/')
else:
    default_dirs = os.listdir('/storage/home/hcoda1/2/psummers8/MITgcmSandbox/DEFAULT_Berg/')
for dir00 in default_dirs:
    if dir00.startswith('.'):
        continue
        
    if OSX == 'Darwin':
        default_dir = '/Users/psummers8/Documents/MITgcm/MITgcm/DEFAULT_Berg/%s/'%dir00
    else:
        default_dir = '/storage/home/hcoda1/2/psummers8/MITgcmSandbox/DEFAULT_Berg/%s/'%dir00    
    default_files = os.listdir(default_dir)
    dst_dir = os.path.join(run_config['run_dir'], dir00)
    
    for file in default_files:

        if file.startswith('.'):
            continue
        else:
            src_fpath = os.path.join(default_dir, file)
            shutil.copy2(src_fpath, dst_dir)
            #print(src_fpath, '>', dst_dir)


# In[4]:


print(run_config['run_dir'])
print(os.listdir(run_config['run_dir']))
setupNotes.write('run directory and subdirectories:\n')
setupNotes.write(str(run_config['run_dir'])+'\n')
setupNotes.write(str(os.listdir(run_config['run_dir']))+'\n')

# In[5]:


# create new analysis sub-dir in your home directory
if OSX == 'Darwin':
    analysis_dir = '/Users/psummers8/Documents/MITgcm/MITgcm/analysis/%s'%run_config['run_name']
else:
    analysis_dir = '/storage/home/hcoda1/2/psummers8/MITgcmSandbox/analysis/%s'%run_config['run_name']
os.makedirs(analysis_dir, exist_ok=True)



# ## Domain  and grid parameters

# In[6]:


secsInDay = 24*60*60
secsInYear = 365*secsInDay


# set domain size
domain_params = {}
domain_params['Lx'] = run_config['Lx_m'] # domain size in x (m)
domain_params['Ly'] = run_config['Ly_m'] # domain size in y (m)
domain_params['L_sponge'] = 5000 # width of eastern sponge layer (m)
domain_params['H'] = 200 # max domain depth (m)


# In[7]:


#---grid parameters ---# 

# NOTE: the only thing you may need to change here is the number of z-grid points)
grid_params = {}

grid_params['nSx'] = 1 # num of tiles per processor in x-direction
grid_params['nSy'] = 1 # num of tiles per processor in y-direction
grid_params['nTx'] = 1 # num of threads per processor in x-direction
grid_params['nTy'] = 1 # num of threads per processor in y-direction
grid_params['OLx'] = 3 # num of overlapping x-gridpoints per tile
grid_params['OLy'] = 3 # num of overlapping y-gridpoints per tile
grid_params['Nr'] = 70 # num of z-grid points

grid_params['nPx'] = run_config['ncpus_xy'][0] #num of processors in x-direction
grid_params['nPy'] = run_config['ncpus_xy'][1] #num of processors in x-direction

# grid_params['nSx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points in sub grid
# grid_params['nSy'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points in sub grid

# grid_params['Nx'] = grid_params['sNx'] * grid_params['nSx'] * grid_params['nPx']
# grid_params['Ny'] = grid_params['sNy'] * grid_params['nSy'] * grid_params['nPy']

grid_params['Nx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points
grid_params['Ny'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points

print("Nx: %s" %grid_params['Nx'])
print("Ny: %s" %grid_params['Ny'])
setupNotes.write("Nx: %s" %grid_params['Nx']+'\n')
setupNotes.write("Ny: %s" %grid_params['Ny']+'\n')

grid_params['sNx'] = grid_params['Nx']/grid_params['nPx']#num of x-gridpoints per tile
grid_params['sNy'] = grid_params['Ny']/grid_params['nPy'] #num of y-gridpoints per tile

print("sNx: %s" %grid_params['sNx'])
print("sNy: %s" %grid_params['sNy'])
setupNotes.write("sNx: %s" %grid_params['sNx']+'\n')
setupNotes.write("sNy: %s" %grid_params['sNy']+'\n')

# NOTE: sNx and sNy should be whole numbers/integers. As long we keep the horizontal resolution,
# domain dimesions, and number of cpus to be multiples of five, we should be ok. 

for key, param  in grid_params.items():
    assert param%1==0, "grid parameter needs to be an integer"
    grid_params[key] = int(param)
    
print(grid_params)
setupNotes.write('Grid parameters\n')
setupNotes.write(str(grid_params)+'\n')
#run_config['grid_params'] = grid_params


# In[8]:


# grid_params cont'd
grid_params['usingCartesianGrid'] = True
grid_params['usingSphericalPolarGrid'] = False 

# horizontal grid spacing
grid_params['delX'] = (domain_params['Lx']/grid_params['Nx'])*np.ones(grid_params['Nx'])
grid_params['delY'] = (domain_params['Ly']/grid_params['Ny'])*np.ones(grid_params['Ny'])


# vertical grid spacing 
# spacing increases with depth---can be modified
# zidx = np.arange(1, grid_params['Nr']+1)
# aa = 10
# dz1 = 2*domain_params['H']/grid_params['Nr']/(aa+1)
# dz2 = aa*dz1
# dz = dz1 + ((dz2-dz1)/2)*(1+np.tanh((zidx-((grid_params['Nr']+1)/2))/aa))
# zz1 = np.append([0], np.cumsum(dz))
# zz = -(zz1[:-1] + np.diff(zz1)/2) # layer midpoints

dz = domain_params['H']/grid_params['Nr']*np.ones(grid_params['Nr']);

grid_params['delZ'] = dz
grid_params['hFacMinDr'] = dz.min()


# In[ ]:





# 
# ## Physical parameters

# In[9]:


#---physical params---#

params01 = {} 

# physical constants
g = 9.81 # acc. due to gravity (m/s**2)
Omega = 2*np.pi*366/365/86400 # planetary rotation rate 
Rp = 6400*1000 # planetary radius (m)
lat_min = -70 # latitude at southern boundary (degrees)
#f0 = 2*Omega*np.sin(np.deg2rad(lat_min)) # coriolis param (1/s)
#beta = (2*Omega*np.cos(np.deg2rad(lat_min))/Rp) # beta param


# momentum scheme
params01['vectorInvariantMomentum'] = True

#Note: here and elsewhere, we need to be explicit about floats vs ints. E.g., use 12.0 to represent float and
# 12 for int

# viscosity parameters
params01['viscA4'] = 0.0000 # Biharmonic viscosity?
params01['viscAz'] = 1.0e-3 # Vertical viscosity
params01['viscAh'] = 2.5e-1 # Vertical viscosity
#params01['viscC2smag'] = 2.2 # ??? viscosity

# advection and time stepping
params01['tempAdvScheme'] = 33 # needs to be int
params01['saltAdvScheme'] = 33 # needs to be int
#params01['tempStepping'] = True
#params01['saltStepping'] = run_config['evolve_salt']
params01['staggerTimeStep'] = True

# diffusivity
#params01['diffK4T'] = 0.0e4 # ?? temp diffusion
params01['diffKhT'] = 2.5e-1 # Horizontal temp diffusion
params01['diffKzT'] = 2.0e-5 # Vertical temp diffusion
#params01['diffK4S'] = 0.0e4 # ?? salt diffusion
params01['diffKhS'] = 2.5e-1 # Horz salt diffusion
#params01['diffKzS'] = 2.0e-5 # Horz salt diffusion


# equation of state
params01['eosType'] = 'JMD95Z'
# params01['eosType'] = 'LINEAR'
# params01['tAlpha'] = 0.4e04
# params01['sBeta'] = 8.0e-4
params01['Tref'] = np.ones(grid_params['Nr'])*0. #ref temp
params01['Sref'] = np.ones(grid_params['Nr'])*34. #ref salt

# boundary conditions
#params01['bottomDragLinear'] = 0.0e-4

params01['no_slip_sides'] = True
params01['no_slip_bottom'] = True
params01['rigidLid'] = False
params01['implicitFreeSurface'] = True
params01['selectAddFluid'] = 1
#params01['implicitViscosity'] = True
#params01['implicitDiffusion'] = True

# physical parameters
#params01['f0'] = -1.36e-4
params01['f0'] = 1.36e-4
params01['beta'] = 0.0e-13
params01['gravity'] = g

# misc
params01['hFacMin'] = 0.05
params01['nonHydrostatic'] = True
params01['readBinaryPrec'] = 64
#params01['plotLevel'] = 0

# AS: Don't allow partial cell height to fall below min grid spacing
#params01['hFacMinDr'] = np.min(dz)


# ## Check for numericl stability?

# ## Numeric solvers and I/O controls

# In[10]:


# numeric solver parameters 

params02 = {}
params02['cg2dMaxIters'] = 300
params02['cg2dTargetResidual'] = 1e-13
params02['cg3dMaxIters'] = 20
params02['cg3dTargetResidual'] = 1e-8

# time stepping parameters 
params03 = {}
params03['nIter0'] = 0
#params03['endTime'] = 864000.0
deltaT = 20.0
params03['abEps'] = 0.1

#if run_config['testing']:
    
params03['chkptFreq'] = 0.0
params03['pChkptFreq'] = 0.0
params03['taveFreq'] = 0.0
params03['dumpFreq'] = 864000.0
params03['taveFreq'] = 0.0
params03['monitorFreq'] = 86400.0
params03['monitorSelect'] = 1


params03['periodicExternalForcing'] = False
params03['ExternForcingPeriod'] = 100.0
params03['ExternForcingCycle'] = 1000.0 


# In[11]:


if run_config['test']:
    nTimeSteps = 10
else:
    nTimeSteps = np.ceil(run_config['ndays']*secsInDay/deltaT)

simTimeAct = nTimeSteps*deltaT

params03['endTime'] = int(params03['nIter0']*deltaT+simTimeAct)
params03['deltaT'] = np.round(deltaT)
grid_params['Nt'] = nTimeSteps


# ## Create 'data' files

# In[12]:


# NOTE: These steps generate the data text files located int the input directory


# In[13]:


# gather params for data file 
params04 = {} #<-- using params04 to be consistent with ordering in Andrew's code
params04['usingCartesianGrid'] = grid_params['usingCartesianGrid']
params04['usingSphericalPolarGrid'] = grid_params['usingSphericalPolarGrid']
params04['delX']  = grid_params['delX']
params04['delY'] = grid_params['delY']
params04['delZ'] = dz


# In[14]:


# get data fnames param
params05 = {}
params05['bathyFile'] ='topog.slope'
params05['hydrogThetaFile'] = 'T.init'
params05['hydrogSaltFile'] = 'S.init'


# In[15]:


data_params = [params01, params02, params03, params04, params05]
rcf.write_data(run_config, data_params, group_name='data', lf=run_config['lf'])


# ## Make SIZE.h File

# In[16]:


# generate size.h file
rcf.createSIZEh(run_config, grid_params)


# ## Specifiy relaxation of temperature and salinity 

# In[17]:


OBCS = {}


# ## Set boundary conditions

# In[ ]:





# ## Specify Diagnostics

# In[18]:


# Here we specify variables that should saved (i.e., written to disk) at various time intervals

# adjust output frequency
if run_config['test']:
    run_config['inst_freq'] = 1 # multiples of timestep
    run_config['tavg_freq'] = 1 # multiples of timestep
    
else:
    run_config['inst_freq'] = 12 # multiples of hours
    run_config['tavg_freq'] = 12 # multiples of hours


#---------specify time averaged fields------#
# NOTE: many more options available see mitgcm docs
diag_fields_avg = [['THETA','SALT','UVEL','WVEL','VVEL'],['UVELSLT ','UVELTH  ','WVELSLT ','WVELTH  '],['BRGfwFlx','BRGhtFlx','BRGmltRt']]
diag_fields_max = 0
diag_fields_avg_name = ['dynDiag','fluxDiag','BRGFlx']
# diag_fields_avg = ['UVEL', 'VVEL', 'WVEL', 'UVELSQ', 'VVELSQ', 'WVELSQ',
#                   'UVELTH', 'VVELTH', 'WVELTH', 'THETA', 'THETASQ',
#                   'PHIHYD', 'LaUH1TH', 'LaVH1TH', 'LaHw1TH','LaHs1TH']

numdiags_avg = len(diag_fields_avg)
numdiags_avg_total = 0
for i in range(len(diag_fields_avg)):
    numdiags_avg_total += len(diag_fields_avg[i])
diag_phase_avg = 0.0

if run_config['test'] == True:
    diag_freq_inst = -run_config['inst_freq']*deltaT # negative values indicate snapshots at given interval
    diag_freq_avg = run_config['tavg_freq']*deltaT # positive values indicate time average over specified interval
else:
    diag_freq_inst = -run_config['inst_freq']*3600
    diag_freq_avg = run_config['tavg_freq']*3600
    
    
diag_params01 = {}
diag_params01['diag_mnc'] = False #<---you would need to modify this if you want netcdf output

for ii in range(numdiags_avg):  
    n = ii+1
    if len(diag_fields_avg[ii]) > diag_fields_max:
        diag_fields_max = len(diag_fields_avg[ii])
    diag_params01['fields(1:%i,%s)'%(len(diag_fields_avg[ii]),n)] ="','".join(diag_fields_avg[ii])
    diag_params01['fileName(%s)'%n] = diag_fields_avg_name[ii]
    diag_params01['frequency(%s)'%n] = diag_freq_avg
    diag_params01['timePhase(%s)'%n] = diag_phase_avg

    
#--------specify instanteous fields (i.e. snapshots)--------#
diag_fields_inst = [['UVEL', 'VVEL', 'WVEL','THETA','SALT']]
diag_fields_names = ['dynDiag']
numdiags_inst = len(diag_fields_inst)
diag_phase_inst = 0.0

for ii in range(numdiags_inst):
    n = numdiags_avg+ii+1
    if len(diag_fields_inst[ii]) > diag_fields_max:
        diag_fields_max = len(diag_fields_inst[ii])
    diag_params01['fields(1:%i,%s)'%(len(diag_fields_inst[ii]),n)] = "','".join(diag_fields_inst[ii])
    diag_params01['fileName(%s)'%n] = diag_fields_names[ii] + '_inst'
    diag_params01['frequency(%s)'%n] = diag_freq_inst
    diag_params01['timePhase(%s)'%n] = diag_phase_inst
    
print(diag_params01)
setupNotes.write('Diagnostic Settings\n')
setupNotes.write(str(diag_params01)+'\n')
Ndiags = n


# In[19]:


diag_params02={}
diag_params = [diag_params01, diag_params02]
rcf.write_data(run_config, diag_params, group_name='diagnostics')


# In[20]:


## create DIAGNOSTICS_SIZE.h
Nlevels = grid_params['Nr']
rcf.createDIAGSIZEh(run_config, Ndiags, Nlevels)


# ## Generate eedata 

# In[21]:


# create eedata
rcf.create_eedata(run_config, grid_params['nTx'], grid_params['nTy'])


# # Generate Input Files

writeFiles = True

def write_bin(fname, data):
    print(fname, np.shape(data))
    setupNotes.write(fname + " " + str(np.shape(data))+'\n')
    if(writeFiles):
        setupNotes.write(fname + " " + str(np.shape(data))+'\n')
        data.astype(">f8").tofile(run_config['run_dir']+'/input/'+fname)
    else:
        setupNotes.write('Not saving\n')
        print('Not saving')


gravity = 9.81
sbeta = 8.0e-4
talpha = 0.4e-4
rho0 = 999.8
T0 = 1
S0 = 34

x = np.zeros([grid_params['Ny'], grid_params['Nx']])
x[:, 0] = run_config['horiz_res_m'] / 2

for i in np.arange(1, grid_params['Nx']):
    x[:,i] = x[:, i - 1] + run_config['horiz_res_m']

z = -np.cumsum(dz)


# In[25]:


# Temperature profile
tcd = 40
Tmin = -1
Tmax = 1
Tc = (Tmax + Tmin) / 2
Trange = Tmax - Tmin
T2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
for j in np.arange(0,grid_params['Ny']):
    T2[:,j] = Tc - Trange / 2 * np.tanh(np.pi * (z + tcd) / tcd)
Tconst = np.zeros([grid_params['Nr'],grid_params['Ny']]) + 0.4
T = T2

Sc = 33.5
Srange = -1
S2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
for j in np.arange(0,grid_params['Ny']):
    S2[:,j] = Sc + Srange / 2 * np.tanh(np.pi * (z + tcd) / tcd)
Sconst = np.zeros([grid_params['Nr'],grid_params['Ny']]) + 35

S = S2

Rref = rho0 * (1 - talpha * (T - T0) + sbeta * (S - S0))

t = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
s = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])

for j in np.arange(0,grid_params['Ny']):
    for k in np.arange(0, grid_params['Nr']):
        t[k, j, :] = t[k, j, :] + T[k,j]
        s[k, j, :] = s[k, j, :] + S[k,j]

write_bin("T.init", t)
write_bin("S.init", s)
write_bin("S.bound", S)
write_bin("T.bound", T)

plt.plot(S - 34, z, 'b', label="Sref - 34")
plt.plot(T, z, 'r', label="Tref")
plt.legend()
plt.savefig("%sinitialTS" % (run_config['run_dir']+'/input/'))
plt.show()


# In[26]:


# Topography
d = np.zeros([grid_params['Ny'], grid_params['Nx']]) - domain_params['H']
d[ 0, :] = 0  # walls of fjord
d[-1, :] = 0
d[: , 0] = 0 #cap west side
write_bin("topog.slope", d)

# Plume
nt = 1 #if variable forcing
runoffVel = np.zeros([grid_params['Ny'],grid_params['Nx'],nt])
runoffRad = np.zeros([grid_params['Ny'],grid_params['Nx'],nt])
plumeMask = np.zeros([grid_params['Ny'],grid_params['Nx']])

# Total runoff (m^3/s)
runoff = 500 

# velocity (m/s) of subglacial runoff
wsg = 1

# ice front location
icefront=1 # adjacent to wall at western end of domain, simulate wall of ice

# plume location
plume_loc = 5

## Define plume-type mask 
# 1 = ice but no plume (melting only)
# 2 = sheet plume (Jenkins)
# 3 = half-conical plume (Morton/Slater)
# 4 = both sheet plume and half-conical plume (NOT YET IMPLEMENTED)
# 5 = detaching conical plume (Goldberg)
# POSITIVE values indicate ice front is orientated north-south
# NEGATIVE values indicate ice front is orientated east-west

# Create virtual ice wall
plumeMask[1:-1,icefront] = 1 
# Located 1 cell in from western boundary (need solid barrier behind), and extending across the fjord with (fjord walls either side)

# Specify discharge location
plumeMask[plume_loc,icefront] = 3 # runoff emerges from centre of grounding line

# specify a runoff velocity of 1 m/s
runoffVel[plume_loc,icefront,:] = wsg

# calculate channel radius
runoffRad[plume_loc,icefront,:] = np.sqrt(2*runoff/(np.pi*wsg))

# Write files
write_bin("runoffVel.bin", runoffVel)
write_bin("runoffRad.bin", runoffRad)
write_bin("plumeMask.bin", plumeMask)

## Boundary conditions

# pre-allocate
EBCu = np.zeros([grid_params['Ny'],grid_params['Nx'],nt])

# Apply barotropic velocity to balance input of runoff
if runoff > 0:
    fjordMouthCrossSection = -np.sum(d[:,-1]) * run_config['horiz_res_m']
    fjordMouthVelocity = runoff/fjordMouthCrossSection
    # Out-of-domain velocity is positive at eastern boundary
    EBCu[:] = fjordMouthVelocity

write_bin("EBCu.bin", EBCu)


# Pace (Stanford)


cluster_params = {}
cluster_params['cluster_name'] = 'PACE'
cluster_params['opt_file'] = 'darwin_amd64_gfortran' #<-- may need to update this at some point
cluster_params['mvapich2_ver'] = '2.3.7' #'4.0.3' 4.1.2'
cluster_params['mvapich2_inc_dir'] = '/usr/local/pace-apps/spack/packages/linux-rhel9-x86_64_v3/gcc-12.3.0/mvapich2-2.3.7-1-qv3gjagtbx5e3rlbdy6iy2sfczryftyt/' 
cluster_params['netcdf_dir'] = ''
cluster_params['use_mpi'] = True

cluster_params['email'] = email

cluster_params['exps_dir'] = run_config['run_dir']
cluster_params['run_dir'] = os.path.join(cluster_params['exps_dir'], run_config['run_name'])
cluster_params['cpus_per_node'] = 10 





# ## Estimate wall clock time
ncpus = run_config['ncpus_xy'][0]*run_config['ncpus_xy'][1]
print('===== Wall Clock Time =====')
estTime = int(grid_params['Ny']) * int(grid_params['Nx']) * int(grid_params['Nr']) * int(grid_params['Nt']) *2e-7
print('Estimated run time is %.2f hours for one CPU' % (estTime/60))
print('Estimated run time is %.2f hours for %i CPUs\n' % (estTime/60/ncpus*1.2,ncpus))

setupNotes.write('===== Wall Clock Time =====\n')
setupNotes.write('Estimated run time is %.2f hours for one CPU\n' % (estTime/60))
setupNotes.write('Estimated run time is %.2f hours for %i CPUs\n' % (estTime/60/ncpus*1.2,ncpus))
comptime_hrs = estTime/60/ncpus*1.2 
setupNotes.close()

shutil.move('setupReport.txt', run_config['run_dir'])


rcf.createSBATCHfile_Sherlock(run_config, cluster_params, walltime_hrs=1.2*comptime_hrs, email=email, mem_GB=1)

print('Done! Remember to make MATLAB files and build before you run the script, building on MPI time is very inefficient')


# ## Next steps
# 
# Once you've successfully set up the model experiment, you will need to do the following:
# 
# + Navigate to the build directory and compile the code (i.e., `bash ../makeBuild ../../..`). This may take a couple minutes and will generates lots of new files in the build directory. The key file is the `mitgcmuv` executable.
# 
# 
# 
# + Next, go to the results directory and submit the job to Sherlock queue or run directly if a test (`bash ../makeRun`). It is recommended that you do this from a login node. This should generate a unique output_xxx file, which logs error messages for each run. You will get emails when the job is released from the queue and when it is completed. Test runs normally take a few minutes. A successful run will generate dozens of output files, including data files like `THETA_inst_000000.data`. 
# 
# 
# To check if the job is running:   `squeue | grep ejdawson`
# 

# In[ ]:




