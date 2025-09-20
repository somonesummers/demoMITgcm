#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import matplotlib.pyplot as plt
import pickle
import importlib
import shutil
import platform
import fileinput
import sys
import glob
import cmocean

OSX = platform.system()

from scipy.interpolate import make_interp_spline

import sys
if OSX == 'Darwin':
    sys.path.append('/Users/psummers8/Documents/MITgcm/MITgcm/elizaScripts/main_scripts')
else:
    sys.path.append('/storage/home/hcoda1/2/psummers8/MITgcmSandbox/elizaScripts/main_scripts')
import build_domain_funcs as build_domain 
import run_config_funcs as rcf # import helpter functions

#Set up new folder, copy over files from DEFAULT directory, will overwrite if existing experiment
makeDirs = True
#Write input files, this lets us update the inputs with a full new run
writeFiles = True

if(makeDirs):
    setupNotes = open("setupReport.txt", "w") 

def setUpPrint(msg):
    print(msg)
    if(makeDirs):
        setupNotes.write(str(msg) + "\n")

# ## Main run configuration
email = 'psummers8@gatech.edu'
# set high level run configurations

setUpPrint('\n====== Welcome to the mélange building script =====')
setUpPrint('\tMaking experiment to compare mélange realizations')
#========================================================================================
#main values to imput 

run_config = {}
grid_params = {}
run_config['ncpus_xy'] = [5, 1] # cpu distribution in the x and y directions
run_config['run_name'] = 'solidTest_m0_b1'
run_config['ndays'] = 5 # simulaton time (days)
run_config['test'] = False # if True, run_config['nyrs'] will be shortened to a few time steps

run_config['horiz_res_m'] = 250 # horizontal grid spacing (m)
run_config['Lx_m'] = 40000 # domain size in x (m)
run_config['Ly_m'] = 5000 + 2 * run_config['horiz_res_m'] # domain size in y (m) with walls
# NOTE: the number of grid points in x and y should be multiples of the number of cpus.

grid_params['Nr'] = 80 # num of z-grid points

# Offshore current =========================
oscStrength = .3 #[m/s] peak strength of offshore current
lengthOffShoreCurrent = 5e3 #width of offshore current [m]
indexOSC = int(lengthOffShoreCurrent/run_config['horiz_res_m'])

# Iceberg configuration =========================
iceBergDepth = 200 # max iceberg depth [meters], used for ICEBERG package
iceExtent = 5000 # [meters] of extent of ice
iceCoverage = 80 # % of ice cover in melange, stay under 97% probably (below hFacMin means no barrier effect)

#========================================================================================
# The rest of this should take care of it self mostly

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
setUpPrint('run_config is %s' %run_config)

#========================================================================================
# Generate new experiment directory and copy over defaults

# create experimentary directory on SCRATCH and copy over default configuration
# NOTE: this step does not overwrite existing directories. 
if(makeDirs):
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
    setUpPrint('run directory and subdirectories:')
    setUpPrint(run_config['run_dir'])
    setUpPrint(os.listdir(run_config['run_dir']))

    # create new analysis sub-dir in your home directory
    # if OSX == 'Darwin':
    #     analysis_dir = '/Users/psummers8/Documents/MITgcm/MITgcm/analysis/%s'%run_config['run_name']
    # else:
    #     analysis_dir = '/storage/home/hcoda1/2/psummers8/MITgcmSandbox/analysis/%s'%run_config['run_name']
    # os.makedirs(analysis_dir, exist_ok=True)
    
secsInDay = 24*60*60
secsInYear = 365*secsInDay


#========================================================================================
# set domain size
setUpPrint('====== Domain Size and Parameters =====')
domain_params = {}
domain_params['Lx'] = run_config['Lx_m'] # domain size in x (m)
domain_params['Ly'] = run_config['Ly_m'] # domain size in y (m)
domain_params['L_sponge'] = 20000 # width of eastern sponge layer (m)
domain_params['H'] = 500 # max domain depth (m)

# NOTE: the only thing you may need to change here is the number of z-grid pointsm, which was set above)

grid_params['nSx'] = 1 # num of tiles per processor in x-direction
grid_params['nSy'] = 1 # num of tiles per processor in y-direction
grid_params['nTx'] = 1 # num of threads per processor in x-direction
grid_params['nTy'] = 1 # num of threads per processor in y-direction
grid_params['OLx'] = 3 # num of overlapping x-gridpoints per tile
grid_params['OLy'] = 3 # num of overlapping y-gridpoints per tile
#grid_params['Nr'] = 70 # num of z-grid points (set above)

grid_params['nPx'] = run_config['ncpus_xy'][0] #num of processors in x-direction
grid_params['nPy'] = run_config['ncpus_xy'][1] #num of processors in x-direction

# grid_params['nSx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points in sub grid
# grid_params['nSy'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points in sub grid

# grid_params['Nx'] = grid_params['sNx'] * grid_params['nSx'] * grid_params['nPx']
# grid_params['Ny'] = grid_params['sNy'] * grid_params['nSy'] * grid_params['nPy']

grid_params['Nx'] = domain_params['Lx']/(run_config['horiz_res_m']) # num of x points
grid_params['Ny'] = domain_params['Ly']/(run_config['horiz_res_m']) # num of y points

setUpPrint("Nx: %s" %grid_params['Nx'])
setUpPrint("Ny: %s" %grid_params['Ny'])

grid_params['sNx'] = grid_params['Nx']/grid_params['nPx']#num of x-gridpoints per tile
grid_params['sNy'] = grid_params['Ny']/grid_params['nPy'] #num of y-gridpoints per tile

setUpPrint("sNx: %s" %grid_params['sNx'])
setUpPrint("sNy: %s" %grid_params['sNy'])

# NOTE: sNx and sNy should be whole numbers/integers. As long we keep the horizontal resolution,
# domain dimesions, and number of cpus to be multiples of five, we should be ok. 

for key, param  in grid_params.items():
    assert param%1==0, "grid parameter needs to be an integer"
    grid_params[key] = int(param)

setUpPrint('Grid parameters')    
setUpPrint(grid_params)

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

#========================================================================================
#Physical parameters

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
#params01['viscA4'] = 0.0000 # Biharmonic viscosity?
params01['viscAz'] = 1.0e-2 # Vertical viscosity
#params01['viscAh'] = 2.5e-1 # Vertical viscosity
params01['viscC2smag'] = 2.2 # ??? viscosity

# advection and time stepping
params01['tempAdvScheme'] = 33 # needs to be int
params01['saltAdvScheme'] = 33 # needs to be int
#params01['tempStepping'] = True
#params01['saltStepping'] = run_config['evolve_salt']
params01['staggerTimeStep'] = True

# diffusivity
#params01['diffK4T'] = 0.0e4 # ?? temp diffusion
params01['diffKhT'] = 1.0 # Horizontal temp diffusion
params01['diffKhS'] = 1.0 # Horz salt diffusion
params01['diffKzT'] = 4.0e-3 # Vertical temp diffusion
params01['diffKzS'] = 4.0e-3 # Vert salt diffusion
#params01['diffK4S'] = 0.0e4 # ?? salt diffusion


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


# ## Check for numericl stability?

# ## Numeric solvers and I/O controls

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
params03['pChkptFreq'] = 864000.0
params03['taveFreq'] = 0.0
params03['dumpFreq'] = 864000.0
params03['taveFreq'] = 0.0
params03['monitorFreq'] = 86400.0
params03['monitorSelect'] = 1


params03['periodicExternalForcing'] = False
params03['ExternForcingPeriod'] = 100.0
params03['ExternForcingCycle'] = 1000.0 

if run_config['test']:
    nTimeSteps = 10
else:
    nTimeSteps = np.ceil(run_config['ndays']*secsInDay/deltaT)

simTimeAct = nTimeSteps*deltaT

params03['endTime'] = int(params03['nIter0']*deltaT+simTimeAct)
params03['deltaT'] = np.round(deltaT)
grid_params['Nt'] = nTimeSteps

# gather params for data file 
params04 = {} #<-- using params04 to be consistent with ordering in Andrew's code
params04['usingCartesianGrid'] = grid_params['usingCartesianGrid']
params04['usingSphericalPolarGrid'] = grid_params['usingSphericalPolarGrid']
params04['delX']  = grid_params['delX']
params04['delY'] = grid_params['delY']
params04['delZ'] = dz


# get data fnames param
params05 = {}
params05['bathyFile'] ='bathymetry.bin'
params05['hydrogThetaFile'] = 'T.init'
params05['hydrogSaltFile'] = 'S.init'

if(makeDirs):
    data_params = [params01, params02, params03, params04, params05]
    #data file
    rcf.write_data(run_config, data_params, group_name='data', lf=run_config['lf'])
    #SIZE file
    rcf.createSIZEh(run_config, grid_params)

#========================================================================================
# Diagnostics

# adjust output frequency
if run_config['test']:
    run_config['inst_freq'] = 1 # multiples of timestep
    run_config['tavg_freq'] = 1 # multiples of timestep
    
else:
    run_config['inst_freq'] = 12 # multiples of hours
    run_config['tavg_freq'] = 1 # multiples of hours


#---------specify time averaged fields------#
# NOTE: many more options available see mitgcm docs
diag_fields_avg = [['THETA','SALT','UVEL','WVEL','VVEL'],
    ['UVELSLT ','UVELTH  ','WVELSLT ','WVELTH  '],
    ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY'],
    ['icefrntW','icefrntT','icefrntS','icefrntR','icefrntM']]
diag_fields_max = 0
diag_fields_avg_name = ['dynDiag','fluxDiag','BRGFlx','plumeDiag']
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
diag_fields_inst = [['THETA','SALT','UVEL','WVEL','VVEL']]
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

setUpPrint('Diagnostic Settings')
setUpPrint(diag_params01)

Ndiags = n

diag_params02={}
diag_params = [diag_params01, diag_params02]


#========================================================================================
# Boundary Conditions
setUpPrint('====== Boundary Conditions =====')

obcs_params01 = {}
obcs_params02 = {}
obcs_params03 = {}

obcs_params01['OB_singleIeast'] = -1
obcs_params01['OB_Jsouth(%i:%i)'%(grid_params['Nx']-indexOSC+1,grid_params['Nx'])] = np.ones(indexOSC,dtype=int)
obcs_params01['OB_Jnorth(%i:%i)'%(grid_params['Nx']-indexOSC+1,grid_params['Nx'])] = -1*np.ones(indexOSC,dtype=int)
obcs_params01['useOBCSsponge'] = False
obcs_params01['useOBCSprescribe']= True
#East
obcs_params01['OBEsFile']='EBCs.bin'
obcs_params01['OBEtFile']='EBCt.bin'
obcs_params01 ['OBEvFile']='EBCv.bin'
#North
obcs_params01['OBNsFile']='NsBCs.bin'  
obcs_params01['OBNtFile']='NsBCt.bin'  
obcs_params01 ['OBNvFile']='NsBCv.bin'
#South
obcs_params01['OBSsFile']='NsBCs.bin'
obcs_params01['OBStFile']='NsBCt.bin'
obcs_params01 ['OBSvFile']='NsBCv.bin'

obcs_params03['spongeThickness'] = int(domain_params['L_sponge'] / run_config['horiz_res_m']) #grid cells
obcs_params03['Urelaxobcsinner'] = 86400.0
obcs_params03['Urelaxobcsbound'] = 3600.0
obcs_params03['Vrelaxobcsinner'] = 86400.0
obcs_params03['Vrelaxobcsbound'] = 3600.0
obcs_params = [obcs_params01, obcs_params02, obcs_params03]
if(makeDirs):
    rcf.write_data(run_config, diag_params, group_name='diagnostics')
    
    ## create DIAGNOSTICS_SIZE.h
    Nlevels = grid_params['Nr']
    rcf.createDIAGSIZEh(run_config, Ndiags, Nlevels)

    # create eedata
    rcf.create_eedata(run_config, grid_params['nTx'], grid_params['nTy'])

    #create data.obcs
    rcf.write_data(run_config, obcs_params, group_name='obcs')

#========================================================================================
#Domain initialization and saving
setUpPrint('====== Domain Initialization =====')

def write_bin(fname, data):
    setUpPrint(fname + str(np.shape(data)))
    if(writeFiles):
        data.astype(">f8").tofile(run_config['run_dir']+'/input/'+fname)
    else:
        setUpPrint('\tNot saving')
#Similar params as fed into MITgcm, but redeclared here
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

y = np.zeros([grid_params['Ny'], grid_params['Nx']])
y[0, :] = run_config['horiz_res_m'] / 2

for j in np.arange(1, grid_params['Ny']):
    y[j,:] = y[j-1, :] + run_config['horiz_res_m']

z = -np.cumsum(dz)


# Topography
sillStart = 45000
sillHeight = 0  #no sill
sillLength = 5000
fjordEnd = int(grid_params['Nx'] - indexOSC)

d = np.zeros([grid_params['Ny'], grid_params['Nx']]) - domain_params['H']
d[x>sillStart] = - domain_params['H'] + (x[x>sillStart]-sillStart) * sillHeight/sillLength
d[x>(sillStart + sillLength)] = sillHeight - domain_params['H']
setUpPrint('fjord end: %i' %fjordEnd)
d[ 0, 1:fjordEnd] = 0  # walls of fjord
d[-1, 1:fjordEnd] = 0
d[: , 0] = 0 #cap west side

plt.figure
plt.plot(x[10,:],d[10,:])
plt.pcolormesh(x,y,d)
plt.colorbar()
if(writeFiles):
    plt.savefig("%sbathymetry" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

write_bin("bathymetry.bin", d)

# Temp/Salt/Vel  initial/boundaries
from scipy import interpolate
t2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
s2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
S2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
T2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
S_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
T_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
V_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
W_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])

# Temp/Salt/Vel  initial/boundaries
from scipy import interpolate
t2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
s2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
S2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
T2 = np.zeros([grid_params['Nr'],grid_params['Ny']])
S_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
T_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
V_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])
W_ns = np.zeros([grid_params['Nr'],(grid_params['Nx'])])

z_tmp =  np.asarray([  0,  10,   50,  100,  200,  300, 500]); #must be increasing, so do depth as positive, see negs later for z[:]
t_tmp =  np.asarray([  1,   1,  1.5,  1.8,  2.1,  2.3, 2.6]);
s_tmp =  np.asarray([ 33,33.2, 33.8, 34.0, 34.3, 34.4,34.6]);
t_int = interpolate.PchipInterpolator(z_tmp, t_tmp)
s_int = interpolate.PchipInterpolator(z_tmp, s_tmp)
for j in np.arange(0,grid_params['Ny']):
    for i in np.arange(0, grid_params['Nx']):
        t2[:, j, i] = t_int(-1 * z[:])
        s2[:, j, i] = s_int(-1 * z[:])

#East BC
for j in np.arange(0,grid_params['Ny']):
    T2[:,j] = t_int(-1 * z[:])
    S2[:,j] = s_int(-1 * z[:])

#BC for V at East side
Ve = np.zeros([grid_params['Nr'],grid_params['Ny']])
Ve[:,:] = oscStrength #[m/s]

#N/S BCs
for i in np.arange(fjordEnd,grid_params['Nx']):
    T_ns[:,i] = t_int(-1 * z[:])
    S_ns[:,i] = s_int(-1 * z[:])
    V_ns[:,i] = oscStrength * (i-fjordEnd)/indexOSC #[m/s] along coast flow


write_bin("T.init", t2)
write_bin("S.init", s2)
write_bin("EBCs.bin", S2)
write_bin("EBCt.bin", T2)
write_bin("EBCs.bin", S2)
write_bin("EBCt.bin", T2)
write_bin("EBCv.bin", Ve)
write_bin("NsBCs.bin", S_ns)
write_bin("NsBCt.bin", T_ns)
write_bin("NsBCv.bin", V_ns)
write_bin("NsBCW.bin", W_ns)

plt.figure()
plt.plot(S2[:,0] - 34, z, 'b', label="Sref - 34")
plt.plot(T2[:,0], z, 'r', label="Tref")
# plt.scatter(s_tmp - 34,-z_tmp,color='b')
# plt.scatter(t_tmp,-z_tmp,color='r')
plt.legend()
if(writeFiles):
    plt.savefig("%sinitialTS" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

# Plume
nt = 1 #if variable forcing
runoffVel = np.zeros([grid_params['Ny'],grid_params['Nx'],nt])
runoffRad = np.zeros([grid_params['Ny'],grid_params['Nx'],nt])
plumeMask = np.zeros([grid_params['Ny'],grid_params['Nx']])

# Total runoff (m^3/s)
runoff = 100

# velocity (m/s) of subglacial runoff
wsg = 1

# ice front location
icefront=1 # adjacent to wall at western end of domain, simulate wall of ice

# plume location
plume_loc = int(np.round(grid_params['Ny']/2))
setUpPrint('Plume Location: %i discharge: %f' %(plume_loc, runoff))
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

plt.figure(1)
plt.pcolormesh(x,y,plumeMask)
plt.colorbar()
if(writeFiles):
    plt.savefig("%splumeMask" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

## Boundary conditions

# pre-allocate
EBCu = np.zeros([grid_params['Nr'],grid_params['Ny']])

# Apply barotropic velocity to balance input of runoff
if runoff > 0:
    fjordMouthCrossSection = -np.sum(d[:,-1]) * run_config['horiz_res_m']
    fjordMouthVelocity = runoff/fjordMouthCrossSection
    # Out-of-domain velocity is positive at eastern boundary
    EBCu[:] = fjordMouthVelocity

write_bin("EBCu.bin", EBCu)

#=======================================================================================
# Make Bergs, now all in python
setUpPrint('====== Making mélange =====')
#Make Masks

nz = grid_params['Nr']
ny = grid_params['Ny']
nx = grid_params['Nx']

deltaZ = dz[0]
deltaY = run_config['horiz_res_m']
deltaX = run_config['horiz_res_m']

bergMask = np.zeros([ny,nx])
driftMask = np.zeros([ny,nx])
meltMask = np.zeros([ny,nx])
barrierMask = np.zeros([ny,nx])
bergConc = np.zeros([ny,nx])
bergMaskNums = np.zeros([ny,nx])
numBergsPerCell = np.zeros([ny,nx],dtype=np.int64)

# Berg parameters
bergType = 1 # 1 = block 2 = cone (not implemented)
alpha = 1.8 # slope of inverse power law size frequency distribution
scaling = 2 # 1 = Sulak 2017 2 = Barker 2004
maxBergDepth = iceBergDepth # (m) - set to zero if 'prescribing' max iceberg width, set at top here
minBergDepth= 20 # (m)
maxBergWidth = 0 # (m) - set to zero if 'prescribing' max iceberg depth
minBergWidth = 20 # (m)

iceExtentIndex = int(np.round(iceExtent/run_config['horiz_res_m']))

# Iceberg mask
bergMask[1:-1,1:iceExtentIndex:1] = 1 # icebergs in inner 5 km, all oriented east-west

# Drift mask, No drift for Melange experiments, but can toggle on here if you want
# driftMask[1:-1,1:iceExtentIndex] = 1 # calculate effect of iceberg drift on melt rates 

# Melt mask, only let bergs melt in this region (make melt water, these don't change size)
meltMask[1:-1,1:iceExtentIndex:1] = 0 # Allow focus on blocking effect only

# Barrier mask
barrierMask[1:-1,1:iceExtentIndex:1] = 1 # make icebergs a physical barrier to water flow
barrierMask[plume_loc,icefront] = 0 #Plume code 

# Iceberg concentration (# of each surface cell that is filled in plan view)
bergConc[1:-1,1:iceExtentIndex:1] = iceCoverage # iceberg concentration set at top

desiredBergArea = np.sum(bergConc/100.0*deltaX*deltaY)
bergMaskArea = np.sum(bergMask*deltaX*deltaY)
setUpPrint('Area where bergs live: ' + str(bergMaskArea) + ' m^2')
setUpPrint('Desired berg area: ' + str(desiredBergArea) + ' m^2')
setUpPrint('Ratio: ' + str(desiredBergArea/bergMaskArea*100) + '%')

setUpPrint('Making bergs to exact size')

bergMaski = 0  #Bad name, but this is the count of cells that will recieve bergs
bergDict = {}

for j in range(ny):
    for i in range(nx):
        if(bergMask[j,i] == 1):
            # print('i,j, bergmask',i,j,bergMask[j,i])
            bergMaski = 1 + bergMaski #Needs to start at 1, as non-bergs will be 0
            bergMaskNums[j,i] = bergMaski #Assign Mask Nums, not random as we'll randomly place bergs in cells
            bergDict[bergMaski] = [j,i] #This lets us do 1-D loops for the whole grid

bergsCell = 25
numberOfBergs = bergMaski * bergsCell
setUpPrint('%i cells with bergs' % bergMaski)

# make bergs, all identical
bergVariation = 50 #[m]
sorted_depth = np.ones(numberOfBergs) * (np.random.normal(maxBergDepth - bergVariation,bergVariation,numberOfBergs))#  STD + 1 is 'max depth'
sorted_depth[sorted_depth < minBergDepth] = minBergDepth + .1 # Code doesn't like clean numbers (ie berg on perfect cell bottom), disallow small/neg bergs 
sorted_width = np.ones(numberOfBergs) * run_config['horiz_res_m'] * np.sqrt(iceCoverage/100/bergsCell) # cannot be 100% full
sorted_length = np.ones(numberOfBergs) * run_config['horiz_res_m'] * np.sqrt(iceCoverage/100/bergsCell) # cannot be 100% full
# assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) # In this script, every berg has a home

# Need to reshuffle as it re-writes this when sorting
# Array for bergs
bergsPerCellLimit = 500
icebergs_depths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_widths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_length = np.zeros([bergMaski,bergsPerCellLimit])  #careful, not plural as to length match

assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) #every Berg has a spot

icebergs_per_cell = np.zeros([bergMaski],dtype=np.int16)
icebergs_area_per_cell = np.zeros([bergMaski])

assignedCell = 0
for i in range(numberOfBergs): 
    j = assignedCell
    bergArea = sorted_width[i] * sorted_length[i]
    loopLimiter = 0
    icebergs_depths[j,icebergs_per_cell[j]] = sorted_depth[i]
    icebergs_widths[j,icebergs_per_cell[j]] = sorted_width[i]
    icebergs_length[j,icebergs_per_cell[j]] = sorted_length[i]
    icebergs_per_cell[j] += 1
    icebergs_area_per_cell[j] += bergArea
    if(icebergs_per_cell[j] == bergsCell): #cell full, go to next
        assignedCell += 1
setUpPrint('Bergs per cell and filled faction at surface for spot check, should be uniform')     
setUpPrint(icebergs_per_cell)
setUpPrint(np.round(icebergs_area_per_cell/(run_config['horiz_res_m']*run_config['horiz_res_m']),2))

# All bergs now sorted 
openFrac = np.zeros([nz,ny,nx])
SA = np.zeros([nz,ny,nx])
SA[:,:,:] = np.nan
cellVolume = deltaX*deltaY*deltaZ

#This loop knows where all bergs are already, different from searching for all bergs across entire grid
for i in range(bergMaski):
    bergCount = icebergs_per_cell[i]
    numBergsPerCell[bergDict[i+1][0],bergDict[i+1][1]] = bergCount
    if(bergCount > 0):
        lengths = icebergs_length[i,icebergs_length[i,:] > 0] #return only non-zeros
        widths = icebergs_widths[i,icebergs_widths[i,:] > 0] #return only non-zeros
        depths = icebergs_depths[i,icebergs_depths[i,:] > 0] #return only non-zeros
        for k in range(nz):
            d_bot = k*deltaZ + deltaZ #bottom of depth bin
            d_top = k*deltaZ
            volume1 = deltaZ * lengths[depths > d_bot] * widths[depths > d_bot]
            SA1 = deltaZ*2*(lengths[depths > d_bot] + widths[depths > d_bot])
            partialFill = (depths < d_bot) & (depths > d_top)
            #partial fill
            volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
            #partial sides
            SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
            #bottom
            SA3 = lengths[partialFill] * widths[partialFill]
            #print(np.sum(volume1), np.sum(volume2))
            openFrac[k,bergDict[i+1][0],bergDict[i+1][1]] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
            SA[k,bergDict[i+1][0],bergDict[i+1][1]] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
    elif(bergCount == 0):
        openFrac[:,bergDict[i+1][0],bergDict[i+1][1]] = 1
        SA[:,bergDict[i+1][0],bergDict[i+1][1]] = 0



# # All bergs now sorted 
# openFrac = np.zeros([nz,ny,nx])
# SA = np.zeros([nz,ny,nx])
# SA[:,:,:] = np.nan
# cellVolume = deltaX*deltaY*deltaZ

# bergsPerCellLimit = 500

# icebergs_depths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
# icebergs_widths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
# icebergs_length2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])

# #Searching for all bergs across entire grid, fill with full iceberg
# for i in range(nx):
#     for j in range(ny):
#         if(bergMask[j,i] == 1):
#             lengths = np.ones(10) * 500 * np.sqrt(.5/10) #return only non-zeros
#             widths = np.ones(10) * 500 * np.sqrt(.5/10) #return only non-zeros
#             depths = np.ones(10) * 200 #return only non-zeros
#             icebergs_depths2D[1:11,j,i] = depths
#             icebergs_widths2D[1:11,j,i] = widths
#             icebergs_length2D[1:11,j,i] = lengths
#             for k in range(nz):
#                 d_bot = k*deltaZ + deltaZ #bottom of depth bin
#                 d_top = k*deltaZ
#                 volume1 = deltaZ * lengths[depths > d_bot] * widths[depths > d_bot]
#                 SA1 = deltaZ*2*(lengths[depths > d_bot] + widths[depths > d_bot])
#                 partialFill = (depths < d_bot) & (depths > d_top)
#                 #partial fill
#                 volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
#                 #partial sides
#                 SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
#                 #bottom
#                 SA3 = lengths[partialFill] * widths[partialFill]
#                 #print(np.sum(volume1), np.sum(volume2))
#                 openFrac[k,j,i] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
#                 SA[k,j,i] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
#             numBergsPerCell[j,i] = 10 #exactly 10

icebergs_depths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
icebergs_widths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
icebergs_length2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])

for k in range(bergMaski):
    j = bergDict[k+1][0]
    i = bergDict[k+1][1]
    icebergs_depths2D[:,j,i] = icebergs_depths[k,:]
    icebergs_widths2D[:,j,i] = icebergs_widths[k,:]
    icebergs_length2D[:,j,i] = icebergs_length[k,:]

# Plots for reference on whats happening berg-wise
fig = plt.figure()
# plt.subplot(2,2,1)
for i in range(bergMaski):
    plt.plot(openFrac[:,bergDict[i+1][0],bergDict[i+1][1]],-np.cumsum(dz),alpha=.5,color='xkcd:gray',linewidth=.5)
plt.plot(np.mean(openFrac[:,bergMask==1],1),-np.cumsum(dz),alpha=1,color='xkcd:black',linewidth=1,linestyle='--',label='Average Bergs')
plt.plot([0,1],[-maxBergDepth,-maxBergDepth],color = 'xkcd:red',linestyle=':', label='Target Max Depth')
plt.plot([1-np.max(bergConc)/100,1-np.max(bergConc)/100],[-nz*deltaZ,0],color = 'xkcd:gray',linestyle=':',label='Target Max Berg Conc')
plt.xlabel('Open Fraction of Cells')
plt.ylabel('Depth [m]')
plt.legend()  #not quite room so off for now

# plt.subplot(2,2,2)
# plt.hist(inversePowerLawDistNumbers_depth,bins = 50)
# plt.ylabel('Count')
# plt.xlabel('Depth [m]')

# plt.subplot(2,2,3)
# plt.hist(inversePowerLawDistNumbers_width,bins = 50)
# plt.ylabel('Count')
# plt.xlabel('Width [m]')

# plt.subplot(2,2,4)
# plt.hist(inversePowerLawDistNumbers_length,bins = 50)
# plt.ylabel('Count')
# plt.xlabel('Length [m]')
# fig.tight_layout()
plt.savefig(run_config['run_dir']+'/input/bergStatistics.png', format='png', dpi=100)
plt.show()

plt.figure()
pc = plt.pcolormesh(np.nanmax(icebergs_depths2D,0),cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.title('Max Iceberg Depths')
plt.show()

fig = plt.figure()
pltHelper = 1-openFrac[0,:,:]
pltHelper[bergMask == 0] = np.nan
plt.subplot(211)
pc = plt.pcolormesh(pltHelper,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Iceberg Cover and Residual')
plt.ylabel('Cell across fjord')
cbar.set_label('iceberg cover')

plt.subplot(212)
pc = plt.pcolormesh(pltHelper - bergConc/100,cmap='cmo.ice')
cbar = plt.colorbar(pc)
pc_min = np.nanmin(pltHelper) * 100
pc_max = np.nanmax(pltHelper) * 100
plt.xlabel("Cell along fjord, Ice coverage is %.2f%% - %.2f%%" % (pc_min, pc_max))
plt.ylabel('Cell across fjord')
cbar.set_label('cover resid')
fig.tight_layout()
plt.savefig(run_config['run_dir']+'/input/bergMap.png', format='png', dpi=200)
plt.show()

fig = plt.figure()
pc = plt.pcolormesh(meltMask,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Melt Mask')
plt.ylabel('Cell across fjord')
plt.xlabel("Cell along fjord")
fig.tight_layout()
plt.savefig(run_config['run_dir']+'/input/meltMask.png', format='png', dpi=200)
plt.show()

## write iceberg txt files
# setUpPrint('Saving text files for bergs...')
# if(writeFiles):
#     for i in range(bergMaski):
#         with open(run_config['run_dir'] + '/input/iceberg_depth_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_depths[i,icebergs_depths[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))
#         with open(run_config['run_dir']+'/input/iceberg_width_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_widths[i,icebergs_widths[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))
#         with open(run_config['run_dir'] + '/input/iceberg_length_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_length[i,icebergs_length[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))


# write global files
write_bin('bergMask.bin',bergMask)
write_bin('bergMaskNums.bin',bergMaskNums)
write_bin('numBergsPerCell.bin',numBergsPerCell)
write_bin('openFrac.bin',openFrac)
write_bin('totalBergArea.bin',SA)
write_bin('meltMask.bin',meltMask)
write_bin('driftMask.bin',driftMask)
write_bin('barrierMask.bin',barrierMask)
write_bin('icebergs_depths.bin',icebergs_depths2D)
write_bin('icebergs_widths.bin',icebergs_widths2D)
write_bin('icebergs_length.bin',icebergs_length2D)



setUpPrint('Berg setup is done.')


#========================================================================================
#Update files in INPUT directory
setUpPrint('====== Cleaning up input/ files =====')
# a bit of a hack to adjust default values just for this exp

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

#turn on ICEBERG if off
replaceAll(run_config['run_dir'] + '/input/data.pkg','ICEBERG=.FALSE.', 'ICEBERG=.TRUE.') 



#========================================================================================
# PACE (GaTech) 
setUpPrint('====== sbatch script and settings =====')

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

#extra run commands for the sbatch script
extraList = []

run_config['extraCommands'] = "".join(extraList)
     

# ## Estimate wall clock time
ncpus = run_config['ncpus_xy'][0]*run_config['ncpus_xy'][1]
setUpPrint('===== Wall Clock Time =====')
estTime = int(grid_params['Ny']) * int(grid_params['Nx']) * int(grid_params['Nr']) * int(grid_params['Nt']) *2e-7
setUpPrint('Estimated run time is %.2f hours for one CPU' % (estTime/60))
setUpPrint('Estimated run time is %.2f hours for %i CPUs\n' % (estTime/60/ncpus*1.2,ncpus))

comptime_hrs = estTime/60/ncpus*1.2 
if(makeDirs):
    setupNotes.close()

if(makeDirs):
    shutil.move('setupReport.txt', run_config['run_dir']+'/input')
    rcf.createSBATCHfile_Sherlock(run_config, cluster_params, walltime_hrs=1.2*comptime_hrs, email=email, mem_GB=1)

print('Done! Remember to build before you run the script, building on MPI time is very inefficient')


