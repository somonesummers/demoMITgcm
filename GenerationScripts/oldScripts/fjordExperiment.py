#!/usr/bin/env python
#Paul Summers Feb 2025, used to generate fjord scale paper 

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
from bisect import bisect_left

OSX = platform.system()

from scipy.interpolate import make_interp_spline

import sys
if OSX == 'Darwin':
    sys.path.append('/Users/psummers8/Documents/MITgcm/MITgcm/elizaScripts/main_scripts')
else:
    sys.path.append('/storage/home/hcoda1/2/psummers8/MITgcmSandbox/elizaScripts/main_scripts')
import build_domain_funcs as build_domain 
import run_config_funcs as rcf # import helpter functions

#Set up new folder
makeDirs = True
#Write input files, this lets us update the inputs with a full new run
writeFiles = True

if(makeDirs):
    setupNotes = open("setupReport.txt", "w") 

def find_closest_indices(sorted_A, sorted_B):
    closest_indices = []
    for a in sorted_A:
        pos = bisect_left(sorted_B, a)  # Find position in B where a would fit
        # Compare neighbors to find the closest
        if pos == 0:
            closest_indices.append(0)
        elif pos == len(sorted_B):
            closest_indices.append(len(sorted_B) - 1)
        else:
            before = pos - 1
            after = pos
            closest_indices.append(before if abs(sorted_B[before] - a) <= abs(sorted_B[after] - a) else after)
    return closest_indices

def setUpPrint(msg):
    print(msg)
    if(makeDirs):
        setupNotes.write(str(msg) + "\n")

# ## Main run configuration
email = 'psummers8@gatech.edu'
# set high level run configurations

briefSummaryOfExp = """Fjord Scale experiment showcasing impact of subglacial hdydro, bergs, and blocking"""


setUpPrint('====== Welcome to the mélange building script =====')
setUpPrint(briefSummaryOfExp)
input("Confirm above is accurate before continuing...")
setUpPrint('\tMaking experiment to compare mélange realizations')
#========================================================================================
#main values to imput 

run_config = {}
grid_params = {}
run_config['ncpus_xy'] = [20,1] # cpu distribution in the x and y directions
run_config['run_name'] = 'fjord_b0_19'
run_config['ndays'] = 200.0 # simulaton time (days)
run_config['test'] = False # if True, run_config['nyrs'] will be shortened to a few time steps

run_config['horiz_res_m'] = 400 # horizontal grid spacing (m)
run_config['Lx_m'] = 80000 # domain size in x (m)
run_config['Ly_m'] = 6000 + 2 * run_config['horiz_res_m'] # domain size in y (m) with walls
# NOTE: the number of grid points in x and y should be multiples of the number of cpus.

grid_params['Nr'] = 50 # num of z-grid points

# Offshore current =========================
oscStrength = 0.12 #[m/s] peak strength of sin forcing current
lengthOffShoreCurrent = 5e3 #width of offshore current [m]
indexOSC = int(lengthOffShoreCurrent/run_config['horiz_res_m']) 

# Iceberg configuration =========================
iceBergDepth = 200 # max iceberg depth [meters], used for ICEBERG package
iceExtent = 15000 # [meters] of extent of ice
iceCoverage = 60 # % of ice cover in melange, stay under 90% ideally
doMelt = 1 # do we actually calculate melt (0/1 = no/yes)
doBlock = 0 # do we actually calculate blocking (0/1 = no/yes)
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
domain_params['L_sponge'] = 8000 # width of sponge layers (m)
domain_params['H'] = 600 # max domain depth (m)

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

# dz = np.ones(nz)*deltaZ

dz_tmp = np.linspace(1,5,grid_params['Nr'])
dz = dz_tmp/np.sum(dz_tmp)*domain_params['H'] 
sum_z = np.cumsum(dz)
print("dz: \n",dz)
print("z: \n",sum_z)

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
# params01['vectorInvariantMomentum'] = True

#Note: here and elsewhere, we need to be explicit about floats vs ints. E.g., use 12.0 to represent float and
# 12 for int

# viscosity parameters
#params01['viscA4'] = 0.0000 # Biharmonic viscosity?
params01['viscAz'] = 1.0e-4 # Vertical viscosity
params01['viscAh'] = 1.0e-3 # Vertical viscosity, this limits our timestep a lot
params01['viscC2smag'] = 2.5 # ??? viscosity

# advection and time stepping
params01['tempAdvScheme'] = 33 # needs to be int
params01['saltAdvScheme'] = 33 # needs to be int
#params01['tempStepping'] = True
params01['saltStepping'] = True
params01['staggerTimeStep'] = True

# diffusivity
#params01['diffK4T'] = 0.0e4 # ?? temp diffusion
params01['diffKhT'] = 1.0e-5 # Horizontal temp diffusion
params01['diffKhS'] = 1.0e-5 # Horz salt diffusion
params01['diffKzT'] = 1.0e-5 # Vertical temp diffusion
params01['diffKzS'] = 1.0e-5 # Vert salt diffusion
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

params01['no_slip_sides'] = False
params01['no_slip_bottom'] = False
params01['rigidLid'] = False
params01['implicitFreeSurface'] = True
params01['implicSurfPress'] = 1.0
params01['implicDiv2DFlow'] = 1.0
params01['selectAddFluid'] = 1
# params01['useRealFreshWaterFlux'] = True
params01['exactConserv'] = True
params01['implicitViscosity'] = True
params01['implicitDiffusion'] = True

# physical parameters
params01['f0'] = 1.37e-4
# params01['f0'] = 0
params01['beta'] = 0.0e-13
params01['gravity'] = g

# misc
params01['hFacMin'] = 0.05
params01['nonHydrostatic'] = False
params01['readBinaryPrec'] = 64
params01['useSmag3D'] = True
params01['smag3D_coeff'] = 1e-4

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
deltaT = 15
params03['abEps'] = 0.1

#if run_config['testing']:
    
params03['chkptFreq'] = 0.0
params03['pChkptFreq'] = 864000.0
params03['taveFreq'] = 0.0
params03['dumpFreq'] = 8640000.0
params03['taveFreq'] = 0.0
params03['monitorFreq'] = 864000.0
params03['monitorSelect'] = 1


params03['periodicExternalForcing'] = True
params03['ExternForcingPeriod'] = 200*86400/25
params03['ExternForcingCycle'] = 200*86400 

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

if(writeFiles):
    data_params = [params01, params02, params03, params04, params05]
    #data file
    rcf.write_data(run_config, data_params, group_name='data', lf=run_config['lf'])
    #SIZE file
    rcf.createSIZEh(run_config, grid_params)

#========================================================================================
# Stability Check 
u_char = .5 #[m/s]
S_adv = 2 * (u_char * deltaT)/run_config['horiz_res_m']  # < 0.5 (Courant–Friedrichs–Lewy)
S_in = params01['f0'] * deltaT # < 0.5 (adams bashforth II)
# S_lh = 8 * params01['viscAh'] * deltaT /(run_config['horiz_res_m']**2) # < 0.6 
S_lv = 4 * params01['viscAz'] * deltaT /(dz[0]**2) # < 0.6 

S_adv_brg = 2 * (u_char * deltaT)/(run_config['horiz_res_m']*(1 - iceCoverage/100))  # < 0.5 (Courant–Friedrichs–Lewy)
S_in = params01['f0'] * deltaT # < 0.5 (adams bashforth II)
S_lv_brg = 4 * params01['viscAz'] * deltaT /((dz[0]*(1 - iceCoverage/100))**2) # < 0.6
setUpPrint('====== Stability Check =====')
setUpPrint("S_adv: <0.5 , S_in: <0.5 , S_lv: <0.6")
setUpPrint("S_adv: %.04f, S_in: %.04f, S_lv: %.04f" % (S_adv, S_in, S_lv))
setUpPrint("S_adv: %.04f, S_in: %.04f, S_lv: %.04f" % (S_adv_brg, S_in, S_lv_brg))

#========================================================================================
# Diagnostics

# adjust output frequency
if run_config['test']:
    run_config['inst_freq'] = 1 # multiples of timestep
    run_config['tavg_freq'] = 1 # multiples of timestep
    
else:
    run_config['inst_freq'] = 120 # multiples of hours
    run_config['tavg_freq'] = 120 # multiples of hours


#---------specify time averaged fields------#
# NOTE: many more options available see mitgcm docs
diag_fields_avg = [['THETA','SALT','UVEL','WVEL','VVEL'],
                    ['BRGfwFlx','BRGhtFlx','BRGmltRt','BRG_TauX','BRG_TauY'],
                    ['icefrntW','icefrntT','icefrntS','icefrntA','icefrntM'],
                    ]
diag_fields_max = 0
diag_fields_avg_name = ['dynDiag','BRGFlx','plumeDiag']
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
# diag_fields_inst = [['THETA','SALT','UVEL','WVEL','VVEL']]
# diag_fields_names = ['dynDiag']
# numdiags_inst = len(diag_fields_inst)
# diag_phase_inst = 0.0

# for ii in range(numdiags_inst):
#     n = numdiags_avg+ii+1
#     if len(diag_fields_inst[ii]) > diag_fields_max:
#         diag_fields_max = len(diag_fields_inst[ii])
#     diag_params01['fields(1:%i,%s)'%(len(diag_fields_inst[ii]),n)] = "','".join(diag_fields_inst[ii])
#     diag_params01['fileName(%s)'%n] = diag_fields_names[ii] + '_inst'
#     diag_params01['frequency(%s)'%n] = diag_freq_inst
#     diag_params01['timePhase(%s)'%n] = diag_phase_inst

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

#time varying forcing
nt = 25

def write_bin(fname, data):
    setUpPrint(fname + " " + str(np.shape(data)))
    if(writeFiles):
        data.astype(">f8").tofile(run_config['run_dir']+'/input/'+fname)
    else:
        setUpPrint('Not saving')
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
fjordEnd = int(grid_params['Nx'] - indexOSC)


d = np.zeros([grid_params['Ny'], grid_params['Nx']]) - domain_params['H']
setUpPrint('fjord end: %i' %fjordEnd)
d[ 0, 1:fjordEnd] = 0  # walls of fjord
d[-1, 1:fjordEnd] = 0
d[: , 0] = 0 #cap west side


plt.figure
plt.plot(x[5,:],d[5,:])
plt.pcolormesh(x,y,d)
plt.colorbar()
if(writeFiles):
    plt.savefig("%sbathymetry" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

write_bin("bathymetry.bin", d)

# Temp/Salt/Vel  initial/boundaries nt is for time varying BCs
from scipy import interpolate
t2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
s2 = np.zeros([grid_params['Nr'],grid_params['Ny'],grid_params['Nx']])
S2 = np.zeros([nt,grid_params['Nr'],grid_params['Ny']])
T2 = np.zeros([nt,grid_params['Nr'],grid_params['Ny']])
S_ns = np.zeros([nt,grid_params['Nr'],(grid_params['Nx'])])
T_ns = np.zeros([nt,grid_params['Nr'],(grid_params['Nx'])])
V_ns = np.zeros([nt,grid_params['Nr'],(grid_params['Nx'])])
W_ns = np.zeros([nt,grid_params['Nr'],(grid_params['Nx'])])

z_tmp =  np.asarray([  0,  600]); #must be increasing, so do depth as positive, see negs later for z[:]
t_tmp =  np.asarray([  1,    3]); #const temp
s_tmp =  np.asarray([ 32,   35]); #linear salt
t_int = interpolate.PchipInterpolator(z_tmp, t_tmp)
s_int = interpolate.PchipInterpolator(z_tmp, s_tmp)
for j in np.arange(0,grid_params['Ny']):
    for i in np.arange(0, grid_params['Nx']):
        t2[:, j, i] = t_int(-1 * z[:])
        s2[:, j, i] = s_int(-1 * z[:])

#East BC
for j in np.arange(0,grid_params['Ny']):
    for k in range(nt):
        T2[k,:,j] = t_int(-1 * z[:])
        S2[k,:,j] = s_int(-1 * z[:])

#BC for V at East side
Ve = np.zeros([nt,grid_params['Nr'],grid_params['Ny']])
Ve[:,:,:] = oscStrength #[m/s]

#N/S BCs
for i in np.arange(fjordEnd,grid_params['Nx']):
    for k in range(nt):
        T_ns[k,:,i] = t_int(-1 * z[:])
        S_ns[k,:,i] = s_int(-1 * z[:])
        V_ns[k,:,i] = oscStrength * (i-fjordEnd)/indexOSC #[m/s] along coast flow

write_bin("T.init", t2)
write_bin("S.init", s2)
write_bin("EBCs.bin", S2)
write_bin("EBCt.bin", T2)
write_bin("EBCv.bin", Ve)
write_bin("NsBCs.bin", S_ns)
write_bin("NsBCt.bin", T_ns)
write_bin("NsBCv.bin", V_ns)
write_bin("NsBCW.bin", W_ns)

plt.figure()
plt.plot(s2[:,1,1] - 34, z, 'b', label="Sref - 34")
plt.plot(t2[:,1,1], z, 'r', label="Tref")
plt.scatter(s_tmp - 34,-z_tmp,color='b')
plt.scatter(t_tmp,-z_tmp,color='r')
plt.legend()
if(writeFiles):
    plt.savefig("%sinitialTS" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

#=======================================================================================
# Plume
setUpPrint('====== Plume =====')
runoffVel = np.zeros([nt,grid_params['Ny'],grid_params['Nx']])
runoffRad = np.zeros([nt,grid_params['Ny'],grid_params['Nx']])
plumeMask = np.zeros([grid_params['Ny'],grid_params['Nx']])

# Total runoff (m^3/s)
# runoff = 350 + 150 * np.sin(np.pi * np.arange(nt)/12.5)
runoff = 500 * np.ones(nt)
setUpPrint('Runoff is:')
setUpPrint(runoff)
# velocity (m/s) of subglacial runoff
wsg = 1 

# ice front location
icefront=1 # adjacent to wall at western end of domain, simulate wall of ice

# plume location
plume_loc = int(np.round(grid_params['Ny']/2))
setUpPrint('Plume Location: %i discharge: %f' %(plume_loc, np.nanmax(runoff)))
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
runoffVel[:,plume_loc,icefront] = wsg

# calculate channel radius
runoffRad[:,plume_loc,icefront] = np.sqrt(2*runoff/(np.pi*wsg))

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

plt.figure()
plt.plot(runoffRad[:,plume_loc,icefront],label='radius')
plt.plot(runoffVel[:,plume_loc,icefront]*.5*np.pi*runoffRad[:,plume_loc,icefront]**2,label='Volume')
plt.legend()
if(writeFiles):
    plt.savefig("%sforcingVelocity" % (run_config['run_dir']+'/input/'))
plt.show()
plt.close()

#=======================================================================================
# Make Bergs, now all in python
setUpPrint('====== Making mélange =====')
#Make Masks

hfacThreshold = .95

nz = grid_params['Nr']
ny = grid_params['Ny']
nx = grid_params['Nx']

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
alpha = 1.9 # slope of inverse power law size frequency distribution
scaling = 1 # 1 = Sulak 2017 2 = Barker 2004
maxBergDepth = iceBergDepth # (m) - set to zero if 'prescribing' max iceberg width, set at top here
minBergDepth= 40 # (m)
maxBergWidth = 0 # (m) - set to zero if 'prescribing' max iceberg depth
minBergWidth = 40 # (m)

iceStart = 1
iceExtentIndex = int(np.round((iceExtent)/run_config['horiz_res_m']))

# Iceberg mask
bergMask[1:-1,iceStart:iceExtentIndex] = 1 # icebergs in inner 5 km, all oriented east-west

# Drift mask, No drift for Melange experiments, but can toggle on here if you want
# driftMask[1:-1,1:iceExtentIndex] = 1 # calculate effect of iceberg drift on melt rates 

# Melt mask, only let bergs melt in this region (make melt water, these don't change size)
meltMask[1:-1,iceStart:iceExtentIndex] = doMelt # Allow focus on blocking effect only

# Barrier mask
barrierMask[1:-1,iceStart:iceExtentIndex] = doBlock # make icebergs a physical barrier to water flow
barrierMask[plume_loc,icefront] = 0 #Plume code struggles with hFac adjustments

# Iceberg concentration (# of each surface cell that is filled in plan view)
bergConc[1:-1,iceStart:iceExtentIndex] = np.linspace(iceCoverage,1,(iceExtentIndex-1)) # iceberg concentration set at top
# bergConc[1:-1,iceStart:iceExtentIndex] = iceCoverage # iceberg concentration set at top

# print(bergConc[1:-1,1:iceExtentIndex])

desiredBergArea = np.sum(bergConc/100.0*deltaX*deltaY)
bergMaskArea = np.sum(bergMask*deltaX*deltaY)
setUpPrint('Area where bergs live: ' + str(bergMaskArea) + ' m^2')
setUpPrint('Desired berg area: ' + str(desiredBergArea) + ' m^2')
setUpPrint('Ratio: ' + str(desiredBergArea/bergMaskArea*100) + '%')

if(scaling == 1): # then use Sulak17 volume-area scaling volume = 6.0*area^1.30
    # assumes volume = L*W*D and W = L/1.62 (Dowdeswell et al 1992)
    if(maxBergWidth==0):
        maxBergWidth = 0.0642449*maxBergDepth**(5/3)
        # minBergWidth = 0.0642449*minBergDepth**(5/3)
    elif(maxBergDepth==0):
        maxBergDepth = 5.19155*maxBergWidth**(5/3)
        minBergDepth = 5.19155*minBergWidth**(5/3)
elif(scaling == 2): # Then use Barker04 width-depth relationship
    # Depth = 2.91*Width^0.71
    if(maxBergWidth==0):
        maxBergWidth = (100*10**(58/71)*maxBergDepth**(100/71)) / (291*291**(29/71))
        #minBergWidth = (100*10**(58/71)*minBergDepth**(100/71)) / (291*291**(29/71))        
    elif(maxBergDepth==0):
        maxBergDepth = 2.91*maxBergWidth^0.71
        minBergDepth = 2.91*minBergWidth^0.71

numberOfBergs = 50 #low start, immediately doubled by scheme below, so guess low, high guesses (300%+) can cause to fail
bergTopArea = 0
areaResidual = 1
# Generate the Inverse Power Law cumulative distribution function
# over the range minBergWidth-maxBergWidth with a slope of alpha.
setUpPrint('Making bergs, this can take a few loops...')
loop_count = 1

np.random.seed(2)
setUpPrint('random seed set, not really random anymore')

while(np.abs(areaResidual) > .005 ): # Create random power dist of bergs, ensure correct surface area
    numberOfBergs = round(numberOfBergs * (1 + areaResidual))  
    setUpPrint('\tnumberOfBergs: ' + str(numberOfBergs))
    x_width = np.arange(minBergWidth, maxBergWidth, (maxBergWidth-minBergWidth)/(numberOfBergs*1e2))
    inversePowerLawPDF_width = ((alpha-1) / minBergWidth) * (x_width/minBergWidth) ** (-alpha)
        # Get the CDF numerically
    inversePowerLawCDF_width = np.cumsum(inversePowerLawPDF_width)
        # Normalize
    inversePowerLawCDF_width = inversePowerLawCDF_width / inversePowerLawCDF_width[-1]
        
        # Generate number_of_bergs uniformly distributed random numbers.
    uniformlyDistributedRandomNumbers = np.random.uniform(0,1,numberOfBergs)
    
    inversePowerLawDistNumbers_width = np.zeros(uniformlyDistributedRandomNumbers.size);

    nearestIndex_width = [0] * uniformlyDistributedRandomNumbers.size

    nearestIndex_width = find_closest_indices(uniformlyDistributedRandomNumbers,inversePowerLawCDF_width)


    inversePowerLawDistNumbers_width = x_width[nearestIndex_width];
    inversePowerLawDistNumbers_length = inversePowerLawDistNumbers_width/1.12 # Widths are bigger 

    randScale = np.random.normal(6,1.22,numberOfBergs)
    randPower = np.random.normal(0.3,0.016,numberOfBergs)
    inversePowerLawDistNumbers_depth = randScale * (inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length) ** randPower * (920/1025)
    inversePowerLawDistNumbers_depth[inversePowerLawDistNumbers_depth < 5.123] = 5.123 #doest like round numbers, cap low end of bergs

    tooWide = np.count_nonzero(inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold))  #disallow completely full cells
    tooLong = np.count_nonzero(inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold))
    inversePowerLawDistNumbers_width[inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(hfacThreshold) # Max width is grid cell (assumed square)
    inversePowerLawDistNumbers_length[inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(hfacThreshold) # Max length is grid cell (assumed square)
    if(tooLong + tooWide > 0):
        setUpPrint('\t\tBergs clipped: %i for width, %i for length' % (tooWide, tooLong))
    
    bergTopArea = sum(inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length)
    areaResidual = (desiredBergArea - bergTopArea)/desiredBergArea
    setUpPrint('\t\t%.2f %% Bergs' % (bergTopArea/bergMaskArea*100))
    setUpPrint('\t\tareaResidual %.2f %%' % (areaResidual * 100))
    loop_count += 1
setUpPrint('====== Success! Found our bergs =====')
setUpPrint('Width min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_width),np.mean(inversePowerLawDistNumbers_width),np.max(inversePowerLawDistNumbers_width)))
setUpPrint('Depth min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth)))
setUpPrint('Total Berg Area %f' % bergTopArea)
setUpPrint('Total Berg fract: %.2f %%' % (bergTopArea/bergMaskArea*100))

# Now we sort these berg into cell, randomly 
bergMaski = 0  #Bad name, but this is the count of cells that will recieve bergs
bergDict = {}

for j in range(ny):
    for i in range(nx):
        if(bergMask[j,i] == 1):
            # print('i,j, bergmask',i,j,bergMask[j,i])
            bergMaski = 1 + bergMaski #Needs to start at 1, as non-bergs will be 0
            bergMaskNums[j,i] = bergMaski #Assign Mask Nums, not random as we'll randomly place bergs in cells
            bergDict[bergMaski] = [j,i] #This lets us do 1-D loops for the whole grid
setUpPrint('%i cells with bergs' % bergMaski)

# Sort my bergs
sorted_indices = np.argsort(-inversePowerLawDistNumbers_depth) # Sort backwards to get descending from big to small bergs 
sorted_depth = inversePowerLawDistNumbers_depth[sorted_indices]
sorted_width = inversePowerLawDistNumbers_width[sorted_indices]
sorted_length = inversePowerLawDistNumbers_length[sorted_indices]
assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) # In this script, every berg has a home

# Array for bergs
bergsPerCellLimit = 500
icebergs_depths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_widths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_length = np.zeros([bergMaski,bergsPerCellLimit])  #careful, not plural as to length match

np.random.seed(2)
assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) #every Berg has a spot

icebergs_per_cell = np.zeros([bergMaski],dtype=np.int16)
icebergs_area_per_cell = np.zeros([bergMaski])

for i in range(numberOfBergs): 
    j = assignedCell[i]
    # print('looking at mask number',j,'at berg',i)
    # print('Berg number', icebergs_per_cell[j],'in this cell')
    bergArea = sorted_width[i] * sorted_length[i]
    loopLimiter = 0
    while(bergArea > (deltaX * deltaY  * (bergConc[bergDict[j+1][0],bergDict[j+1][1]]/100) - icebergs_area_per_cell[j])): #if above 'full', pick random new cell, accept with decreasing probability
        j_old = j
        j = np.random.randint(0,bergMaski)
        loopLimiter += 1
        if((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY) < hfacThreshold - .01): #only consider accepting if under 95
            odds = np.abs(np.random.normal(0,.5,1))  #randomly accepts those that are big in overfull cells, but at decreasing frequency
            overFull = ((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY)*100 - bergConc[bergDict[j+1][0],bergDict[j+1][1]])
            if(odds > overFull):
                 # print('accepting overfull')
                 assignedCell[i] = j  #if we it a shuffling critera
                 break
        if(loopLimiter > bergMaski*20): #eventually we have to force some in
            indexesAllowed = np.where((deltaX * deltaY * hfacThreshold - icebergs_area_per_cell)  > bergArea)
            randi = np.random.randint(0,len(indexesAllowed[0]))
            j = indexesAllowed[0][randi]
            assignedCell[i] = j  #if we it a shuffling critera, must line up for calculation below
            # setUpPrint('\t Randomly missed, will force into cell with room: %i' % j)
            if((np.min(icebergs_area_per_cell) + bergArea)/(deltaX * deltaY) > hfacThreshold):
                setUpPrint('WARNING cell very full: %.2f%%' %((np.min(icebergs_area_per_cell) + bergArea)*100/(deltaX * deltaY)))
            break
     
    icebergs_depths[j,icebergs_per_cell[j]] = sorted_depth[i]
    icebergs_widths[j,icebergs_per_cell[j]] = sorted_width[i]
    icebergs_length[j,icebergs_per_cell[j]] = sorted_length[i]
    icebergs_per_cell[j] += 1
    # icebergs_area_per_cell[j] = np.sum(icebergs_widths[j,:]*icebergs_length[j,:])
    icebergs_area_per_cell[j] += bergArea
setUpPrint('Bergs per cell and filled faction at surface for spot check')     
setUpPrint(icebergs_per_cell)
setUpPrint(np.round(icebergs_area_per_cell/(deltaX*deltaY),2))
setUpPrint('Max fill is: %.2f%%' % (np.nanmax(icebergs_area_per_cell/(deltaX*deltaY))*100))

# All bergs now sorted 
openFrac = np.zeros([nz,ny,nx])
SA = np.zeros([nz,ny,nx])
SA[:,:,:] = np.nan

#This loop knows where all bergs are already, different from searching for all bergs across entire grid
for i in range(bergMaski):
    bergCount = icebergs_per_cell[i]
    numBergsPerCell[bergDict[i+1][0],bergDict[i+1][1]] = bergCount
    if(bergCount > 0):
        lengths = icebergs_length[i,icebergs_length[i,:] > 0] #return only non-zeros
        widths = icebergs_widths[i,icebergs_widths[i,:] > 0] #return only non-zeros
        depths = icebergs_depths[i,icebergs_depths[i,:] > 0] #return only non-zeros
        for k in range(nz):
            cellVolume = deltaX*deltaY*dz[k]
            d_bot = sum_z[k] #bottom of depth bin
            d_top = sum_z[k] - dz[k]
            volume1 = dz[k] * lengths[depths > d_bot] * widths[depths > d_bot]
            SA1 = dz[k]*2*(lengths[depths > d_bot] + widths[depths > d_bot])
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

# Plots for reference on whats happening berg-wise
fig = plt.figure()
plt.subplot(2,2,1)
for i in range(bergMaski):
    plt.plot(openFrac[:,bergDict[i+1][0],bergDict[i+1][1]],-sum_z,alpha=.5,color='xkcd:gray',linewidth=.5)
plt.plot(np.mean(openFrac[:,bergMask==1],1),-sum_z,alpha=1,color='xkcd:black',linewidth=1,linestyle='--',label='Average Bergs')
plt.plot([0,1],[-maxBergDepth,-maxBergDepth],color = 'xkcd:red',linestyle=':', label='Target Max Depth')
plt.plot([1-np.max(bergConc)/100,1-np.max(bergConc)/100],[-domain_params['H'],0],color = 'xkcd:gray',linestyle=':',label='Target Max Berg Conc')
plt.xlabel('Open Fraction of Cells')
plt.ylabel('Depth [m]')
# plt.legend()  #not quite room so off for now

plt.subplot(2,2,2)
plt.hist(inversePowerLawDistNumbers_depth,bins = 50)
plt.ylabel('Count')
plt.xlabel('Depth [m]')

plt.subplot(2,2,3)
plt.hist(inversePowerLawDistNumbers_width,bins = 50)
plt.ylabel('Count')
plt.xlabel('Width [m]')

plt.subplot(2,2,4)
plt.hist(inversePowerLawDistNumbers_length,bins = 50)
plt.ylabel('Count')
plt.xlabel('Length [m]')
fig.tight_layout()
if(writeFiles):
    plt.savefig(run_config['run_dir']+'/input/bergStatistics.png', format='png', dpi=200)
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
if(writeFiles):
    plt.savefig(run_config['run_dir']+'/input/bergMap.png', format='png', dpi=200)
plt.show()

fig = plt.figure()
pc = plt.pcolormesh(meltMask,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Melt Mask')
plt.ylabel('Cell across fjord')
plt.xlabel("Cell along fjord")
fig.tight_layout()
if(writeFiles):
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

icebergs_depths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
icebergs_widths2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])
icebergs_length2D = np.zeros([bergsPerCellLimit,grid_params['Ny'],grid_params['Nx']])

for k in range(bergMaski):
    j = bergDict[k+1][0]
    i = bergDict[k+1][1]
    icebergs_depths2D[:,j,i] = icebergs_depths[k,:]
    icebergs_widths2D[:,j,i] = icebergs_widths[k,:]
    icebergs_length2D[:,j,i] = icebergs_length[k,:]

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

#turn on ICEBERG if off, turn ICEPLUME off
if(makeDirs):
    replaceAll(run_config['run_dir'] + '/input/data.pkg','ICEBERG=.FALSE.', 'ICEBERG=.TRUE.') 
    # replaceAll(run_config['run_dir'] + '/input/data.pkg','ICEPLUME=.TRUE.', 'ICEPLUME=.FALSE.') 


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
estTime = int(grid_params['Ny']) * int(grid_params['Nx']) * int(grid_params['Nr']) * int(grid_params['Nt']) *2e-8
setUpPrint('Estimated run time is %.2f hours for one CPU' % (estTime/60))
setUpPrint('Estimated run time is %.2f hours for %i CPUs\n' % (estTime/60/ncpus*1.2,ncpus))

comptime_hrs = estTime/60/ncpus*1.2 

if(makeDirs):
    if os.path.isfile(run_config['run_dir']+'/input/setupReport.txt'):   
        os.remove(run_config['run_dir']+'/input/setupReport.txt')
        setUpPrint('previous setupReport.txt deleted in '+ run_config['run_dir']+'/input/')
    shutil.move('setupReport.txt', run_config['run_dir']+'/input')
    shutil.copy('fjordExperiment.py', run_config['run_dir']+'/input/buildScript.py')
    replaceAll(run_config['run_dir']+'/input/buildScript.py','makeDirs = True', 'makeDirs = False') 
    rcf.createSBATCHfile_Sherlock(run_config, cluster_params, walltime_hrs=1.2*comptime_hrs, email=email, mem_GB=1)
    setupNotes.close()
    print('Done! Remember to build before you run the script, building on MPI time is very inefficient')
elif(writeFiles):
    print("Done! You shouldn't have to rebuild as we only changed run time options here")
else:
    print('Nothing was saved, I hope you liked the pretty plots at least')
