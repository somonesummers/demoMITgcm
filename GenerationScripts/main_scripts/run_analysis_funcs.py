import numpy as np
import os
import importlib
import matplotlib.pylab as plt
#import MITgcmutils # see https://mitgcm.readthedocs.io/en/latest/utilities/utilities.html
#from IPython.core.debugger import set_trace # for debugging
import pickle
import run_config_funcs as rcf

output_dir = '../output'
os.makedirs(output_dir, exist_ok=True)


def get_exp_params_py(exp_path, iter0=0, add_hFac=True, add_forcing=True, 
                      vname='THETA', periodic_forcing_tlen=365):
    
    """
    Function to import run parameters and forcing
    """
    
    #exp_path = os.path.join(exp_root_dir, exp_name)
    exp_name = exp_path.split('/')[-1]
    results_path = os.path.join(exp_path, 'results')
    input_path = os.path.join(exp_path, 'input')
    
    # open run config file
    with open(os.path.join(input_path, 'params.p'), 'rb') as f:
        run_config = pickle.load(f)
       
    # import parameters (Note: using params structure from previous code for convenience)
    params = {}
    params['Nx'] = run_config['grid_params']['Nx']
    params['Ny'] = run_config['grid_params']['Ny']
    params['Nr'] = run_config['grid_params']['Nr']

    params['delX'] = run_config['grid_params']['delX']
    params['delY'] = run_config['grid_params']['delY']
    params['delR'] = run_config['grid_params']['delR']


    params['Lx'] = run_config['domain_params']['Lx']
    params['Ly'] = run_config['domain_params']['Ly']
    params['H'] = run_config['domain_params']['H']

    params['xx'] = np.cumsum(np.append([0], params['delX'][:-1])) + params['delX']/2 - params['Lx']/2
    params['yy'] = np.cumsum(np.append([0], params['delY'][:-1])) + params['delY']/2 
    params['zz'] = -(np.cumsum(np.append([0], params['delR'][:-1])) + params['delR']/2) 

    params_IO = run_config['mitgcm_params']['data'][2]
    params['deltaT'] = params_IO['deltaT']
    
    params['dumpIters'] = rcf.getSavedIters(vname, results_path)
    params['dumpIters'] = params['dumpIters'][params['dumpIters'] >= iter0]
    #print(iters)

    params['startTime'] = params['dumpIters']*params_IO['deltaT']
    params['nDumps'] = len(params['dumpIters'])
    params['nTimeSteps'] = np.ceil((params_IO['endTime'] - params['startTime'])/params['deltaT'])
    

    params['rho0'] = 1000 # kg/m3
    params['Cp'] = 4e3  # J/kg/K

    
    # extract dump freq (varies based on variable)
    diags = run_config['mitgcm_params']['data_diagnostics'][0]
    
    #print(diags)

    # map field name to dump freq
    vble_dump_freq = {}
    for key in diags:

        if key.startswith('fileName'):
            fname = diags[key]

            if fname not in vble_dump_freq:
                fnum = int(key[len('fileName')+1: -1])
                #print(fnum)
                vble_dump_freq[fname] = abs(diags['frequency(%i)' %fnum])
    
    params['dumpFreq'] = vble_dump_freq[vname]
    
    # pass mitgcm params and grid_params
    params['mitgcm_params'] = run_config['mitgcm_params']
    params['grid_params'] = run_config['grid_params']
    
    if add_hFac:
        params['DRF'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'DRF'))
        params['hFacS'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacS'))
        params['hFacW'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacW'))
        params['hFacC'] = MITgcmutils.mds.rdmds(os.path.join(results_path, 'hFacC'))

    # tack on exp_name and directory paths
    params['exp_name'] = exp_name
    params['exp_path'] = exp_path
    params['input_path'] = input_path
    params['results_path'] = results_path
    
    
    # get forcing and bathy files names
    input_fnames = run_config['mitgcm_params']['data'][4]
    params['bathyFile'] = input_fnames['bathyFile']
    params['zonalWindFile'] = input_fnames['zonalWindFile']
    params['surfQfile'] = input_fnames['surfQfile']
   
    # read in forcing data
    forcing = {}
    
    bathy_path = os.path.join(input_path, params['bathyFile'])
    forcing['bathy'] = np.fromfile(bathy_path, dtype='>f8').reshape(params['Ny'], params['Nx'])
    
    params03 = params['mitgcm_params']['data'][2]
    
    if 'periodicExternalForcing' in params03 and params03['periodicExternalForcing']:
        periodic_forcing = True
        forcing_shape = (periodic_forcing_tlen, params['Ny'], params['Nx'])
    else:
        forcing_shape = (params['Ny'], params['Nx'])

    if add_forcing:
        wind_path = os.path.join(input_path, params['zonalWindFile'])
        forcing['zonalWind'] = np.fromfile(wind_path, dtype='>f8').reshape(forcing_shape)

        wind_path = os.path.join(input_path, params['surfQfile'])
        forcing['surfQ'] = np.fromfile(wind_path, dtype='>f8').reshape(forcing_shape)
        
        if 'EmPmRFile' in input_fnames:
            params['EmPmRFile'] = input_fnames['EmPmRFile']
            emp_path = os.path.join(input_path, params['EmPmRFile'])
            forcing['EmPmR'] = np.fromfile(emp_path, dtype='>f8').reshape(forcing_shape)

        
    else:
        forcing['zonalWind'] = np.zeros(forcing_shape)
        forcing['surfQ'] = np.zeros(forcing_shape)
        forcing['EmPmR'] = np.zeros(forcing_shape)

    
    return params, forcing
    
    

    
def computeDomainInt(exp_path, vname, zr=[], xr=[], yr=[], tstep=1, load_previous=True, 
                     suppress_err_msg=True):
    
    #home_dir = '/home/earlew/research/scripts/MITgcm_py/'
    #output_dir = os.path.join(home_dir, 'output/')
   
        
    
    exp_name = exp_path.split('/')[-1]
        
    #load control params
    params, forcing = get_exp_params_py(exp_path, vname=vname)
    ##-1
    dumpIters = params['dumpIters']
    
#     print(nDumps)
#     print(dumpIters)
    
    zz = params['zz']
    xx = params['xx']
    yy = params['yy']
    
    secsInYr = 36000*24*365
    
    
    def get_domain_subset(sr, ss, dim):
        
        if len(sr)==2:
            assert sr[-1]>=sr[0]
            si = np.logical_and(ss>=sr[0], ss<sr[-1])
            if dim=='z':
                s_str = "_%s=%s-%sm" %(dim, np.abs(sr).min(), np.abs(sr).max())
            else:
                s_str = "_%s=%s-%skm" %(dim, np.abs(sr).min()/1e3, np.abs(sr).max()/1e3)
        else:
            si = np.logical_and(ss>=ss.min(), ss<=ss.max())
            s_str = ''
            
        return si, s_str


    zi, z_str = get_domain_subset(zr, zz, 'z')
    yi, y_str = get_domain_subset(yr, yy, 'y')
    xi, x_str = get_domain_subset(xr, xx, 'x')
        
    s_str = z_str+y_str+x_str
    # Get domain volume
    DX = params['delX'][xi, np.newaxis, np.newaxis]
    DY = params['delY'][np.newaxis, yi, np.newaxis]
    DZ = params['delR'][np.newaxis, np.newaxis, zi]
    DV = DX*DZ*DY*params['hFacC'][zi, :, :][:, yi, :][:, :, xi].T
    vol = np.sum(DV)
    
    #print(vol)
    
    rho0 = 1025 # kg/m3


            
    # pre-allocate containers and counters
    tt = []
    iters_list = [] # logs dumped files actually found

    vble_avg = []
    pdata = {'tt': [], 'vble_avg':[], 'dumpIters': []}
    fname = '%s_%s_avg%s_tseries_tstep%s.p'%(exp_name, vname, s_str, tstep)

    # load previous output
    output_fpath = os.path.join(output_dir, fname)
    
    n0 = 0
    if load_previous:
        try:
            pdata = pickle.load(open(output_fpath, 'rb'))
            n0 = np.argmin(np.abs(pdata['dumpIters'][-1]-dumpIters))+1 # set n0 so loop starts from the last saved dumpIter
            assert n0<=nDumps
        except FileNotFoundError:
            pass

    if suppress_err_msg:
        print_err = False
    else:
        print_err = True
        
    cc = 0 # counter
    nDumps = int(params['nDumps'])
    
    for n in range(n0, nDumps):          
        if n==0 or (n%tstep)==0 or n==(nDumps-1):

            try:
                vble = rdmdsWrapper(os.path.join(params['exp_path'], 'results/%s' %vname), dumpIters[n])
                
                # Calculate domain-mean potential temperature
                vble_avg.append(np.sum(vble[zi, :, :][:, yi, :][:, :, xi].T*DV)/vol)
                      
                tt.append(dumpIters[n]*params['deltaT'])
                iters_list.append(dumpIters[n])
                    
                # Increment counter
                cc += 1

            except OSError as e:
                if print_err:
                    print(str(e))
                    print("skipping over missing files...")
                    print_err = False
                continue

    if cc==0:
        print("Warning: no new data loaded")
    else:
        print("%s updated to year %s" %(exp_name, np.round(tt[-1]/secsInYr)))

    # save output
    for ii in range(len(tt)):
        pdata['tt'].append(tt[ii])
        pdata['dumpIters'].append(iters_list[ii])

        pdata['vble_avg'].append(vble_avg[ii])

    pickle.dump(pdata, open(output_fpath, 'wb')) 
    
    return pdata


def get_units(vname):

    if 'THETA' in vname:
        units = '$^{circ}$C'
    elif vname in ['SALT']:
        units = 'PSU'
    elif vname in ['UVEL', 'VVEL', 'WVEL'] or 'VEL' in vname or vname.startswith('SPD'):
        units = 'm/s'
    elif vname in ['SSH']:
        units = 'm'
    elif vname in ['BT_psi']:
        units = 'Sv'
    
    elif vname in ['SALT', 'S', 'SALT_inst']:
        units = 'PSU'
    
    else:
        units = ''
    
    #debug_here()
    return units



def get_slice_mean(vble, rr, vname, params, dim='x'):
    
    """
    Function that computes the mean of a slice in the x, y, or z direction.
    This function supercedes get_layer_mean() and get_zonal_mean().
    """
    
    assert dim in ['x', 'y', 'z'], "invalid coordinates"

    dim2 = dim+dim
    ii = np.logical_and(params[dim2]>=rr[0], params[dim2]<=rr[-1])
    
    # using if statements here is a bit clumsy
    # having vble be an xarray with labelled dimensions/coordinates would make this a lot cleaner
    if dim=='x':
        vble_avg = vble[ii, :, :].mean(axis=0)
    elif dim=='y':
        vble_avg = vble[:, ii, :].mean(axis=1)
    elif dim=='z':
        vble_avg = vble[:, :, ii].mean(axis=2)
        
    return vble_avg





    
