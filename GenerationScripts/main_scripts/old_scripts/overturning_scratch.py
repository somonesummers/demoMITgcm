import numpy as np
import os
import matplotlib.pylab as plt
import helper_functions as hf
import imp
import xarray as xr
from numba import jit

from IPython.core.debugger import set_trace
debug_here = set_trace

home_dir = '/home/earlew/research/scripts/MITgcm_py/'
plot_dir = os.path.join(home_dir, 'plots_temporary/')
output_dir = os.path.join(home_dir, 'psi_output/')

imp.reload(hf)


def calcOverturning_v2(exp_name, tr, Nurser_Lee_remap=False):
    
    """
    Function to calculate overturning using the MITgcm layers package.
    
    EW attempt to make this function from scratch (sort of)
    """
    
    params,_ = hf.get_exp_params(exp_name)
    
    
    # Density bins for MOC calculation
    ptlevs = params['layers_bounds']
    Npt = len(params['layers_bounds'])-1
    
    hFacS = params['hFacS'].T
#     Nx = params['Nx']
#     Ny = params['Ny']
#     Nr = params['Nr']
#     delR = params['delR']
#     zz = -np.cumsum((params['delR'] + np.append([0], params['delR'][:params['Nr']-1]))/2) # same as params[zz]?
    
    ## Load flux output  
    print("Loading time-averaged layer flux output...")
    vh_iso = hf.readIters(vname='LaVH1TH', tr=tr, params=params, mask_zeros=False) # product of layer thickness and velocity
    
    # get residual overturning 
    print("Computing residual overturning in density space...")
    dx2d  = params['delX'][:, np.newaxis]
    Npt = len(params['layers_bounds'])-1
    psi_res_pt = np.zeros((params['Ny'], Npt))
    vh_iso_rev = vh_iso[:, :, ::-1] # to ensure we integrate from top to bottom?
    for ii in range(1, Npt):
        psi_res_pt[:, ii] = psi_res_pt[:, ii-1] + np.sum(vh_iso_rev[:, :, ii-1]*dx2d, axis=0)

    psi_res_pt= psi_res_pt[:, ::-1]
    debug_here()
    
    # map streamfunction to z-coordinates
    print("Mapping residual overturning to z-coordinates...")
    h_iso = hf.readIters(vname='LaHs1TH', tr=tr, params=params, mask_zeros=False) # thickness (sums to ocean depth)
    
    # get isopycnal depths
    if Nurser_Lee_remap:
        ###------Alternative remapping method (Nurser and Lee 2004)----------#
        print("Using Nurser and Lee (2004)")
        # calculate x-z area below each isopycnal layer at each latitude
        h_iso_xzArea = (h_iso*params['delX'][:, np.newaxis, np.newaxis]).sum(axis=0).cumsum(axis=1)

        # calculate x-z ocean area below at each depth layer at each latitude
        z_xzArea_temp = (hFacS*delR[np.newaxis, np.newaxis, :]*params['delX'][:, np.newaxis, np.newaxis]).sum(axis=0)
        z_xzArea = z_xzArea_temp[:, ::-1].cumsum(axis=1) #delR goes from top to bottom. Need to integrate from bottom to top

        # interpolate to get isopycnal depths based on x-z area bounded by each isopycnal
        from scipy.interpolate import interp1d
        z_pt = np.nan*np.zeros(h_iso_xzArea.shape)
        for jj in range(Ny):
            z_intp = interp1d(z_xzArea[jj, :], zz[::-1], bounds_error=False)
            z_pt[jj, :] = z_intp(h_iso_xzArea[jj, :].filled())
         
        remap_str = '_Nurser_Lee'
        print(z_pt[10, :5])
    else:
        
        z_pt = np.zeros(h_pt_xtavg.shape)
        for m in range(Npt):
            z_pt[:, m] = -np.sum(h_pt_xtavg[:,:m-1], 1)
            
        remap_str = ''
        print(z_pt[10, :5])
    
    
    
    # With the isopycnal depth, do the actual remapping
    psi_res_z = np.zeros((params['Ny'], params['Nr']))
    for jj in range(params['Ny']):
        intp = interp1d(z_pt[jj, :], psi_res_pt[jj, :], bounds_error=False)
        psi_res_z[jj, :] = intp(params['zz'])
        
        
    #print("Loading other time-averaged output...")
    
    #theta = hf.readIters(vname='THETA', tr=tr, params=params, mask_zeros=False)
    #vvel = hf.readIters(vname='VVEL', tr=tr, params=params, mask_zeros=False)
    
    ## shift potential temperature to v-gridpoints  
#     theta_vgrid = np.zeros(theta.shape)*np.nan
#     theta_vgrid[:, 1:params['Ny'], :] = 0.5*(theta[:, :params['Ny']-1, :] + theta[:, 1:params['Ny']])
    
    
#     vflux = vvel*h_iso
    
#     ## compute mean flux within each density layer
#     vvel_dz = vvel*hFacS*params['delR'][np.newaxis, np.newaxis, :] # multply vvel by wet cell height
    
#     vflux_mean = np.zeros(vflux.shape)
    
#     # compute volume fluxes that fall outside ptlevs range
#     ii_above_range = theta>ptlevs[Npt] 
#     vflux_mean[:, :, Npt-1] = vflux_mean[:, :, Npt-1] + np.sum(vdz*ii_above_range, 2)
#     ii_below_range = pt_f<=ptlevs[1]
#     vflux_mean[:, :, 0] = vflux_mean[:, :, 0] + np.sum(vdz*ii_below_range, 2) 
    
#     # compute layer volume flux
#     print("computing vflux_mean (layer volume flux)...")
#     for ii in range(1, Npt-1):
#         tmp = np.logical_and(theta>ptlevs[m], theta<=ptlevs[m+1]) # find theta layer
#         vflux_mean[:, :, ii] = vflux_mean[:, :, ii] + np.sum(vdz*tmp, 2)
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    