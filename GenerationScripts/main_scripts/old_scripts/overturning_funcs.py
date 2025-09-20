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

@jit(nopython=True)
def get_interp_grid(ffac, hFacS, hFacS_f, Nx, Ny, Nr, Nrf, delRf, delR):
    
    ZZ = np.zeros((Nx, Ny, Nr))
    ZZ_f = np.zeros((Nx, Ny, Nrf))
    DZ = np.zeros((Nx, Ny, Nr))
    DZ_f = np.zeros((Nx, Ny, Nrf))
    PP = np.zeros((Nx, Ny, Nrf))
    
    ZZ[:, :, 0] = -delR[0]*hFacS[:,:,0]/2
    for k in range(1, Nr):
        ZZ[:,:,k] = ZZ[:,:,k-1] - 0.5*delR[k-1]*hFacS[:,:,k-1] - 0.5*delR[k]*hFacS[:,:,k]
        
    ZZ_f[:, :, 0] = -delRf[0]*hFacS[:,:,0]/2
    for k in range(1, Nrf):
        ZZ_f[:,:,k] = ZZ_f[:,:,k-1] - 0.5*delRf[k-1]*hFacS_f[:,:,k-1] - 0.5*delRf[k]*hFacS_f[:,:,k]
                       
    for k in range(Nr):
        DZ[:, :, k] = delR[k]
        
    for k in range(Nrf):
        DZ_f[:, :, k] = delRf[k]
    
#     for k in range(Nr):
#         PP[:, :, k] = -delR[k] # NOTE: Looks like PP = -DZ
    
#     dZZ = ZZ_f.mean()-ZZ.mean()
#     assert np.abs(dZZ)<1e-10, "means are not consistent"
    
    #debug_here()
        
    
    ### Create matrices for vertical interpolation
    #print("creating matrices for vertical integration...")
    k_p = np.zeros((Nx, Ny, Nrf)) # k: vertical index. k_p = k+1, k_n = k-1
    k_n = np.zeros((Nx, Ny, Nrf))
    w_n = np.zeros((Nx, Ny, Nrf)) # w: weights?
    w_p = np.zeros((Nx, Ny, Nrf))
    is_wet_col = np.zeros((Nx, Ny))
    
    # Need to be careful here since we are dealing with different index conventions.
    for ii in range(Nx):
        for jj in range(Ny):
            # indices of the lowest cells
            hFacS_ij = hFacS[ii, jj, :]
            kmax = int(np.sum(hFacS_ij!=0))
            #kmax_f = ffac*kmax
            is_wet_col[ii, jj] = kmax!=0
            
#             if is_wet_col[ii, jj] == 0:
#                 continue
#             else:
#                 kmax = kmax - 1
            
#             assert kmax<=Nr, "kmax too large"
#             assert kmax>=0, "kmax too small"
            
            for k in range(Nrf):
                kk = k + 1
                # Previous and next interpolation indices 
                k_p[ii, jj, k] = np.ceil((kk)/ffac-0.5) # ceil() makes things a bit tricky. Can't just offset by 1.
                k_n[ii, jj, k] = k_p[ii, jj, k] + 1
                
                # print(k_p[ii, jj, k])
                # Fine grid cell is above highest coarse grid cell, so fine grid
                # gamma will just be set equal to uppermost coarse grid gamma (EW: not sure what gamma is referring to here)
                if k_p[ii, jj, k] <= 0:
                    
                    k_p[ii, jj, k] = 1
                    w_p[ii, jj, k] = 0
                    w_n[ii, jj, k] = 1
                    
                else:
                    # Fine grid cell is below lowest coarse grid cell, so fine grid 
                    # gamma will just be set equal to lowermost coarse grid gamma
                    if k_n[ii, jj, k] > kmax:
#                         while k_n[ii, jj, k] > kmax:
#                             k_n[ii, jj, k] = k_n[ii, jj, k] -1

#                         k_p[ii, jj, k] = k_n[ii, jj, k]-1
                        
                        k_n[ii, jj, k] = kmax
                        w_n[ii, jj, k] = 0
                        w_p[ii, jj, k] = 1
                       
                            
                    else:
                        # Otherwise set weights to interpolate linearly between neighboring
                        # coarse-grid gammas
                        k_n_ijk = int(k_n[ii, jj, k])-1
                        k_p_ijk = int(k_p[ii, jj, k])-1
                        assert (k_n_ijk-k_p_ijk)==1, "value error"
                        
                        w_p[ii, jj, k] = (ZZ[ii, jj, k_n_ijk]-ZZ_f[ii, jj, k])/(ZZ[ii, jj, k_n_ijk]-ZZ[ii, jj, k_p_ijk])
                        
                        err = 1e-10
                        assert w_p[ii, jj, k]>=0-err and w_p[ii, jj, k]<=1+err, "w_p has incorrect size"
                        
                        if w_p[ii, jj, k]>=0-err and w_p[ii, jj, k]<0:
                            w_p[ii, jj, k] = 0
                        
                        w_n[ii, jj, k] = 1 - w_p[ii, jj, k]
                        
    return w_p, w_n, k_p, k_n, is_wet_col, DZ_f 

@jit(nopython=True, parallel=True)
def calc_psi(Npt, pt_f, ptlevs, vdz, vflux_mean):
    for m in range(1, Npt-1):
        tmp = np.logical_and(pt_f>ptlevs[m], pt_f<=ptlevs[m+1])
        vflux_mean[:, :, m] = vflux_mean[:, :, m] + np.sum(vdz*tmp, 2)
        
    return vflux_mean

@jit(nopython=True, parallel=True)
def xint_psi(Nx, delX, vflux_tavg, vflux_mean, vflux_xint, vflux_mean_xint):
    for ii in range(Nx):
        vflux_xint = vflux_xint + delX[ii]*vflux_tavg[ii, :, :]
        vflux_mean_xint = vflux_mean_xint + delX[ii]*vflux_mean[ii, :, :]
        
    return vflux_xint, vflux_mean_xint


def calcOverturning(exp_name, tr, Nurser_Lee_remap=False):
    """
    Function to calculate overturning using the MITgcm layers package.
    This is python translation of calcOverturning.m written by Andrew Stewart.
    
    tr: array defining the time range to do calculation
    """
    
    params,_ = hf.get_exp_params(exp_name)
    
    # Density bins for MOC calculation
    ptlevs = params['layers_bounds']
    Npt = len(ptlevs)-1
    
    # Frequency of diagnostic output
    dumpFreq = np.abs(params['diag_frequency']['THETA_inst'])
    nDumps = np.round(params['nTimeSteps']*params['deltaT']/dumpFreq)
    dumpIters = np.round(np.arange(1, nDumps+1)*dumpFreq/params['deltaT'])
    dumpIters = dumpIters[dumpIters >= params['nIter0']]
    nDumps = len(dumpIters)
    
    # Create a finer vertical grid (NOTE: might be able to do this using interp)
    ffac = 5  # amount by which to increase vertical resolution
    Nrf = ffac*params['Nr']
    delRf = np.zeros(Nrf)
    for n in range(params['Nr']):
        for m in range(ffac):
            delRf[(n)*ffac+m] = params['delR'][n]/ffac
            
    zz = -np.cumsum((params['delR'] + np.append([0], params['delR'][:params['Nr']-1]))/2)
    zz_f = -np.cumsum((delRf + np.append([0], delRf[:Nrf-1]))/2)
    
    ### create finer hFacS
    hFacS_f = np.zeros((params['Nx'], params['Ny'], Nrf))
    hFacS = params['hFacS'].T
    for k in range(params['Nr']):
        hFacS_k = hFacS[:,:,k]
        hFacS_f[:, :, ffac*(k): ffac*(k+1)] =  hFacS_k[:, :, np.newaxis]
     
    dhFacS = hFacS_f.mean()-hFacS.mean()
    assert np.abs(dhFacS)<1e-10, "hFacS_f and hFacS means are not consistent"
    
    print("creating grid for vertical positions...")
    
    #---------------------------#
    Nx = params['Nx']
    Ny = params['Ny']
    Nr = params['Nr']
    delR = params['delR']
    w_p, w_n, k_p, k_n,is_wet_col, DZ_f  = get_interp_grid(ffac, hFacS, hFacS_f, Nx, Ny, Nr, Nrf, delRf, delR)
    
    
    # adjust indices to comply with python convention
    k_p = k_p - 1
    k_n = k_n - 1
    k_p[k_p<0] = 0
    k_n[k_n<0] = 0
    
    #debug_here()
    err = 1e-10
    assert np.min(k_n) >=0 and np.max(k_n) <params['Nr'], "Value error. min(k_n)=%s, max(k_n)=%s" %(np.min(k_n), np.max(k_n))
    assert np.min(k_p) >=0 and np.max(k_p) <params['Nr'], "Value error. min(k_p)=%s, max(k_p)=%s" %(np.min(k_p), np.max(k_p))
    assert np.min(w_p) >=0 and np.max(w_p) <=1+err, "Value error. min(w_p)=%s, max(w_p)=%s" %(np.min(w_p), np.max(w_p))
    #---------------------------#
    
    
    ### Load flux output  
    print("Loading flux output...")
    vflux_tavg = hf.readIters(vname='LaVH1TH', tr=tr, params=params, mask_zeros=False)
    h_pt_tavg = hf.readIters(vname='LaHs1TH', tr=tr, params=params, mask_zeros=False)
    pt_tavg = hf.readIters(vname='THETA', tr=tr, params=params, mask_zeros=False)
    vvel_tavg = hf.readIters(vname='VVEL', tr=tr, params=params, mask_zeros=False)
    
    debug_here()
    
    ### Interpolate potential temperature to v-gridpoints  
    pt_v = np.zeros(pt_tavg.shape)*np.nan
    pt_v[:, 1:params['Ny'], :] = 0.5*(pt_tavg[:, :params['Ny']-1, :] + pt_tavg[:, 1:params['Ny']])
    
    ### Interpolate onto a finer grid 
    vvel_f = np.zeros((params['Nx'], params['Ny'], Nrf))
    pt_f = np.nan*np.zeros((params['Nx'], params['Ny'], Nrf)) 
    
    print("Interpolating temp and vvel to finer grid...")
    if ffac == 1:
        # Shortcut if fine grid resolution = coarse grid resolution
        vvel_f = vvel_tavg     
        pt_f = pt_v
        
    else:
        # Velocity uniform throughout each coarse grid cell to preserve mass conservation
        for k in range(params['Nr']):
            vvel_tavg_k = vvel_tavg[:,:,k]
            vvel_f[:, :, ffac*k:ffac*(k+1)] = vvel_tavg_k[:, :, np.newaxis]
         
        # print(vvel_f.mean())
        # debug_here()
        # Linearly interpolate density
        for ii in range(params['Nx']):
            for jj in range(params['Ny']):
                if is_wet_col[ii, jj]:
                    k_p_ij = k_p[ii, jj, :].astype(int)
                    k_n_ij = k_n[ii, jj, :].astype(int)
                    pt_f[ii, jj,:] = w_p[ii, jj,:]*pt_v[ii, jj, k_p_ij] + w_n[ii, jj, :]*pt_v[ii, jj, k_n_ij]

    
    ### Calculate mean fluxes within mean density surfaces
    print("Computing streamfunction...")
    vflux_mean = np.zeros(vflux_tavg.shape)
    vdz = vvel_f*hFacS_f*DZ_f # Eulerian overturning?
    
    pt_f_i1 = pt_f>ptlevs[Npt] # note: this is a boolean array but numpy flips values to int during calculations
    vflux_mean[:, :, Npt-1] = vflux_mean[:, :, Npt-1] + np.sum(vdz*pt_f_i1, 2)
    pt_f_i2 = pt_f<=ptlevs[1]
    vflux_mean[:, :, 0] = vflux_mean[:, :, 0] + np.sum(vdz*pt_f_i2, 2) 
    
    #debug_here()
    
    #---------------------------#
    vflux_mean = calc_psi(Npt, pt_f, ptlevs, vdz, vflux_mean)
    #---------------------------#
#     for m in range(1, Npt-1):
#         tmp = np.logical_and(pt_f>ptlevs[m], pt_f<=ptlevs[m+1])
#         vflux_mean[:, :, m] = vflux_mean[:, :, m] + np.sum(vdz*tmp, 2)
        
    ### Zonally integrate meridional fluxes
    vflux_xint = np.zeros((params['Ny'], Npt))
    vflux_mean_xint = np.zeros((params['Ny'], Npt))
    #---------------------------#
    vflux_xint, vflux_mean_xint = xint_psi(params['Nx'], params['delX'], vflux_tavg, vflux_mean, vflux_xint, vflux_mean_xint)
    #---------------------------#
    
#     for ii in range(params['Nx']):
#         vflux_xint = vflux_xint + params['delX'][ii]*vflux_tavg[ii, :, :]
#         vflux_mean_xint = vflux_mean_xint + params['delX'][ii]*vflux_mean[ii, :, :]
                              
    ### Sum fluxes to obtain streamfunction (go from m2/s to m3/s)
    psi_pt = np.zeros((params['Ny'], Npt+1))
    psi_mean_pt = np.zeros((params['Ny'], Npt+1))
    for m in range(Npt):
        psi_pt[:, m] = np.sum(vflux_xint[:, m:Npt], 1)
        psi_mean_pt[:, m] = np.sum(vflux_mean_xint[:, m:Npt], 1)
        
    
    psi_pt = psi_pt/1e6 # units now in Sv
    psi_mean_pt = psi_mean_pt/1e6
    psi_eddy_pt = psi_pt - psi_mean_pt
    
    ### Calculate mean density surface heights
    h_pt_xtavg = np.nanmean(h_pt_tavg, axis=0)
    
    
    if Nurser_Lee_remap:
        ###------Alternative remapping method (Nurser and Lee 2004)----------#
        print("Mapping to z-coordinates using Nurser and Lee (2004)")
        # calculate x-z area below each isopycnal layer at each latitude
        h_pt_tavg_xzArea = (h_pt_tavg*params['delX'][:, np.newaxis, np.newaxis]).sum(axis=0).cumsum(axis=1)

        # calculate x-z ocean area below at each depth layer at each latitude
        z_xzArea_temp = (hFacS*delR[np.newaxis, np.newaxis, :]*params['delX'][:, np.newaxis, np.newaxis]).sum(axis=0)
        z_xzArea = z_xzArea_temp[:, ::-1].cumsum(axis=1) #delR goes from top to bottom. Need to integrate from bottom to top

        # interpolate to get isopycnal depths based on x-z area bounded by each isopycnal
        from scipy.interpolate import interp1d
        z_pt = np.nan*np.zeros(h_pt_xtavg.shape)
        for jj in range(Ny):
            z_intp = interp1d(z_xzArea[jj, :], zz[::-1], bounds_error=False)
            z_pt[jj, :] = z_intp(h_pt_tavg_xzArea[jj, :].filled())
         
        remap_str = '_Nurser_Lee'
        print(z_pt[10, :5])
    else:
        
        z_pt = np.zeros(h_pt_xtavg.shape)
        for m in range(Npt):
            z_pt[:, m] = -np.sum(h_pt_xtavg[:,:m-1], 1)
            
        remap_str = ''
        print(z_pt[10, :5])
        

    ### Calculate zonal-mean potential temperature
    pt_xtavg = np.nanmean(pt_tavg, axis=0)
    pt_f_xtavg = np.nanmean(pt_f, axis=0)
    
    
    ### Convert to z-coordinates by mapping the streamfunction at each temp
    ### level to the mean height of that density surface
    print("mapping streamfunction back to z-coords...")
    psi_z = np.nan*np.ones((params['Ny'], Nrf))
    psi_mean_z = np.nan*np.ones((params['Ny'], Nrf))
    psi_eddy_z = np.nan*np.ones((params['Ny'], Nrf))
    
    debug_here()
    # EW remake
    for jj in range(params['Ny']):
        
        
        pt_intp = interp1d(z_pt[jj, :], psi_pt[jj, :], bounds_error=False)
        #psi_z[jj, :] = 
    
    
    for jj in range(params['Ny']):
        for kk in range(Nrf):
            
            # Density lies in the lowest bin
            # EW: to handle cases where density/temp is less than lowest bin
            if pt_f_xtavg[jj, kk] < ptlevs[1]:
                psi_z[jj, kk] = psi_pt[jj, 0]      
                psi_mean_z[jj, kk] = psi_mean_pt[jj, 0]   
                psi_eddy_z[jj, kk] = psi_eddy_pt[jj, 0] 
                continue
            
            # Density lies in the highest bin
            # EW: to handle cases where density/temp is greater than highest bin
            if pt_f_xtavg[jj, kk] > ptlevs[Npt]:
                psi_z[jj, kk] = psi_pt[jj, Npt]      
                psi_mean_z[jj, kk] = psi_mean_pt[jj, Npt]   
                psi_eddy_z[jj, kk] = psi_eddy_pt[jj, Npt] 
                continue 
                
                
            # Density lies in an intermediate bin, so find the bin and assign
            # the overturning streamfunction via linear interpolation
            for m in range(1, Npt-1):
                if pt_f_xtavg[jj, kk] < ptlevs[m+1]:
                    pt_n = ptlevs[m+1]
                    pt_p = ptlevs[m]
                    wp = (pt_n-pt_f_xtavg[jj, kk])/(pt_n-pt_p)
                    wn = 1 - wp
                    psi_z[jj, kk] = wp*psi_pt[jj, m] + wn*psi_pt[jj, m+1];
                    psi_mean_z[jj, kk] = wp*psi_mean_pt[jj, m] + wn*psi_mean_pt[jj, m+1];
                    psi_eddy_z[jj, kk] = wp*psi_eddy_pt[jj, m] + wn*psi_eddy_pt[jj, m+1];
                    break

    
    
    #debug_here()
    # Finally, save output as xarray dataset
    """
      save(fullfile(outdir,[expname,'_MOC_pt.mat']), ...
    'xx','yy','zz','zz_f','hFacS_f','delRf','ptlevs', ... 
    'vvel_tavg','vvel_f','pt_tavg','pt_f', ...  
    'vflux_tavg','vflux_m', ...
    'psi_pt','psim_pt','psie_pt', ...
    'psi_z','psim_z','psie_z');
    
    """
    
    ptlevs_mean = ptlevs[:-1] + (ptlevs[1:]-ptlevs[:-1])/2
    
    psi_ds = xr.Dataset({'hFacS_f': (['xx', 'yy', 'zz_f'], hFacS_f), 'delRf': (['zz_f'], delRf),
                         'psi_pt': (['yy', 'r_pt'], psi_pt), 'psi_mean_pt': (['yy', 'r_pt'], psi_mean_pt), 
                         'psi_eddy_pt': (['yy', 'r_pt'], psi_eddy_pt), 'psi_z': (['yy', 'zz_f'], psi_z),
                         'psi_mean_z': (['yy', 'zz_f'], psi_mean_z), 'psi_eddy_z': (['yy', 'zz_f'], psi_eddy_z), 
                         'z_pt': (['yy', 'ptlevs_mean'], z_pt),
                         'xx': params['xx'], 'yy': params['yy'], 'zz': params['zz'], 
                         'zz_f': zz_f, 'r_pt': ptlevs, 'r_pt_mean': ptlevs_mean, })
    
    #'vvel_tavg': (['xx', 'yy', 'zz'], vvel_tavg), 'vvel_f': (['xx', 'yy', 'zz_f'], vvel_f),
    #'vflux_tavg': (['xx', 'yy', 'r_pt_mean'], vflux_tavg), 
    #'pt_tavg': (['xx', 'yy', 'zz'], pt_tavg), 'pt_f': (['xx', 'yy', 'zz_f'], pt_f),
    # 'vflux_mean': (['xx', 'yy', 'r_pt_mean'], vflux_mean), 
    
    fname = params['exp_name'] + remap_str + '_yrs_%s-%s.nc' %(tr[0], tr[-1])
    psi_ds.to_netcdf(os.path.join(output_dir, fname))
    #debug_here()
    
    

    
    
    
    
    
    
    
    
    
    
    
            
    ### create grid of actual vertical positions accounting for partial cells
#     print("creating grid for vertical positions...")
#     ZZ = np.zeros((params['Nx'], params['Ny'], params['Nr']))
#     ZZ_f = np.zeros((params['Nx'], params['Ny'], Nrf))
#     DZ = np.zeros((params['Nx'], params['Ny'], params['Nr']))
#     DZ_f = np.zeros((params['Nx'], params['Ny'], Nrf))
#     PP = np.zeros((params['Nx'], params['Ny'], Nrf))
    
#     ZZ[:, :, 0] = -params['delR'][0]*hFacS[:,:,0]/2
#     for k in range(1, params['Nr']):
#         ZZ[:,:,k] = ZZ[:,:,k-1] - 0.5*params['delR'][k-1]*hFacS[:,:,k-1] - 0.5*params['delR'][k]*hFacS[:,:,k]
        
#     ZZ_f[:, :, 0] = -delRf[0]*hFacS[:,:,0]/2
#     for k in range(1, Nrf):
#         ZZ_f[:,:,k] = ZZ_f[:,:,k-1] - 0.5*delRf[k-1]*hFacS_f[:,:,k-1] - 0.5*delRf[k]*hFacS_f[:,:,k]
                       
#     for k in range(params['Nr']):
#         DZ[:, :, k] = params['delR'][k]
        
#     for k in range(Nrf):
#         DZ_f[:, :, k] = delRf[k]
    
#     for k in range(params['Nr']):
#         PP[:, :, k] = -params['delR'][k] # NOTE: Looks like PP = -DZ
        
    
#     ### Create matrices for vertical interpolation
#     print("creating matrices for vertical integration...")
#     k_p = np.zeros((params['Nx'], params['Ny'], Nrf)) # k: vertical index. k_p = k+1, k_n = k-1
#     k_n = np.zeros((params['Nx'], params['Ny'], Nrf))
#     w_n = np.zeros((params['Nx'], params['Ny'], Nrf)) # w: weights?
#     w_p = np.zeros((params['Nx'], params['Ny'], Nrf))
#     is_wet_col = np.zeros((params['Nx'], params['Ny']))
    
#     # EW: not quite sure what's happening below
#     # Need to be careful here since we are dealing with different index conventions.
#     for ii in range(params['Nx']):
#         for jj in range(params['Ny']):
#             # indices of the lowest cells
#             hFacS_ij = hFacS[ii, jj, :]
#             kmax = int(np.sum(hFacS_ij[hFacS_ij!=0]))-1 # need to offset index to comply with Python indexing
#             kmax_f = ffac*kmax
#             is_wet_col[ii, jj] = kmax!=0
            
#             #print(kmax)
#             if kmax>=(ZZ.shape[2]):
#                 #print(kmax)
#                 print("AHHHHHHHHHH!!!")
#                 debug_here()
#                 #pass
            
#             for k in range(Nrf):
#                 # Previous and next interpolation indices (TODO: check that these slight mods are consistent with OG code)
#                 k_p[ii, jj, k] = np.ceil(k/ffac-0.5)
#                 k_n[ii, jj, k] = k_p[ii, jj, k] + 1
                
#                 # print(k_p[ii, jj, k])
#                 # Fine grid cell is above highest coarse grid cell, so fine grid
#                 # gamma will just be set equal to uppermost coarse grid gamma (EW: not sure what gamma is referring to here)
#                 if k_p[ii, jj, k] <= 0:
#                     k_p[ii, jj, k] = 1
#                     w_p[ii, jj, k] = 0
#                     w_n[ii, jj, k] = 1
                    
#                 else:
#                     # Fine grid cell is below lowest coarse grid cell, so fine grid 
#                     # gamma will just be set equal to lowermost coarse grid gamma
#                     if  k_n[ii, jj, k] > kmax:
#                         k_p[ii, jj, k] = int(kmax)
#                         w_p[ii, jj, k] = 0
#                         w_n[ii, jj, k] = 1
                        
#                     else:
#                         # Otherwise set weights to interpolate linearly between neighboring
#                         # coarse-grid gammas
#                         k_n_ijk = k_n[ii, jj, k].astype(int)
#                         k_p_ijk = k_p[ii, jj, k].astype(int)
#                         #debug_here()
#                         w_p[ii, jj, k] = (ZZ[ii, jj, k_n_ijk]-ZZ_f[ii, jj, k])/(ZZ[ii, jj, k_n_ijk]-ZZ[ii, jj, k_p_ijk])
#                         w_n[ii, jj, k] = 1 - w_p[ii, jj, k]
                        
#                         #debug_here()
                            
    
    
    
    
    
    
    

