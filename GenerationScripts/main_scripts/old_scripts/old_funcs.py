import numpy as np
import os
import matplotlib.pylab as plt
import imp
import xarray as xr
import matplotlib.gridspec as gridspec
import seawater as sw
import pickle
#import matplotlib.animation as animation
from scipy import stats, interpolate

import helper_functions as hf
import plot_functions as pf

from IPython.core.debugger import set_trace
debug_here = set_trace

imp.reload(hf)
imp.reload(pf)

plot_dir = 'plots_temporary/'
output_dir = 'output/'


def plot_bathy_2D(exp_names, exp_names_alias=[], fz=14, save_plots=True):
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    axes = axes.flatten()
    for ii, exp_name in enumerate(exp_names):
        
        plt.sca(axes[ii])
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii]  
        else:
            exp_name_str = exp_names_alias[ii]
        
        params, forcing = hf.get_exp_params(exp_name)

        xx = params['xx']/1e3
        yy = params['yy']/1e3
        ZZ = forcing['bathy']
        ZZ[:, -2:] = -4000
        
        XX, YY = np.meshgrid(xx, yy)
            
        plt.contourf(XX, YY, ZZ, 30, cmap=plt.cm.inferno, vmin=-4000, vmax=0)
        plt.colorbar()
        plt.xlabel('X (km)', fontsize=fz)
        plt.ylabel('Y (km)', fontsize=fz)
        plt.title(exp_name_str, fontsize=fz)
        
    
    plt.subplots_adjust(hspace=0.4)
    
    if save_plots:
        fname = 'bathy_comp_2D.pdf'
        plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')
        
        
def plot_bathy_tau_surfq(exp_name, exp_name_alias, fz=14, xlim1=[-10, 10], xlim2=[-0.15, 0.15], save_plots=True):
    
    # load surface forcing and bathy
    params, forcing = hf.get_exp_params(exp_name)
    yy = params['yy']/1e3
    xx = params['xx']/1e3
    tau_xmean = forcing['zonalWind'].mean(axis=1)
    q_xmean = forcing['surfQ'].mean(axis=1)
    ZZ = forcing['bathy']
    ZZ[:, -2:] = -4000
    
    fig = plt.figure(figsize=(15, 6))
    gs = gridspec.GridSpec(1, 3)
    #gs = fig.add_gridspec(1, 3)
    
    
    # plot surface forcing
    ax1 = plt.subplot(gs[0, 0])
    #ax1.set_title('Surface Forcing', fontsize=fz)
    
    col1 = 'red'
    ax1 = plt.gca()
    ax1.plot(q_xmean, yy, linewidth=2, color=col1)
    ax1.set_xlabel("Ocean heat loss (W/m$^2$)", fontsize=fz, color=col1)
    ax1.tick_params(axis='x', labelcolor=col1)
    ax1.set_xlim(xlim1)
    ax1.grid(True)
    ax1.tick_params(axis='both', which='major', labelsize=fz)
    
    ax1b = ax1.twiny()
    col2 = 'blue'
    ax1b.plot(tau_xmean, yy, linewidth=2, color=col2)
    ax1b.set_xlabel("Zonal wind stress (N/m$^2$)", fontsize=fz, color=col2)
    ax1b.tick_params(axis='x', labelcolor=col2)
    ax1b.set_xlim(xlim2)
    ax1b.tick_params(axis='both', which='major', labelsize=fz)
    
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    
    axes = [ax1, ax1b]
    for ii in range(2):
        axes[ii].set_ylabel("Y (km)", fontsize=fz)
        
        axes[ii].set_ylim(0, 2500)
        axes[ii].tick_params(axis='both', which='major', labelsize=fz)
        axes[ii].axvline(0, linestyle='--', color='k')
        
    # plot bathy
    if len(exp_name_alias)==0:
        exp_name_str = exp_name 
    else:
        exp_name_str = exp_name_alias
        
    ax2 = plt.subplot(gs[0, 1:])
    plt.sca(ax2)
    ZZ[-1:, :] = np.nan
    XX, YY = np.meshgrid(xx, yy)

    im = plt.contourf(XX, YY, ZZ, 30, cmap=plt.cm.inferno, vmin=-4000, vmax=0)
    cb = plt.colorbar(im)
    cb.set_label('Depth (m)', fontsize=fz)
    plt.xlabel('X (km)', fontsize=fz)
    #plt.ylabel('Y (km)', fontsize=fz)
    plt.title("%s bathymetry" %exp_name_str, fontsize=fz)
    ax2.tick_params(axis='both', which='major', labelsize=fz)
    plt.ylim(0, 2500)
    
    plt.subplots_adjust(wspace=0.4)
    
    if save_plots:
        fname = '%s_surface_forcing_bathy.pdf' %exp_name
   
    
def plot_bathy(exp_names, tr, xr, zr, exp_names_alias=[], ylim=[], fz=14, cmap=plt.cm.viridis, save_plots=True):
        
    params, forcing = hf.get_exp_params(exp_name)
    
    # This import registers the 3D projection, but is otherwise unused.
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
    from mpl_toolkits import mplot3d
    # from mayavi import mlab
    
    
    fig = plt.figure(figsize=(15, 10))
    ax = fig.gca(projection='3d')
    
    XX, YY = np.meshgrid(params['xx']/1e3, params['yy']/1e3)
    bathy = forcing['bathy'].T
    
    bathy[:, -2:] = -4000
    # bathy[-1, :] = -4000
    # print(bathy.shape)
    # Plot the surface.
    #surf = ax.plot_surface(XX,  YY, bathy, cmap=plt.cm.inferno, vmin=-4000, vmax=0, rstride=2, cstride=10, edgecolor='none')
    #surf = ax.contourf3D(XX, YY, bathy, 50, cmap=plt.cm.inferno)
    ss = 5
    surf = ax.plot_trisurf(XX.flatten()[::ss],  YY.flatten()[::ss], bathy.flatten()[::ss], 
                       cmap=plt.cm.inferno, vmin=-4000, vmax=0,  edgecolor='none')

    #mlab.surf(XX, YY, bathy, warp_scale='auto')

    # Customize the z axis.
#     ax.set_zlim(-1.01, 1.01)
#     ax.zaxis.set_major_locator(LinearLocator(10))
#     ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    # Add a color bar which maps values to colors.
    cb = fig.colorbar(surf, shrink=0.4, aspect=10)
    cb.set_label("Depth (m)", fontsize=fz)
    
    ax.set_ylabel("\n X (km)", fontsize=fz)
    ax.set_xlabel("\n Y (km)", fontsize=fz)
    ax.set_zlabel("\n Z (m)", fontsize=fz)
    #ax.set_title("%s bathymetry" %exp_name, fontsize=fz+2)
    ax.view_init(30, -70)
    
    plt.show()
    #mlab.show()
    
    fname = '3d_test.pdf'
    #plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')        
        
        
        
        
        
        