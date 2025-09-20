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

fz = 14


def plot_gyre_temp_change(theta_ds, dep_r, ds_type, gyres_coords, params, fz=14, xlim=[]):
    
    if ds_type=='inst':
        plot_sm = True
        lw = 0.5
        msz=4
    elif ds_type=='avg':
        plot_sm = False
        lw=2
        msz=8
        
    secsInYr = 3600*24*365 
    tt = theta_ds['time']*params['deltaT']/secsInYr
    
    plt.figure(figsize=(10, 5))
    
    for gyre_name in gyres_coords:
        coords = gyres_coords[gyre_name]
        
        zz = np.abs(theta_ds['Z'].values)
        theta_gyre = theta_ds.isel(XC=np.logical_and(theta_ds['XC']>coords['x'][0], 
                                                     theta_ds['XC']<coords['x'][-1]),
                                      YC=np.logical_and(theta_ds['YC']>coords['y'][0], 
                                                        theta_ds['YC']<coords['y'][-1]), 
                                      Z=np.logical_and(zz>=dep_r[0], zz<dep_r[-1]))


        theta_gyre_xymean = theta_gyre.mean(dim=['XC', 'YC', 'Z'])
        theta_gyre_xymean_sm = theta_gyre_xymean.rolling(time=3, center=True).mean()

        
        p = plt.plot(tt, theta_gyre_xymean-theta_gyre_xymean[0], '-', marker='o', lw=lw, markersize=msz, label=gyre_name)
        if plot_sm:
            plt.plot(tt, theta_gyre_xymean_sm-theta_gyre_xymean[0], '-',color=p[0].get_color(), lw=2)


    if len(xlim)==2:
        plt.xlim(xlim)
    plt.ylabel("Temperature ($^{\circ}$C)", fontsize=fz)
    plt.xlabel('Time (year)', fontsize=fz)
    plt.title("%s-%s m temperature change" %(dep_r[0], dep_r[1]), fontsize=fz+2)
    plt.grid(True)
    plt.legend(loc=0, fontsize=fz)
    plt.grid(True)
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)
    plt.axhline(0, linestyle='--', color='k', lw=2)

def save_zonal_sect(theta_ds, exp_name, clvls, movie_dirname, x_range=[], cmap=plt.cm.viridis, 
                    vname='THETA', dpi=100, bg_col='k', testing=True):
    
    movie_dir = 'movies/%s/' %movie_dirname
    os.makedirs(movie_dir, exist_ok=True)
    
    params, _ = hf.get_exp_params(exp_name, vname=vname)
    
    if len(x_range)==2:
        theta_sect = theta_ds.isel(XC=np.logical_and(theta_ds['XC']>x_range[0], theta_ds['XC']<x_range[-1]))
        theta_sect_mean = theta_sect.mean(dim=['XC'])
    else:
        theta_sect_mean = theta_ds.mean(dim=['XC'])

    secsInYr = 3600*24*365 
    tvec = theta_ds['time']*params['deltaT']/secsInYr
    
    vmin = clvls.min()
    vmax = clvls.max()
    
    clines = clvls[::8]
    
#     yy = theta_ds['YC']/1e3
#     zz = theta_ds['Z']
    
    ZZ, YY = np.meshgrid(params['zz'], params['yy']/1000)
    
    for tt in range(len(tvec)): #range(len(tvec)):
        fig = plt.figure(figsize=(1000/dpi, 600/dpi), dpi=dpi)
        
        if tt%5==0:
            print("t = %s/%s" %(tt, len(tvec)))
                  
        theta_sect_mean_tt = np.ma.masked_where(theta_sect_mean[tt,:, :]==0, theta_sect_mean[tt,:, :])
        theta_sect_mean_anom = theta_sect_mean_tt - theta_sect_mean[0,:, :]
        #theta_sect_mean_anom = np.ma.masked_where(theta_sect_mean_anom==0, theta_sect_mean_anom)
        #debug_here()
        im = plt.contourf(YY, ZZ, theta_sect_mean_anom.T, clvls, cmap=cmap, extend='both')
        cb = plt.colorbar(im, extend='both', ticks=clines)
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        plt.title("Temperature anomaly (year %s)" %(np.round(tvec[tt].values)), fontsize=fz+2)
        
        plt.xlabel("Y (km)", fontsize=fz)
        plt.ylabel("Z (m)", fontsize=fz)
        plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        plt.gca().set_facecolor(bg_col)
        
        fpath = os.path.join(movie_dir, 'frame_%05d.png' %tt)
        plt.savefig(fpath, dpi=dpi)
        
        plt.close('all')
        
        if testing and tt>=5:
            break
    
    
def save_layer_avg_movie(vble_zavg_ds, exp_name, clvls, movie_dirname, cmap=plt.cm.viridis, 
                         vname='', vble_units='', dpi=100, bg_col='k', testing=True, 
                         time_units='years', skip_first_step=True, mask=False, cb_extend='both', max_frames=200):   
    
    movie_dir = 'movies/%s/' %movie_dirname
    os.makedirs(movie_dir, exist_ok=True)
    
    params, _ = hf.get_exp_params(exp_name, vname='THETA')
    
    secsInDay = 3600*24
    secsInYr = 3600*24*365
    
    if time_units=='days':
        tvec = (vble_zavg_ds['time'] - vble_zavg_ds['time'][0])*params['deltaT']/secsInDay 
    elif time_units=='years':
        tvec = vble_zavg_ds['time']*params['deltaT']/secsInYr
    
    vmin = clvls.min()
    vmax = clvls.max()
    
    XX, YY = np.meshgrid(params['xx']/1000, params['yy']/1000)
    
    if len(mask)==len(vble_zavg_ds[0, :, :].values):
        
        use_mask=True
    else:
        use_mask=False
    
    for tt in range(len(tvec)): #range(len(tvec)):
        
        if tt==0 and skip_first_step:
            t_off = 0
            continue
        else:
            t_off = 1

            
        fig = plt.figure(figsize=(1200/dpi, 600/dpi), dpi=dpi)
        if tt%5==0:
            print("t = %s/%s" %(tt, len(tvec)))
            
        if use_mask:
            vble_ma = np.ma.masked_where(mask, vble_zavg_ds[tt, :, :].values)
            
        else:
            vble_ma = np.ma.masked_where(vble_zavg_ds[tt, :, :].values==0, vble_zavg_ds[tt, :, :].values)
            
            
        im = plt.contourf(XX, YY, vble_ma, clvls, vmin=vmin, vmax=vmax, cmap=cmap, extend=cb_extend)
        plt.xlabel("X (km)", fontsize=fz)
        plt.ylabel("Y (km)", fontsize=fz)
        plt.gca().set_facecolor(bg_col) 
        plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        plt.xticks(np.arange(-2000, 2001, 500))
        plt.yticks(np.arange(0, 2501, 500))
        
                   
        
        plt.suptitle("%s (day %i)" %(vname,tvec[tt]+t_off), fontsize=fz+4, y=0.92)
        plt.subplots_adjust(hspace=0.3)
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        cb = plt.colorbar(im, cax=cbar_ax, extend=cb_extend)
        cb.set_label('%s (%s)' %(vname, vble_units), fontsize=fz)
        #plt.sca(axes[ii])
        
        fpath = os.path.join(movie_dir, 'frame_%05d.png' %tt)
        plt.savefig(fpath, dpi=dpi)
        
        plt.close('all')
        
        if testing and tt>=5:
            break
            
        if tt>=max_frames:
            
            print("Maximum frames saved. Exiting...")
            break
            
    
    
    
def save_layer_avg_movie_combo(vble_zavg_ds_list, exp_name, clvls, movie_dirname, cmap=plt.cm.viridis, 
                         vname='', vble_units='', dpi=100, bg_col='k', testing=True, 
                         time_units='years', skip_first_step=True, mask=False, cb_extend='both', max_frames=200):   
    
    #movie_dir = 'movies/%s/' %movie_dirname
    os.makedirs(movie_dirname, exist_ok=True)
    
    params, _ = hf.get_exp_params(exp_name, vname='VVEL')
    
    secsInDay = 3600*24
    secsInYr = 3600*24*365
    
    if time_units=='days':
        tvec = (vble_zavg_ds_list[0]['time'] - vble_zavg_ds_list[0]['time'][0])*params['deltaT']/secsInDay 
    elif time_units=='years':
        tvec = vble_zavg_ds_list[0]['time']*params['deltaT']/secsInYr
    
    vmin = clvls.min()
    vmax = clvls.max()
    
    XX, YY = np.meshgrid(params['xx']/1000, params['yy']/1000)
    
    if len(mask)==len(vble_zavg_ds_list[0][0, :, :].values):
        
        use_mask=True
    else:
        use_mask=False
    
    for tt in range(len(tvec)): #range(len(tvec)):
        
        if tt==0 and skip_first_step:
            t_off = 0
            continue
        else:
            t_off = 1

        if tt%5==0:
            print("t = %s/%s" %(tt, len(tvec)))
            
        
        fig, axes = plt.subplots(2, 1, figsize=(1200/dpi, 800/dpi), dpi=dpi)
        
        for ii, vble_zavg_ds in enumerate(vble_zavg_ds_list):
            
            plt.sca(axes[ii])
            if use_mask:
                vble_ma = np.ma.masked_where(mask, vble_zavg_ds[tt, :, :].values)
            else:
                vble_ma = np.ma.masked_where(vble_zavg_ds[tt, :, :].values==0, vble_zavg_ds[tt, :, :].values)

            im = plt.contourf(XX, YY, vble_ma, clvls, vmin=vmin, vmax=vmax, cmap=cmap, extend=cb_extend)
            plt.xlabel("X (km)", fontsize=fz)
            plt.ylabel("Y (km)", fontsize=fz)
            plt.gca().set_facecolor(bg_col) 
            plt.gca().tick_params(axis='both', which='major', labelsize=fz)
            plt.xticks(np.arange(-2000, 2001, 500))
            plt.yticks(np.arange(0, 2501, 500))

        plt.suptitle("%s (day %i)" %(vname,tvec[tt]+t_off), fontsize=fz+4, y=0.92)
        plt.subplots_adjust(hspace=0.3)
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        cb = plt.colorbar(im, cax=cbar_ax, extend=cb_extend)
        cb.set_label('%s (%s)' %(vname, vble_units), fontsize=fz)
        #plt.sca(axes[ii])
        
        fpath = os.path.join(movie_dirname, 'frame_%05d.png' %tt)
        #print(fpath)
        plt.savefig(fpath, dpi=dpi)
        
        plt.close('all')
        
        if testing and tt>=5:
            break
            
        if tt>=max_frames:
            
            print("Maximum frames saved. Exiting...")
            break   
    
    
    
    
    
    
    
    
    
    
    
    
    
