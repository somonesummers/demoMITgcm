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
import run_utils as ru


from IPython.core.debugger import set_trace
debug_here = set_trace

imp.reload(hf)
imp.reload(pf)
imp.reload(ru)


home_dir = '/home/earlew/research/scripts/MITgcm_py/'
plot_dir = os.path.join(home_dir, 'plots_temporary/')
output_dir = os.path.join(home_dir, 'output/')

def plot_layer_avg_anom(exp_names, vname, t0, t1, zr, exp_names_alias='', clvls=[], clvls_ano=[], fz=14, 
                        cmap=plt.cm.viridis, lcol='w', bg_col='k',add_clines=True, save_plots=False, cls=4,
                        compareFirstLast=True):
    
    """
    function to plot and compare layer-averaged quantities at two different time instances
    """
    
    printLoadStatus = False
    cbar_ext = 'both'
    
    if len(exp_names)==1:
        fig, axes = plt.subplots(3, 1, figsize=(12, 15))
        axes = axes[:, np.newaxis]
    else:
        fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    
    for ii, exp_name in enumerate(exp_names):
        
        exp_name_alias = exp_names_alias[ii]
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)

        if compareFirstLast:
            print("comparing first and last years. Ignoring t0, t1")
            yrs = ru.getSavedYrs(exp_name)
            yrs = np.round(yrs).astype(int)
            
            tr0 = [yrs[0], yrs[0]]
            tr1 = [yrs[-1], yrs[-1]]
            
        else:
            tr0 = [t0, t0]
            tr1 = [t1, t1]  
      
        
        
        title_str = 'Mean %s' %vname
        
        # load output
        if vname.lower() == 'ssh':

            # load variable at t=t0
            vble_t0 = hf.readIters(vname='PHIHYD', tr=tr0, params=params)
            vble_t0 = vble_t0/params['gravity']
            theta = hf.readIters(vname='THETA', tr=tr0, params=params, printStatus=printLoadStatus)
            vble_t0 = np.ma.masked_where(theta.mask, vble_t0)

            # load variable at t=t1
            vble_t1 = hf.readIters(vname='PHIHYD', tr=tr1, params=params)
            vble_t1 = vble_t1/params['gravity']
            vble_t1 = np.ma.masked_where(theta.mask, vble_t1)
            
            vble_t0_zavg = vble_t0[...,0] # not really a zavg. Just trying to keep the names consistent
            vble_t1_zavg = vble_t1[...,0]
            
        elif vname == 'BT_psi':
            theta = hf.readIters(vname='THETA', tr=tr0, params=params, printStatus=printLoadStatus)
            vble_t0_zavg = hf.calc_BT_streamfunction(exp_name, tr0)
            vble_t1_zavg = hf.calc_BT_streamfunction(exp_name, tr1)
            
            vble_t0_zavg = np.ma.masked_where(theta[:, :, 0].mask, vble_t0_zavg[1:, 1:]/1e6)
            vble_t1_zavg = np.ma.masked_where(theta[:, :, 0].mask, vble_t1_zavg[1:, 1:]/1e6)
            
            #vble_t0_zavg = vble_t0_zavg[1:, 1:]/1e6
            #vble_t1_zavg = vble_t1_zavg[1:, 1:]/1e6
            
        elif vname in ['SPD', 'SPD_inst']:
            if vname.endswith('inst'):
                ext = '_inst'
            else:
                ext = ''
                
            params, _ = hf.get_exp_params(exp_name, vname='UVEL'+ext)
            
            vble_t0a = hf.readIters(vname='UVEL'+ext, tr=tr0, params=params, printStatus=printLoadStatus)
            vble_t0b = hf.readIters(vname='VVEL'+ext, tr=tr0, params=params, printStatus=printLoadStatus)
            vble_t0 = np.ma.sqrt(vble_t0a**2 + vble_t0b**2)
            vble_t0_zavg = hf.get_layer_mean(vble_t0, zr, vname, params)
            
            vble_t1a = hf.readIters(vname='UVEL'+ext, tr=tr1, params=params, printStatus=printLoadStatus)
            vble_t1b = hf.readIters(vname='VVEL'+ext, tr=tr1, params=params, printStatus=printLoadStatus)
            vble_t1 = np.ma.sqrt(vble_t1a**2 + vble_t1b**2)
            vble_t1_zavg = hf.get_layer_mean(vble_t1, zr, vname, params)   
            
            cbar_ext = 'max'
            
        else:

            vble_t0 = hf.readIters(vname=vname, tr=tr0, params=params)
            vble_t1 = hf.readIters(vname=vname, tr=tr1, params=params)
            
            print(vble_t0.shape)
            
            vble_t0_zavg = hf.get_layer_mean(vble_t0, zr, vname, params)
            vble_t1_zavg = hf.get_layer_mean(vble_t1, zr, vname, params)
            
            title_str = 'Mean %s-%sm %s' %(zr[0], zr[1], vname)
            
            
        vble_zavg_ano = vble_t1_zavg - vble_t0_zavg

        # plot vble at t0, t1, and their difference

        YY,XX = np.meshgrid(params['yy']/1000, params['xx']/1000)

        if len(clvls)==0:
            clvls = np.linspace(-0.9*vble_t0_zavg.max(), 0.9*vble_t0_zavg.max(), 25)

        if len(clvls_ano)==0:
            clvls_ano = np.linspace(-0.9*vble_zavg_ano.max(), 0.9*vble_zavg_ano.max(), 25)

        if len(exp_name_alias)==0:
            exp_name_str = exp_name 
        else:
            exp_name_str = exp_name_alias

        clines = clvls[::cls]

        plt.sca(axes[0, ii])
        im = plt.contourf(XX, YY, vble_t0_zavg, clvls, cmap=cmap, extend=cbar_ext)
        cb = plt.colorbar(im, extend=cbar_ext, ticks=clines)
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        plt.title("%s\n%s at year %s" %(exp_name_str, title_str, t0), fontsize=fz+2)
        plt.text(0.95*XX.min(), 0.9*YY.max(), "(A)", fontsize=fz+4, color='k')

        plt.sca(axes[1, ii])
        im = plt.contourf(XX, YY, vble_t1_zavg, clvls, cmap=cmap, extend=cbar_ext)
        cb = plt.colorbar(im, extend=cbar_ext, ticks=clines)
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        plt.title("%s at year %s" %(title_str, t1), fontsize=fz+2)
        plt.text(0.95*XX.min(), 0.9*YY.max(), "(B)", fontsize=fz+4, color='k')

        plt.sca(axes[2, ii])
        im = plt.contourf(XX, YY, vble_zavg_ano, clvls_ano, cmap=plt.cm.bwr, extend='both')
        cb = plt.colorbar(im, extend='both')
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        plt.text(0.95*XX.min(), 0.9*YY.max(), "(C)", fontsize=16, color='k')
        plt.title("Year %s $-$ Year %s" %(t1, t0), fontsize=fz+2)

        if add_clines:
            cs = axes[0, ii].contour(XX, YY, vble_t0_zavg, clines, colors=lcol, linewidths=1)
            plt.gca().clabel(cs, inline=1, fontsize=fz-2, fmt='%.1f')

            cs = axes[1, ii].contour(XX, YY, vble_t1_zavg, clines, colors=lcol, linewidths=1)
            plt.gca().clabel(cs, inline=1, fontsize=fz-2, fmt='%.1f')


    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    for ax in axes.flatten():
        plt.sca(ax)
        # add labels
        plt.xlabel("X (km)", fontsize=fz)
        plt.ylabel("Y (km)", fontsize=fz)
        # set label size
        plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        # change BG color
        ax.set_facecolor(bg_col)

    if save_plots:
        if vname in ['SSH']:
            fname = '%s_diff_yrs%s-%s.pdf' %(vname, tr1, tr0)
        else:
            fname = 'mean_%s-%sm_%s_diff_yrs%s-%s.pdf' %(zr[0], zr[1], vname, tr1, tr0)

        plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')
        
        
        
def plot_zonal_avg_anom(exp_names, vname, t0, t1, xr, exp_names_alias=[], clvls=[], clvls_ano=[], fz=14, cls=4, cls_cbar=4,
                        cmap=plt.cm.viridis, lcol='w', bg_col='k',add_clines=True, show_diff=True, save_plots=False,
                        add_title=True, fig_w=12, add_text_label=True, compareFirstLast=True):


    if show_diff:
        nr = 3
        figH = 15
    else:
        nr = 2
        figH = 10
        
    if len(exp_names)==1:
        fig, axes = plt.subplots(nr, 1, figsize=(fig_w, figH))
        axes = axes[:, np.newaxis]
    else:
        fig, axes = plt.subplots(nr, 2, figsize=(fig_w+3, figH))
    
    for ii, exp_name in enumerate(exp_names):
        
        #exp_name_alias = exp_names_alias[ii]
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)
        
        
        if compareFirstLast:
            print("comparing first and last years. Ignoring t0, t1")
            yrs = ru.getSavedYrs(exp_name)
            yrs = np.round(yrs).astype(int)
            
            tr0 = [yrs[0], yrs[0]]
            tr1 = [yrs[-1], yrs[-1]]
            
        else:
            tr0 = [t0, t0]
            tr1 = [t1, t1]  
      
        
        vble_t0 = hf.readIters(vname=vname, tr=tr0, params=params)
        vble_t1 = hf.readIters(vname=vname, tr=tr1, params=params)
        
        vble_t0_xavg = hf.get_zonal_mean(vble_t0, xr, vname, params)
        vble_t1_xavg = hf.get_zonal_mean(vble_t1, xr, vname, params)
        vble_xavg_ano = vble_t1_xavg - vble_t0_xavg
        
        #title_str = 'Mean %s-%sm %s' %(xr[0]/1e3, xr[1]/1e3, vname)
        title_str = '%s' %(vname)
        
        ZZ, YY = np.meshgrid(params['zz'], params['yy']/1000)
        xr_mean = np.array(xr).mean()/1e3

        if len(clvls)==0:
            clvls = np.linspace(-0.9*vble_t0_xavg.max(), 0.9*vble_t0_xavg.max(), 25)

        if len(clvls_ano)==0:
            clvls_ano = np.linspace(-0.9*vble_xavg_ano.max(), 0.9*vble_xavg_ano.max(), 25)

        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 

        clines = clvls[::cls]
        cbar_ticks = clvls[::cls_cbar]
        
        plt.sca(axes[0, ii])
        im = plt.contourf(YY, ZZ, vble_t0_xavg, clvls, cmap=cmap, extend='both')
        cb = plt.colorbar(im, extend='both', ticks=cbar_ticks)
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        if add_title:
            plt.title("%s\n%s at year %s" %(exp_name_str, title_str, t0), fontsize=fz+2)
        
        #plt.text(0.02*YY.max(), 0.95*ZZ.min(), "(A)", fontsize=fz+4, color='k')

        plt.sca(axes[1, ii])
        im = plt.contourf(YY, ZZ, vble_t1_xavg, clvls, cmap=cmap, extend='both')
        cb = plt.colorbar(im, extend='both', ticks=cbar_ticks)
        cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        if add_title:
            plt.title("%s at year %s" %(title_str, t1), fontsize=fz+2)
        #plt.text(0.02*YY.max(), 0.95*ZZ.min(), "(B)", fontsize=fz+4, color='k')

        if show_diff:
            plt.sca(axes[2, ii])
            im = plt.contourf(YY, ZZ, vble_xavg_ano, clvls_ano, cmap=plt.cm.seismic, extend='both')
            cb = plt.colorbar(im, extend='both')
            cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
            cb.ax.tick_params(labelsize=fz) 
            #plt.text(0.02*YY.max(), 0.95*ZZ.min(), "(C)", fontsize=fz+4, color='k')
            #plt.title("Year %s $-$ Year %s" %(t1, t0), fontsize=fz+2)
            

        if add_clines:
            cs = axes[0, ii].contour(YY, ZZ, vble_t0_xavg, clines, colors=lcol, linewidths=0.75)
            plt.gca().clabel(cs, inline=1, fontsize=fz-3, fmt='%.1f')

            cs = axes[1, ii].contour(YY, ZZ, vble_t1_xavg, clines, colors=lcol, linewidths=0.75)
            plt.gca().clabel(cs, inline=1, fontsize=fz-3, fmt='%.1f')
        
        
        if add_text_label:
            axes[0, ii].text(0.02*YY.max(), 0.95*ZZ.min(), "(A)", fontsize=fz+4, color='k')
            axes[1, ii].text(0.02*YY.max(), 0.95*ZZ.min(), "(B)", fontsize=fz+4, color='k')
            axes[2, ii].text(0.02*YY.max(), 0.95*ZZ.min(), "(C)", fontsize=fz+4, color='k')
            axes[2, ii].set_title("B $-$ A", fontsize=fz+2)
            
        
        plt.subplots_adjust(hspace=0.4, wspace=0.3)
        
        
        

    for ax in axes.flatten():
        plt.sca(ax)
        # add labels
        plt.xlabel("Y (km)", fontsize=fz)
        plt.ylabel("Z (m)", fontsize=fz)
        # set label size
        plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        # change BG color
        ax.set_facecolor(bg_col)
        
    if save_plots:
        fname = 'zonal_mean_%s_x%skm_diff_yrs%s-%s.pdf' %(vname, xr_mean, t1, t0)
        plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')

        
def plot_psi_z_comp(exp_names, t0, t1,  exp_names_alias=[], vname='psi_z', clvls=[], clvls_ano=[], 
                    fz=14, ymax=2400, t_offset=[0, 0], add_isotherms=False, isotherms=[3, 4, 5], 
                    save_plots=False):
    
    vname_list = ['psi_z', 'psi_mean_z', 'psi_eddy_z']
    if vname not in vname_list:
        print("vname should be in %s" %vname_list)
        return
    
    if len(exp_names)==1:
        fig, axes = plt.subplots(3, 1, figsize=(12, 15))
        axes = axes[:, np.newaxis]
    else:
        fig, axes = plt.subplots(3, 2, figsize=(15, 15))
        
    psi_output_dir = os.path.join(home_dir, 'psi_output/')
    
    for ii, exp_name in enumerate(exp_names):
        
        exp_name_alias = exp_names_alias[ii]
        
        #load params
        params, forcing = hf.get_exp_params(exp_name)
        
        bathy = forcing['bathy']
        x_idx = np.argmin(np.abs(params['xx']-0))
        depth_profile = bathy[:, x_idx]
    
        tr0_r = [t0, t0]
        tr1_r = [t1, t1]  
        
        # load overturning data
        #outdir = 'psi_output'
        fname = params['exp_name'] + '_yrs_%s-%s.nc' %(t0-t_offset[0], t0)
        ds0 = xr.open_dataset(os.path.join(psi_output_dir, fname))
        
        fname = params['exp_name'] + '_yrs_%s-%s.nc' %(t1-t_offset[1], t1)
        ds1 = xr.open_dataset(os.path.join(psi_output_dir, fname))
        
        ZZ_psi, YY_psi = np.meshgrid(ds0['zz_f'], params['yy']/1000)
        
        # adjust psi (copied from AS code).
        # EW: better to do all this in the original calculations
        ds_list = [ds0, ds1] 
        for ds in ds_list:
            for jj in range(params['Ny']):
                # Adjust height of top point
                ZZ_psi[jj, 0] = 0

                # Calculate depth of bottom cell  
                hFacS_col = ds['hFacS_f'][0, jj, :]
                kmax = len(hFacS_col[hFacS_col>0])
                if kmax >0:
                    tmp_val = ds0['delRf']*hFacS_col
                    ZZ_psi[jj, kmax-1] = -tmp_val.sum()

                # Force streamfunction to be zero at boundaries
                ds['psi_z'][jj, 0] = 0
                ds['psi_mean_z'][jj, 0] = 0
                ds['psi_eddy_z'][jj, 0] = 0

                if kmax>0:
                    ds['psi_z'][jj, kmax-1] = 0
                    ds['psi_mean_z'][jj, kmax-1] = 0
                    ds['psi_eddy_z'][jj, kmax-1] = 0   

                    ds['psi_z'][jj, kmax-1:] = np.nan
                    ds['psi_mean_z'][jj, kmax-1:] = np.nan
                    ds['psi_eddy_z'][jj, kmax-1:] = np.nan  
        
        if len(clvls)==0:
            clvls = np.linspace(-0.9*ds0['psi_z'].max(), 0.9*ds0['psi_z'].max(), 25)
         
        psi_diff = ds1[vname]-ds0[vname]
        if len(clvls_ano)==0:
            clvls_ano = np.linspace(-0.9*psi_diff.max(), 0.9*psi_diff.max(), 25)
             
        #vble_list = ['psi_z', 'psi_mean_z', 'psi_eddy_z']
        #vname_list = ['$\Psi_r$', '$\Psi_m$', '$\Psi_e$']
        
        vname_Tex = {'psi_z':'$\Psi_r$', 'psi_mean_z':'$\Psi_m$', 'psi_eddy_z':'$\Psi_e$'}
        
        if len(exp_name_alias)==0:
            exp_name_str = exp_name 
        else:
            exp_name_str = exp_name_alias
            
        if add_isotherms:
            x_r = [-2000e3, 2000e3]
            clines_theta = np.arange(1, 5.01, 1)
            
            #params, _ = hf.get_exp_params(exp_name[0])
            ZZ,YY = np.meshgrid(params['zz'], params['yy']/1000)
            theta_ctrl = hf.readIters(vname='THETA', tr=tr0_r, params=params)
            theta_xavg_ctrl = hf.get_zonal_mean(theta_ctrl, x_r, 'THETA', params)

            theta_ptrb = hf.readIters(vname='THETA', tr=tr1_r, params=params)
            theta_xavg_ptrb = hf.get_zonal_mean(theta_ptrb, x_r, 'THETA', params)
            
            iso_col = 'DarkGreen' #w
            
            
        
        cl_min = clvls.min()
        step = abs(cl_min)/2
        clines = np.arange(cl_min, abs(cl_min), step)
    
        cbar_ticks = clvls[::4]

        #------------ Plot Psi t=t0 ----------------#
        plt.sca(axes[0, ii])
        im = plt.contourf(YY_psi, ZZ_psi, ds0[vname], clvls, cmap=plt.cm.bwr, extend='both')
        plt.contour(YY_psi, ZZ_psi, ds0[vname], [0], colors='k', linewidths=2)
        plt.xlim(right=ymax)
        #plt.xlabel('y (km)', fontsize=fz)
        plt.ylabel('Depth (m)', fontsize=fz)
        axes[0, ii].set_facecolor('k')  
        
        
        plt.title("%s\n %s at year %s" %(exp_name_str, vname_Tex[vname], t0), fontsize=fz+2)
        plt.plot(params['yy']/1e3, depth_profile, '--', lw=2, color='k' )
        #plt.text(0.01*YY_psi.max(), 0.95*ZZ_psi.min(), "%s$_1$" %(vname), fontsize=fz+4, color='w')
        #plt.grid(True)
        
        cb = plt.colorbar(im, extend='both', ticks=cbar_ticks)
        cb.set_label('$\Psi$ (Sv)', fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        
        if add_isotherms:
            #debug_here()
            cs = plt.contour(YY, ZZ, theta_xavg_ctrl, isotherms, colors=iso_col, linewidths=1)
            #cs = plt.contour(YY, ZZ, theta_xavg_ctrl, clines_theta, colors='k', linewidths=0.5)
            plt.gca().clabel(cs, inline=1, fontsize=fz-4, fmt='%.1f')
            

        #------------ Plot Psi t=t1 ----------------#
        plt.sca(axes[1, ii])
        im = plt.contourf(YY_psi, ZZ_psi, ds1[vname], clvls, cmap=plt.cm.bwr, extend='both')
        plt.contour(YY_psi, ZZ_psi, ds1[vname], [0], colors='k', linewidths=2)
        plt.xlim(right=ymax)
        #plt.xlabel('y (km)', fontsize=fz)
        plt.ylabel('Depth (m)', fontsize=fz)
        axes[1, ii].set_facecolor('k')  
        plt.title("%s at year %s" %(vname_Tex[vname], t1), fontsize=fz+2)
        plt.plot(params['yy']/1e3, depth_profile, '--', lw=2, color='k' )
        
        #plt.text(0.01*YY_psi.max(), 0.95*ZZ_psi.min(), "%s$_2$" %(vname), fontsize=fz+4, color='w')
        cb = plt.colorbar(im, extend='both', ticks=cbar_ticks)
        cb.set_label('$\Psi$ (Sv)', fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        
        if add_isotherms:
            cs = plt.contour(YY, ZZ, theta_xavg_ptrb, isotherms, colors=iso_col, linewidths=1)
            plt.gca().clabel(cs, inline=1, fontsize=fz-4, fmt='%.1f')

        #------------ Plot difference ----------------#
        plt.sca(axes[2, ii])
        im = plt.contourf(YY_psi, ZZ_psi, psi_diff, clvls_ano, cmap=plt.cm.bwr, extend='both')
        plt.xlim(right=ymax)
        plt.xlabel('y (km)', fontsize=fz)
        plt.ylabel('Depth (m)', fontsize=fz)
        cb.set_label('$\Psi$ (Sv)', fontsize=fz)
        cb.ax.tick_params(labelsize=fz)
        #plt.grid(True)
        #
        axes[2, ii].set_facecolor('k')  
        #axes[2, ii].axvline(1000, linestyle=':', linewidth=2, color='k')
        plt.title("Year %s $-$ Year %s" %(t1, t0), fontsize=fz+2)
        
        cb = plt.colorbar(extend='both')
        cb.set_label('$\Psi$ (Sv)', fontsize=fz)
        cb.ax.tick_params(labelsize=fz) 
        
        #text_str = "%s$_2$ $-$ %s$_1$" %(vname, vname)
        #plt.text(0.01*YY_psi.max(), 0.95*ZZ_psi.min(), text_str, fontsize=fz+4, color='w')
    
    axes = axes.flatten()    
    for ii,ax in enumerate(axes):
        ax.axvline(1250, linestyle=':', linewidth=1, color='k')
        ax.axvline(1750, linestyle=':', linewidth=1, color='k')
        
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
#     if save_plots:
#         fname = 'zonal_mean_psi_z_comp_yrs%s-%s_pm%s.pdf' %(tr[0], tr[1], plot_mode)
#         plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')                    


def ZonalDepthMeanEvol(exp_names, vname, tr, xr=[], zr=[], exp_names_alias=[], ylim=[], ylim_ano=[], 
                       add_zero_line=False, fz=14, add_title=True, add_text_label=True):
    
    # WARNING: xr argument conflicts with xarray alias
    
    fig, axes = plt.subplots(2, len(exp_names), figsize=(15, 10))
    if len(exp_names)==1:
        axes = axes[:, np.newaxis]
        
    
    
    secsInYear = 3600*24*365
    Yap_km = 1250
    Ysa_km = 1750
    
    if len(zr)==2:
        title_str = 'Mean %s-%sm %s' %(abs(zr[1]), abs(zr[0]), vname)
    else:
        title_str = 'Mean %s' %(vname)
        
    for ii, exp_name in enumerate(exp_names):
        
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)
        
        # load vble
        vble_tr, tvec = hf.readIters(vname=vname, tr=tr, params=params, returnTimeAvg=False)
        
        # compute zonal
        if len(xr)==0:
            vble_tr_xavg = vble_tr.mean(axis=0)
        else:
            xx = params['xx']
            xi = np.logical_and(xx>=xr[0], xx<=xr[-1])
            vble_tr_xavg = vble_tr[xi, ...].mean(axis=0)
            
        
        if len(zr)==0:
            vble_tr_xzavg = vble_tr_xavg.mean(axis=1)
        else:
            zz = params['zz']
            zi = np.logical_and(zz>=zr[0], zz<=zr[-1])
            vble_tr_xzavg = vble_tr_xavg[:, zi, :].mean(axis=1)
            

        # plot time means
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
            
            
        
        tvec_yrs = np.round(tvec/secsInYear)
        nt = vble_tr_xzavg.shape[1]
        cols = plt.cm.rainbow(np.linspace(0,1,nt))
        vble_tr_xzavg_ano = vble_tr_xzavg-np.atleast_2d(vble_tr_xzavg[:, 0]).T
        for tt in range(nt):
            if tt==0 or tt==(nt-1):
                axes[0, ii].plot(params['yy']/1e3, vble_tr_xzavg[:, tt], color=cols[tt], 
                                 lw=2, label="year %i" %tvec_yrs[tt])
                
                axes[1, ii].plot(params['yy']/1e3, vble_tr_xzavg_ano[:, tt], color=cols[tt], 
                                 lw=2, label="year %i" %tvec_yrs[tt])
            else:
                axes[0, ii].plot(params['yy']/1e3, vble_tr_xzavg[:, tt], color=cols[tt], lw=2)
                axes[1, ii].plot(params['yy']/1e3, vble_tr_xzavg_ano[:, tt], color=cols[tt], lw=2)
        
    
        
        axes[0, ii].legend(loc=0, fontsize=fz-4)
        
        if add_title:
            axes[0, ii].set_title("%s\n%s" %(exp_name_str, title_str), fontsize=fz+2)
            axes[1, ii].set_title("%s anomaly" %(title_str), fontsize=fz+2)
        
        if len(ylim)==2:
            axes[0, ii].set_ylim(ylim)
            
        if len(ylim_ano)==2:
            axes[1, ii].set_ylim(ylim_ano)
            
        if add_zero_line:
            axes[0, ii].axhline(0, linestyle='--', linewidth=2, color='k')  
        
        
        for ax in axes[:, ii]:
            plt.sca(ax)
            plt.grid(True)
            plt.ylabel('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
            plt.xlabel('Y (km)', fontsize=fz)
            plt.xlim(0, 2500)
            ax.axvline(Yap_km, linestyle=':', linewidth=2, color='k')  
            ax.axvline(Ysa_km, linestyle=':', linewidth=2, color='k')  
        

    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    

def VTransProfileEvol(exp_names, tr, vname='VVEL', xr=[], y_km=1000, exp_names_alias=[], xlim=[], xlim_ano=[], 
                       add_zero_line=False, fz=14):
    
    fig, axes = plt.subplots(2, len(exp_names), figsize=(15, 10))
    if len(exp_names)==1:
        axes = axes[:, np.newaxis]
        
    
    secsInYear = 3600*24*365
    Yap_km = 1250
    Ysa_km = 1750
    
    if len(xr)==0:
        title_str = 'Meridional transport at Y=%s km' %(y_km)
    else:
        title_str = 'Meridional transport at Y=%s km' %(y_km)
        
    for ii, exp_name in enumerate(exp_names):
        
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)
        
        # load vble
        vble_tr, tvec = hf.readIters(vname=vname, tr=tr, params=params, mask_zeros=False, returnTimeAvg=False)
        
        # get zonal slice
        y_idx = np.argmin(np.abs(params['yy']-y_km*1e3))
        vble_tr_y = vble_tr[:, y_idx, :, :]
        
        # compute volume flux
        hFacS_y = params['hFacS'][:, y_idx, :].T
        vh_y = vble_tr_y*hFacS_y[:, :, np.newaxis]*params['delR'][np.newaxis,:, np.newaxis]
        
        # compute zonal integral
        if len(xr)==0:
            vtrans_profile = np.ma.sum(vh_y*params['delX'][:, np.newaxis, np.newaxis], axis=0)
        else:
            xx = params['xx']
            xi = np.logical_and(xx>=xr[0], xx<=xr[-1])
            vhx_y = vh_y*params['delX'][:, np.newaxis, np.newaxis]
            vtrans_profile = vhx_y[xi, ...].sum(axis=0)
            

        # plot time means
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
            

        tvec_yrs = np.round(tvec/secsInYear)
        nt = vtrans_profile.shape[1]
        cols = plt.cm.rainbow(np.linspace(0,1,nt))
        
        vtrans_profile_Sv = vtrans_profile/1e6
        vtrans_profile_Sv_ano = vtrans_profile_Sv-np.atleast_2d(vtrans_profile_Sv[:, 0]).T
        
        zz = params['zz']
        
        zz_bins = np.arange(0, 4001, 200)
        zz_bins_c = zz_bins[:-1] + np.diff(zz_bins)/2
        
        #print(vble_tr_xavg_profile_ano.shape)
        
        for tt in range(nt):
            
            vtrans_Sv_binsum = stats.binned_statistic(np.abs(zz), vtrans_profile_Sv[:, tt].compressed(), 
                                                      statistic='sum', bins=zz_bins)[0]
            #vtrans_Sv_binsum = stats.binned_statistic(zz, vtrans_profile_Sv[:, tt], statistic='sum', bins=zz_bins)[0]
            
            if tt==0:
                ref = vtrans_Sv_binsum
                
            vtrans_Sv_binsum_ano = vtrans_Sv_binsum-ref

            if tt==0 or tt==(nt-1):
                axes[0, ii].plot(vtrans_Sv_binsum, -zz_bins_c,'-', color=cols[tt], lw=2, label="year %i" %tvec_yrs[tt])
                #axes[1, ii].plot(vtrans_Sv_binsum_ano, -zz_bins_c, color=cols[tt], lw=2, label="year %i" %tvec_yrs[tt])
            else:
                axes[0, ii].plot(vtrans_Sv_binsum, -zz_bins_c,'-', color=cols[tt], lw=2)
                #axes[1, ii].plot(vtrans_Sv_binsum_ano, -zz_bins_c,'-', color=cols[tt], lw=2)
                
            if tt==(nt-1):
                axes[1, ii].barh(-zz_bins_c, vtrans_Sv_binsum_ano, align='center', height=np.diff(zz_bins)) 
                #debug_here()
                
            
        axes[0, ii].set_title("%s\n%s" %(exp_name_str, title_str), fontsize=fz+2)
        axes[0, ii].legend(loc=0, fontsize=fz-4)
        axes[1, ii].set_title("Year %s $-$ Year %s" %(tr[1], tr[0]), fontsize=fz+2)
        
        
        if len(xlim)==2:
            axes[0, ii].set_xlim(xlim)
            
        if len(xlim_ano)==2:
            axes[1, ii].set_xlim(xlim_ano)
            
        if add_zero_line:
            axes[0, ii].axvline(0, linestyle='--', linewidth=2, color='k')  

        
        for ax in axes[:, ii]:
            plt.sca(ax)
            plt.grid(True)
            #plt.xlabel('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
            plt.xlabel(r'Meridional transport (10$^{6}$ m$^3$/s)', fontsize=fz)
            plt.ylabel('Depth (m)', fontsize=fz)
            plt.ylim(zz.min(), zz.max())
        

    plt.subplots_adjust(hspace=0.4, wspace=0.3)  
    
    
def resTransportProfile(exp_names, t0, t1, vname='psi_z', x_r=[], y_km=1000, exp_names_alias=[], 
                        xlim=[], xlim_ano=[], add_zero_line=False, fz=14):

    import warnings
    warnings.filterwarnings("ignore", category=FutureWarning)
    
    vname_list = ['psi_z', 'psi_mean_z', 'psi_eddy_z']
    vname_Tex = {'psi_z':'$\Psi_r$', 'psi_mean_z':'$\Psi_m$', 'psi_eddy_z':'$\Psi_e$'}
    
    if vname not in vname_list:
        print("vname should be in %s" %vname_list)
        
    else:
        
        title_str = '%s transport at Y=%s km' %(vname_Tex[vname], y_km)
    
    fig, axes = plt.subplots(2, len(exp_names), figsize=(15, 12))
    if len(exp_names)==1:
        axes = axes[:, np.newaxis]
    
    for ii, exp_name in enumerate(exp_names):
        
        exp_name_alias = exp_names_alias[ii]
        
        #load params
        params, _ = hf.get_exp_params(exp_name)
        tr0_r = [t0, t0]
        tr1_r = [t1, t1]  
        
        # load overturning data
        outdir = 'psi_output'
        fname = params['exp_name'] + '_yrs_%s-%s.nc' %(t0, t0)
        ds0 = xr.open_dataset(os.path.join(outdir, fname))
        
        fname = params['exp_name'] + '_yrs_%s-%s.nc' %(t1, t1)
        ds1 = xr.open_dataset(os.path.join(outdir, fname))
        
        ZZ_psi, YY_psi = np.meshgrid(ds0['zz_f'], params['yy']/1000)
        
        # adjust psi (copied from AS code).
        # EW: better to do all this in the original calculations
        ds_list = [ds0, ds1] 
        for ds in ds_list:
            for jj in range(params['Ny']):
                # Adjust height of top point
                ZZ_psi[jj, 0] = 0

                # Calculate depth of bottom cell  
                hFacS_col = ds['hFacS_f'][0, jj, :]
                kmax = len(hFacS_col[hFacS_col>0])
                if kmax >0:
                    tmp_val = ds0['delRf']*hFacS_col
                    ZZ_psi[jj, kmax-1] = -tmp_val.sum()

                # Force streamfunction to be zero at boundaries
                ds['psi_z'][jj, 0] = 0
                ds['psi_mean_z'][jj, 0] = 0
                ds['psi_eddy_z'][jj, 0] = 0

                if kmax>0:
                    ds['psi_z'][jj, kmax-1] = 0
                    ds['psi_mean_z'][jj, kmax-1] = 0
                    ds['psi_eddy_z'][jj, kmax-1] = 0   

#                     ds['psi_z'][jj, kmax-1:] = np.nan
#                     ds['psi_mean_z'][jj, kmax-1:] = np.nan
#                     ds['psi_eddy_z'][jj, kmax-1:] = np.nan  


        # get zonal slice
        y_idx = np.argmin(np.abs(params['yy']-y_km*1e3))
        psi_t0_Y = ds0[vname][y_idx, :]
        psi_t1_Y = ds1[vname][y_idx, :]
        
        p0 = interpolate.interp1d(ds0['zz_f'], psi_t0_Y)
        p1 = interpolate.interp1d(ds0['zz_f'], psi_t1_Y)
        
        z_hr = np.linspace(ds0['zz_f'][0], ds0['zz_f'][-1], 1000)
        psi_t0_Y_hr = p0(z_hr)
        psi_t1_Y_hr = p1(z_hr)
        
        #debug_here()
        vflux_t0_Y = -np.gradient(psi_t0_Y_hr, z_hr, axis=0)
        vflux_t1_Y = -np.gradient(psi_t1_Y_hr, z_hr, axis=0)
        
        
        vTrans_t0_Y = vflux_t0_Y*abs(np.diff(z_hr).mean())
        vTrans_t1_Y = vflux_t1_Y*abs(np.diff(z_hr).mean())
        
        print((vTrans_t1_Y-vTrans_t0_Y).sum())
        
        #debug_here()
        zz = z_hr #ds0['zz_f']
        
        
        zz_bins = np.arange(0, 4001, 500)
        zz_bins_c = zz_bins[:-1] + np.diff(zz_bins)/2
        dz_bin = np.abs(np.diff(zz_bins))
        
        
#         psi_t0_Y = np.ma.masked_invalid(vflux_t0_Y).filled(0)
#         psi_t1_Y = np.ma.masked_invalid(vflux_t1_Y).filled(0)

        
        
        vTrans_t0_Y_binsum = stats.binned_statistic(np.abs(zz), vTrans_t0_Y, statistic='sum', bins=zz_bins)[0]
        vTrans_t1_Y_binsum = stats.binned_statistic(np.abs(zz), vTrans_t1_Y, statistic='sum', bins=zz_bins)[0]
        vTrans_t1_Y_binsum_ano = vTrans_t1_Y_binsum - vTrans_t0_Y_binsum
        
        # plot stuff
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
        
#         axes[0, ii].plot(vTrans_t0_binsum, -zz_bins_c,'-', lw=2, label="year %i" %tvec_yrs[tt])
#         axes[0, ii].plot(vTrans_t1_binsum, -zz_bins_c,'-', lw=2, label="year %i" %tvec_yrs[tt])
        
        bf = 0.4
        axes[0, ii].barh(-zz_bins[1:]+dz_bin*bf, vTrans_t0_Y_binsum, align='edge', height=dz_bin*bf, label="year %i" %t0) 
        axes[0, ii].barh(-zz_bins[1:], vTrans_t1_Y_binsum, align='edge', height=dz_bin*bf, label="year %i" %t1) 
        
        #debug_here()
        
        axes[1, ii].barh(-zz_bins_c, vTrans_t1_Y_binsum_ano, align='center', height=np.diff(zz_bins)) 
        
        axes[0, ii].set_title("%s\n%s" %(exp_name_str, title_str), fontsize=fz+2)
        axes[0, ii].legend(loc=1, fontsize=fz-2)
        axes[1, ii].set_title("Year %s $-$ Year %s" %(t1, t0), fontsize=fz+2)
        
        if len(xlim)==2:
            axes[0, ii].set_xlim(xlim)
            
        if len(xlim_ano)==2:
            axes[1, ii].set_xlim(xlim_ano)
            
        if add_zero_line:
            axes[0, ii].axvline(0, linestyle='--', linewidth=2, color='k')  

        
        for ax in axes[:, ii]:
            plt.sca(ax)
            plt.grid(True)
            #plt.xlabel('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
            plt.xlabel(r'Meridional transport (10$^{6}$ m$^3$/s)', fontsize=fz)
            plt.ylabel('Depth (m)', fontsize=fz)
            plt.ylim(zz.min(), zz.max())
        

    plt.subplots_adjust(hspace=0.3, wspace=0.3)  
        

def plotSSHdiff(exp_names, tr, xr=[-800, 1000], exp_names_alias=[], fz=14):
    
#     fig, axes = plt.subplots(2, len(exp_names), figsize=(15, 10))
    fig, axes = plt.subplots(1,1, figsize=(12, 8))
    if len(exp_names)==1:
        axes = axes[:, np.newaxis]
        
    title_str = 'Meridional SSH difference'
    
    cols = ['cornflowerblue', 'red']
    secsInYear = 3600*24*365
    
    for ii, exp_name in enumerate(exp_names):
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
        
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname='SSH')
        
        # load variable 
        phyd_tr, tvec = hf.readIters(vname='PHIHYD', tr=tr, params=params, mask_zeros=True, returnTimeAvg=False)
        ssh_tr = phyd_tr[:, :, 0, :]/params['gravity']
#         theta = hf.readIters(vname='THETA', tr=[tr[0], tr[0]], params=params)
#         ssh_tr = np.ma.masked_where(theta.mask, ssh_tr)
        ssh_tr = np.ma.masked_invalid(ssh_tr)
        
        tvec_yrs = np.round(tvec/secsInYear)
        
        # compute meridional difference
        if len(xr)==0:
            ssh_xmean_tr = ssh_tr.mean(axis=0)
        else:
            xr_m = np.array(xr)*1000
            xi = np.logical_and(params['xx']>=xr_m[0], params['xx']<xr_m[1])
            ssh_xmean_tr = ssh_tr[xi, :, :].mean(axis=0)
        
        yrS = np.array([500, 1250])*1e3
        yrN = np.array([1750, 2500])*1e3
  
        yiN = np.logical_and(params['yy']>=yrN[0], params['yy']<yrN[1])
        yiS = np.logical_and(params['yy']>=yrS[0], params['yy']<yrS[1])
        #yiS = params['yy']<y0_km*1e3
        
        ssh_N_xymean_tr = ssh_xmean_tr[yiN, :].mean(axis=0)
        ssh_S_xymean_tr = ssh_xmean_tr[yiS, :].mean(axis=0)
        
        ssh_diff_xymean_tr = ssh_N_xymean_tr-ssh_S_xymean_tr
        #tvec_yr = np.arange(len(ssh_N_xymean_tr))*5
        
        #debug_here()
        plt.plot(tvec_yrs, ssh_N_xymean_tr, '--s', color=cols[ii], lw=1, label="%s (north)" %exp_name_str)
        plt.plot(tvec_yrs, ssh_S_xymean_tr, '-->', color=cols[ii], lw=1, label="%s (south)" %exp_name_str)
        plt.plot(tvec_yrs, ssh_diff_xymean_tr, '-', color=cols[ii], lw=2, label="N-S")
        
    
    plt.grid(True)
    plt.xlabel("Years", fontsize=fz)
    plt.ylabel("SSH (m)", fontsize=fz)
    plt.legend(loc=0, fontsize=fz-2, ncol=2)
    plt.axhline(0, linestyle='--', color='k')
    plt.title(title_str, fontsize=fz)
    plt.ylim(-1, 2)
    
    
def plotZonalTransport(exp_names, tr, vname='UVEL', x_km=-1000, yr=[], exp_names_alias=[], fz=14):
    
    
    fig, axes = plt.subplots(1,1, figsize=(12, 8))
    if len(exp_names)==1:
        axes = np.array([axes])
        
    #title_str = 'Zonal transport at X = %i km' %x_km
    
    cols = ['cornflowerblue', 'red', 'green']
    secsInYear = 3600*24*365
    
    for ii, exp_name in enumerate(exp_names):
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
        
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)
    
        # load variable 

        uvel_tr, tvec = hf.readIters(vname=vname, tr=tr, params=params, mask_zeros=True, returnTimeAvg=False)
        tvec_yrs = np.round(tvec/secsInYear)

        
        # get zonal slice
        x_idx = np.argmin(np.abs(params['xx']-x_km*1e3))
        uvel_xx = uvel_tr[x_idx, ...]
        
        # TODO: implement yr
        
        # compute volume flux
        hFacW_x = params['hFacW'][:, :, x_idx].T
        dA_3d = hFacW_x[:, :, np.newaxis]*params['delR'][np.newaxis,:, np.newaxis]*params['delY'][:, np.newaxis, np.newaxis]
        utrans = np.sum(np.sum(uvel_xx*dA_3d, axis=0), axis=0)/1e6
        
        #print(utrans)
        
        # plot
        plt.plot(tvec_yrs, utrans, '-o', color=cols[ii], lw=2, label="%s" %exp_name_str)
        
        
    plt.grid(True)
    plt.xlabel("Years", fontsize=fz)
    plt.ylabel("Zonal transport (Sv)", fontsize=fz)
    plt.legend(loc=0, fontsize=fz-2, ncol=1)
    plt.axhline(0, linestyle='--', color='k')
    #plt.title(title_str, fontsize=fz+2)
    #tick_params(axis='both', which='major', labelsize=fz)
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)
    #plt.ylim(-1, 2)
       
        
def plotRegionalEvol(exp_names, tr, vname='TEMP', xr_km=[], yr_km=[], zr=[], exp_names_alias=[], xpos=0.5, ypos=0.05, fz=14):
    
    
    fig, axes = plt.subplots(1,1, figsize=(12, 8))
    if len(exp_names)==1:
        axes = axes[:, np.newaxis]
        
    #title_str = 'Mean upper %s m %s across %s$<x<$%s km and %s$<x<$%s km' %(
    
    cols = ['cornflowerblue', 'red']
    secsInYear = 3600*24*365
    
    for ii, exp_name in enumerate(exp_names):
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii] 
        else:
            exp_name_str = exp_names_alias[ii] 
        
        
        #load params
        params, _ = hf.get_exp_params(exp_name, vname=vname)
    
        # load variable 
        vble_tr, tvec = hf.readIters(vname=vname, tr=tr, params=params, mask_zeros=True, returnTimeAvg=False)
        tvec_yrs = np.round(tvec/secsInYear)
        
        # get regional subset
        xr = np.array(xr_km)*1e3
        yr = np.array(yr_km)*1e3
        if len(xr_km)==2:
            xi = np.logical_and(params['xx']>=xr[0], params['xx']<xr[1])
        else:
            xi = np.ones(len(params['xx']), dtype=bool)
            
        if len(yr_km)==2:
            yi = np.logical_and(params['yy']>=yr[0], params['yy']<yr[1])
        else:
            yi = np.ones(len(params['yy']), dtype=bool)
         
        zr = np.array(zr)
        if len(zr)==2:
            assert zr[0]<zr[1]
            assert np.all(zr>=0)
            zz2 = np.abs(params['zz'])
            zi = np.logical_and(zz2>=zr[0], zz2<zr[1])
        else:
            zi = np.ones(len(params['zz']), dtype=bool)
            
        vble_sub = vble_tr[xi, :, :, :][:, yi, :, :][:, :, zi, :]
        shp = vble_sub.shape
        vble_submean = vble_sub.reshape((-1, shp[3])).mean(0)
        # TODO: implement yr
        
        # compute volume flux
#         hFacW_x = params['hFacW'][:, :, x_idx].T
#         dA_3d = hFacW_x[:, :, np.newaxis]*params['delR'][np.newaxis,:, np.newaxis]*params['delY'][:, np.newaxis, np.newaxis]
#         utrans = np.sum(np.sum(uvel_xx*dA_3d, axis=0), axis=0)/1e6
        
        # plot
        plt.plot(tvec_yrs, vble_submean-vble_submean[0], '-o', color=cols[ii], lw=2, label="%s" %exp_name_str)
        
        
    plt.grid(True)
    plt.xlabel("Years", fontsize=fz)
    plt.ylabel("Mean %s anomaly (%s)" %(vname, hf.get_units(vname)), fontsize=fz)
    plt.legend(loc=1, fontsize=fz-2, ncol=1)
    plt.axhline(0, linestyle='--', color='k')
    text_str = 'Region subset:\n%s<=X<%s km\n%s<=Y<%s km\n%s<=Z<%s m' %(xr_km[0], xr_km[1], yr_km[0], yr_km[1], zr[0], zr[1])
    plt.text(xpos, ypos, text_str, fontsize=fz-2, horizontalalignment='center', 
             transform = plt.gca().transAxes, backgroundcolor='w')
    #plt.title(title_str, fontsize=fz+2)
    #tick_params(axis='both', which='major', labelsize=fz)
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)        
        

def plot_surf_wmt(exp_names, tr, use_inst_temp=False,  exp_names_alias=[], xlims=[], show_diff=True, fz=14):  
    
    
#     exp_name_ctrl = exp_names[0]
#     exp_name_prtb = exp_names[1]
    
    vname = 'THETA'
    
    rho0 = 1035 # kg/m3 (from Abernathey et al. 2016)
    cp = 3994 # J/kg/K (from Abernathey et al. 2016)
    s_ref = 35 # psu
    
    plt.figure(figsize = (10, 8))
    
    t0, t1 = tr[:]
    
    wmt_list = []
    for ii, exp_name in enumerate(exp_names):
    
        #load control params
        params, forcing = hf.get_exp_params(exp_name)
        surfQ = forcing['surfQ']

        # load temp data
        if use_inst_temp:
            theta_t0 = hf.readIters(vname='THETA_inst', tr=[t0, t0], params=params)
            theta_t1 = hf.readIters(vname='THETA_inst', tr=[t1, t1], params=params)
            inst_str = '_instT'
        else:
            theta_t0 = hf.readIters(vname='THETA', tr=[t0, t0], params=params)
            theta_t1 = hf.readIters(vname='THETA', tr=[t1, t1], params=params)
            inst_str = ''
            
        theta_surf_t0 = theta_t0[:, :, 0].T
        theta_surf_t1 = theta_t1[:, :, 0].T

        # get coordinate data
        dy = params['delY']
        dx = params['delX']
        dz0 = params['delR'][0]
        hFacS0 = params['hFacS'][0,...]
        DZ0 = dz0*hFacS0
        DA = dy[:, np.newaxis]*dx[np.newaxis, :]
        #debug_here()
        #dX, dY = np.meshgrid(xx, yy)
        
        # plot transformation pt-space   
        if len(exp_names_alias)==len(exp_names):
            exp_names_str = exp_names_alias[ii]
        else:
            exp_names_str = exp_names[ii]  

        # Density bins for MOC calculation
        ptlevs = params['layers_bounds']
        Npt = len(ptlevs)-1
        
        # pre-allocate
        wmt_surf = np.zeros((Npt, 2))
        kk = 0
        for theta_surf in [theta_surf_t0, theta_surf_t1]:
            for jj in range(Npt):

                isel = np.logical_and(theta_surf>=ptlevs[jj], theta_surf<ptlevs[jj+1])

                # compute surface water mass transformation
                DA_pt = DA[isel].flatten()
                DZ0_pt = DZ0[isel].flatten()
                surfQ_pt = surfQ[isel].flatten()
                theta0_pt = theta_surf[isel].flatten() 

                DZ0_pt = np.ma.masked_where(DZ0_pt==0, DZ0_pt)
                dpt_qsurf =  -surfQ_pt/(rho0*cp*DZ0_pt) # flip sign since Q is defined as positive upward
                #alpha = sw.alpha(s_ref, theta0_pt, 0, pt=True)
                #drho_dt = dtheta_dt*alpha

                dpt = ptlevs[jj+1]-ptlevs[jj]
                wmt_surf[jj, kk] = np.ma.sum(DA_pt*dpt_qsurf*DZ0_pt)/dpt      

            ptlevs_2 = ptlevs[:-1] + np.diff(ptlevs)/2
            plt.plot(ptlevs_2, wmt_surf[:, kk]/1e6, '-', linewidth=2, label="Year: %s" %tr[kk])
            #wmt_list.append(wmt_surf/1e6)
    
            kk = kk + 1
    
    plt.gca().set_xticks(np.arange(ptlevs[0], ptlevs[-1]+0.01, 1))
    plt.gca().axhline(0, linestyle=':', linewidth=2, color='k')
    plt.grid(which='major', alpha=0.6)
    plt.legend(loc=0, fontsize=fz-3) 
    if len(xlims)==2:
        plt.xlim(xlims)
    plt.ylabel("Tranformation rate (Sv)", fontsize=fz)
    plt.xlabel("Potential temperature ($^{\circ}$C)", fontsize=fz)
    if use_inst_temp:
        plt.title("Surface water-mass transformation rate (INST)", fontsize=fz)
    else:
        plt.title("Surface water-mass transformation rate (%s)" %exp_names_str, fontsize=fz)

    #debug_here()    
    fname = 'surf_wmt%s_%s-%s.pdf' %(inst_str, tr[0], tr[1])
    #plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')         
        