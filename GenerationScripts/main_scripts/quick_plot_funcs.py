import numpy as np
import os
import matplotlib.pylab as plt
import helper_functions as hf
import imp
import xarray as xr
import matplotlib.gridspec as gridspec
import seawater as sw
import pickle


plot_dir = '../plots'
os.makedirs(plot_dir, exist_ok=True)

output_dir = '../output'
os.makedirs(output_dir, exist_ok=True)


def plot_forcing(exp_path, fz=14, xlim1=[-10, 10], xlim2=[-0.15, 0.15]):
    
    # load surface forcing and bathy
    params, forcing = hf.get_exp_params_py(exp_path)
    yy = params['yy']/1e3
    xx = params['xx']/1e3
    tau_xmean = forcing['zonalWind'].mean(axis=1)
    #q_xmean = forcing['surfQ'].mean(axis=1)
    xc = -500e3
    xidx = np.argmin(np.abs(xc-params['xx']))
    q_mean = forcing['surfQ'].mean(axis=1)
    q0 = forcing['surfQ'][:, 0]
    q_aabw = forcing['surfQ'][:, xidx]
    
    fig = plt.figure(figsize=(12, 6))
    
    # plot surface forcing
    col1 = 'red'
    ax1 = plt.subplot(121)
    ax1.plot(q_mean, yy, linewidth=2, color=col1, label='zonal mean')
    ax1.plot(q_aabw, yy, '--', linewidth=2, color=col1, label='bottom water region')
    ax1.plot(q0, yy, ':', linewidth=2, color=col1, label="elsewhere")
    plt.legend(loc=0, fontsize=fz-4)
    ax1.set_xlabel("Ocean heat loss (W/m$^2$)", fontsize=fz, color=col1)
    ax1.tick_params(axis='x', labelcolor=col1)
    ax1.set_xlim(xlim1)
    
    ax1.tick_params(axis='both', which='major', labelsize=fz)
    
    #ax1b = ax1.twiny()
    ax1b = plt.subplot(122)
    col2 = 'blue'
    ax1b.plot(tau_xmean, yy, linewidth=2, color=col2)
    ax1b.set_xlabel("Zonal wind stress (N/m$^2$)", fontsize=fz, color=col2)
    ax1b.tick_params(axis='x', labelcolor=col2)
    ax1b.set_xlim(xlim2)
    ax1b.tick_params(axis='both', which='major', labelsize=fz)
    
    axes = [ax1, ax1b]
    for ii in range(2):
        axes[ii].grid(True)
        axes[ii].set_ylabel("Y (km)", fontsize=fz)
        axes[ii].set_ylim(0, 2500)
        axes[ii].tick_params(axis='both', which='major', labelsize=fz)
        axes[ii].axvline(0, linestyle='--', color='k')
   
    plt.subplots_adjust(wspace=0.4)
    
    
def plot_tau_curl_bathy(exp_path, fz=14, xlim1=[-0.2, 0.2], xlim2=[-6, 6]):
    
    # plot wind stress curl
    params, forcing = hf.get_exp_params_py(exp_path)

    yy = params['yy']/1e3
    xx = params['xx']/1e3
    ZZ = forcing['bathy']
    tau_xmean = forcing['zonalWind'].mean(axis=1)
    curl_tauxmean = -np.gradient(tau_xmean, yy*1e3)
    ZZ[-1:, :] = np.nan
    XX, YY = np.meshgrid(xx, yy)

    #debug_here()
    fig = plt.figure(figsize=(15, 6))
    gs = gridspec.GridSpec(1, 8)
    ax1a = plt.subplot(gs[0, :2])
    
    plt.sca(ax1a)
    col1 = 'blue'
    ax1a.plot(tau_xmean, yy, linewidth=2, color=col1)
    ax1a.set_xlabel("Zonal wind stress (N/m$^2$)", fontsize=fz, color=col1)
    ax1a.set_ylabel('Y (km)', fontsize=fz)
    ax1a.tick_params(axis='x', labelcolor=col1)
    ax1a.set_xlim(xlim1)
    plt.grid(True)
    ax1a.tick_params(axis='both', which='major', labelsize=fz)
    
    
    ax1b = ax1a.twiny()
    col2 = 'orangered'
    plt.plot(curl_tauxmean*1e7, yy, color=col2)
    plt.axvline(0, linestyle='--', color='k')
    plt.ylim(0, 2500)
    
    ax1b.tick_params(axis='x', labelcolor=col2)
    ax1b.set_xlim(xlim2)
    plt.xlabel('wind stress curl ($10^{-7}$ N/m$^3$)', fontsize=fz, color=col2)
    plt.ylabel('Y (km)', fontsize=fz)

    ax2 = plt.subplot(gs[0, 2:])
    im = plt.contourf(XX, YY, ZZ, 30, cmap=plt.cm.inferno, vmin=-4000, vmax=0)
    cb = plt.colorbar(im)
    cb.set_label('Depth (m)', fontsize=fz)
    plt.xlabel('X (km)', fontsize=fz)
    plt.ylabel('Y (km)', fontsize=fz)
    #plt.title("%s bathymetry" %exp_name_str, fontsize=fz)
    ax2.tick_params(axis='both', which='major', labelsize=fz)
    plt.ylim(0, 2500)

    plt.subplots_adjust(wspace=2)
    
    
    
def plot_mean_series(exp_root_dir, exp_names_list, vname='THETA_inst',
                     exp_names_alias=[], zr=[], xr=[], yr=[], tstep=1, fz=14, 
                     load_previous=True, save_plots=False, time_units='years', xlim=[], 
                     xstep=1, ylim=[], linestyle='-'):
    
    """
    This should supercede  plot_mean_temp_series()
    
    Function to compute the domain mean time series of a variable. This is useful for checking equilibration.
    
    Input:
    exp_root_dir: path to experiments directory
    exp_names_list: list of experiment names 
    vname: name of variable to compute mean [THETA_inst]
    exp_names_alias: list of aliases for the experiment names []
    xr, yr, zr: lists defining min/max limits for control volume. Full ranges are used if list is empty. [], [], []
    tstep: int that specifies how subsample output. 1 to take every output, 2 to take every other output, etc. [1]
    fz: fontsize for plots [14]
    load_previous: bool specifying whether to load previous calculation [True]
    save_plots:bool specifying whether to save plots [False]
    time_units: str specifying xaxis units (secs, minutes, days, months, [years])
    xlim: list specifying min/max of xaxis  []
    xstep: year increments for xticks []
    ylim: list specifying min/max of yaxis  []
    
    """
    
    secsInMins = 60
    secsInhrs= 60*secsInMins
    secsInDays = 24*secsInhrs
    secsInMons = 30*secsInDays
    secsInYrs = 365*secsInDays
    
    units_conv = {'minutes':secsInMins, 'hours':secsInhrs, 'days':secsInDays, 
                  'months':secsInMons, 'years':secsInYrs}

    units_list = list(units_conv.keys())
    
    plt.figure(figsize=(10, 6))
    
    print("computing domain-average %s. May take a while..." %vname)
    cols = ['cornflowerblue', 'orange', 'green']
    xmax = []
    show_ctrl = False
    
    for ii, exp_name in enumerate(exp_names_list):
        
        exp_path_ii = os.path.join(exp_root_dir, exp_name)
        
        if ii==0:
            params, _ = hf.get_exp_params_py(exp_path_ii, vname=vname)
        
        pdata = hf.computeDomainInt(exp_path_ii, vname=vname,zr=zr, xr=xr, yr=yr, tstep=tstep, 
                                    load_previous=load_previous)

        if len(exp_names_alias)==len(exp_names_list):
            exp_name_str = exp_names_alias[ii]
        else:
            exp_name_str = exp_name
            
            
        tvec = pdata['tt']
        for jj in range(len(units_list)):
            tvec = np.array(pdata['tt'])/units_conv[units_list[jj]]
            if np.max(tvec)<1000:
                break
            
            #print(np.max(tvec))
            
        

        # plot time series
        plt.plot(tvec, pdata['vble_avg'], linestyle, label=exp_name_str, color=cols[ii], linewidth=3)
        plt.legend(loc=0, fontsize=fz-2)
        plt.ylabel("%s (%s)" %(vname, hf.get_units(vname)), fontsize=fz)
        plt.xlabel(units_list[jj], fontsize=fz)

        # set label size
        plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        

            
    if len(ylim)==2:
        plt.ylim(ylim)
    
    # plt.grid(which='both')
    plt.grid(which='minor', alpha=0.4)
    plt.grid(which='major', alpha=0.6)

    def get_subset_str(sr, dim='z'):

        if len(sr)==2:
            if dim=='z':
                rnge = [np.abs(sr).min(), np.abs(sr).max()]
                s_str = "%s=%s -- %s m" %(dim, *rnge)
            else:
                rnge = sr.min()/1e3, sr.max()/1e3
                s_str = "%s=$%s$ -- $%s$ km" %(dim, *rnge)

        else:
            s_str = s_str = "%s = --" %(dim)
#             dd = dim*2 # duplicates str
#             if dim=='z':
#                 rnge = np.round([np.abs(params[dd]).min()/1e3, np.abs(params[dd]).max()/1e3])
#                 s_str = "%s=[%s$ -- $%s] m" %(dim, *rnge)
#             else:
#                 rnge = np.round([params[dd].min()/1e3, params[dd].max()/1e3])
#                 s_str = "%s=[$%s$ -- $%s$]km" %(dim,*rnge)
                

        return s_str


    z_str = get_subset_str(zr, dim='z')
    y_str = get_subset_str(yr, dim='y')
    x_str = get_subset_str(xr, dim='x')
    sub_str = "_".join([z_str, y_str, x_str])
                
    if len(sub_str)==0:
        title_str = 'Domain mean %s' %vname
    else:
        combo_str = ", ".join([z_str, y_str, x_str])
        title_str = 'mean %s\n(%s)'%(vname, combo_str)
        
        
    plt.title("%s" %title_str, fontsize=14)
    
    
    if save_plots:
        fname = 'mean_%s_tseries%s.pdf' %(vname, sub_str)
        plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')

        
        
def plot_layer_avg_anom(exp_path, vname, t0, zr, exp_name_alias='', clvls=[], fz=14, 
                        cmap=plt.cm.viridis, lcol='w', bg_col='k', add_clines=True, 
                        save_plots=False, cls=4):
    
    """
    function to plot and compare layer-averaged quantities at two different time instances
    
    Input:
    
    
    """
    
    secsInYear = 3600*24*365

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    title_str = 'Mean %s' %vname
    
    exp_name = exp_path.split('/')[-1]
    
    # load output
    params, _ = hf.get_exp_params_py(exp_path, vname=vname)
    
    tyears = params['dumpIters']*params['deltaT']/secsInYear
    
    if t0<0:
        tr = [tyears[t0], tyears[t0]]
    else:
        tr = [t0, t0]

    vble_t0 = hf.readIters(vname=vname, tr=tr, params=params)
    vble_t0_zavg = hf.get_slice_mean(vble_t0, zr, vname, params, dim='z')
    
    #print(vble_t0.mean())

    title_str = 'Mean %s-%sm %s' %(zr[0], zr[1], vname)


    YY,XX = np.meshgrid(params['yy']/1000, params['xx']/1000)

    if len(clvls)==0:
        clvls = np.linspace(-0.9*vble_t0_zavg.max(), 0.9*vble_t0_zavg.max(), 25)

    if len(exp_name_alias)==0:
        exp_name_str = exp_name 
    else:
        exp_name_str = exp_name_alias

    clines = clvls[::cls]
    
    tr_rnd = np.round(tr, decimals=2)

    plt.sca(ax)
    im = plt.contourf(XX, YY, vble_t0_zavg, clvls, cmap=cmap, extend='both')
    cb = plt.colorbar(im, extend='both', ticks=clines)
    cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
    cb.ax.tick_params(labelsize=fz) 
    plt.title("%s\n%s at year %s" %(exp_name_str, title_str, tr_rnd[0]), fontsize=fz+2)
    plt.text(0.95*XX.min(), 0.9*YY.max(), "(A)", fontsize=fz+4, color='k')

    if add_clines:
        cs = ax.contour(XX, YY, vble_t0_zavg, clines, colors=lcol, linewidths=1)
        plt.gca().clabel(cs, inline=1, fontsize=fz-2, fmt='%.1f')

    plt.xlabel("X (km)", fontsize=fz)
    plt.ylabel("Y (km)", fontsize=fz)
    
    # set label size
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)
    
    # change BG color
    ax.set_facecolor(bg_col)
        
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

def plot_zonal_sect_anom(exp_path, vname, t0, xr, exp_name_alias='', clvls=[], add_clines=True, 
                         fz=14, cmap=plt.cm.viridis, plot_MLD=False):
    
    """
    function to plot zonally-averaged quantities
    """
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    #load control params
    params, forcing = hf.get_exp_params_py(exp_path, vname=vname)
   
    tr = [t0, t0]
    
    # load output
    if vname.lower() == 'ssh':
        raise ValueError("%s not supported for section plots." %vname)
    else:
        vble = hf.readIters(vname=vname, tr=tr, params=params)

    vble_xavg = hf.get_zonal_mean(vble, xr, vname, params)
    
    
#     # compute MLD
#     if plot_MLD:
        
#         print("Warning: MLD function uses a temperature threshold. Need to update!")
        
#         if vname.endswith('_inst'):
#             suffix = '_inst'
#         else:
#             suffix = ''
            
#         print(suffix)
        
#         temp = hf.readIters(vname='THETA'+suffix, tr=tr, params=params)
#         mld_xavg = hf.compute_MLD_fast(temp, params, xr=xr, dT=0.2, get_xavg=True)
        
        
    #fig, axes = plt.subplots(1, 1, figsize=(12, 15))

    ZZ,YY = np.meshgrid(params['zz'], params['yy']/1000)
    xr_mean = np.array(xr).mean()/1e3
    
    #debug_here()
    
    if len(clvls)==0:
        clvls = np.linspace(-0.9*vble_xavg.max(), 0.9*vble_xavg.max(), 25)
        
    if len(exp_name_alias)>0:
        exp_name_str = exp_name_alias
    else:
        exp_name_str = exp_name
        
        
        
    plt.sca(ax)
    im = plt.contourf(YY, ZZ, vble_xavg, clvls, cmap=cmap, extend='both')
    
    #if vname in ['THETA', 'THETA_inst']:
    if add_clines:
        clines = clvls[::2]
        cs = plt.contour(YY, ZZ, vble_xavg, clines, colors='k', linewidths=1)
        plt.gca().clabel(cs, inline=1, fontsize=fz-2, fmt='%.1f')
        plt.gca().axvline(1000, linestyle='--', linewidth=1.5, color='orange')
     
    if plot_MLD:
        plt.plot(params['yy']/1000, mld_xavg, lw=2, color='red')
        
     
    tr_rnd = np.round(tr, decimals=2)
        
    cb = plt.colorbar(im, extend='both', ticks=clines[::2])
    cb.set_label('%s (%s)' %(vname, hf.get_units(vname)), fontsize=fz)
    cb.ax.tick_params(labelsize=fz) 
    plt.title("Mean %s between $%i<X<%i$ km for years %s-%s\n%s" 
              %(vname, xr[0]/1e3, xr[-1]/1e3, tr_rnd[0], tr_rnd[-1], exp_name_str), fontsize=fz+2)
    

    plt.sca(ax)
    # add labels
    plt.xlabel("Y (km)", fontsize=fz)
    plt.ylabel("Depth (m)", fontsize=fz)

    # set label size
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)

    # change BG color
    ax.set_facecolor('k') 
    

