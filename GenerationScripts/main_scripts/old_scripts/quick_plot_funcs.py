import numpy as np
import os
import matplotlib.pylab as plt
import helper_functions as hf
import imp
import xarray as xr
import matplotlib.gridspec as gridspec
import seawater as sw
import pickle

from IPython.core.debugger import set_trace
debug_here = set_trace

imp.reload(hf)

# plot_dir = 'plots_temporary/'
# output_dir = 'output/'

# TODO: import these directories



def plot_bathy(exp_names, exp_names_alias=[], fz=14):
    
#     analysis_dir = '/scratch/users/earlew/research/modeling/MITgcm_PG/analysis/'
#     output_dir = os.path.join(analysis_dir, 'output/')
#     plot_dir = os.path.join(analysis_dir, 'plots/')
    
    
    nexp = len(exp_names)
    if nexp <=3:
        fig, axes = plt.subplots(nexp, 1, figsize=(8, 10))
    else:
        fig, axes = plt.subplots(2, nexp//3, figsize=(12, 15))
        
    # ---------- plot bathy ---------- #
    for ii, exp_name in enumerate(exp_names):
        params, forcing = hf.get_exp_params_py(exp_name)
        yy = params['yy']/1e3
        xx = params['xx']/1e3
        ZZ = forcing['bathy']
        ZZ[:, -2:] = -4000

        #ax2 = plt.subplot(gs[0, 4:])
        #ax_list.append(ax2)
        plt.sca(axes[ii])
        ZZ[-1:, :] = np.nan
        XX, YY = np.meshgrid(xx, yy)
        
        if len(exp_names_alias)==0:
            exp_name_str = exp_names[ii]
        else:
            exp_name_str = exp_names_alias[ii]

        im = plt.contourf(XX, YY, ZZ, 30, cmap=plt.cm.inferno, vmin=-4000, vmax=0)
        if ii==nexp-1:
            plt.xlabel('X (km)', fontsize=fz)
        plt.title("%s" %exp_name_str, fontsize=fz)
        #axes[ii].set_yticklabels([])
        axes[ii].tick_params(axis='both', which='major', labelsize=fz)
        plt.ylim(0, 2500)

    plt.subplots_adjust(hspace=0.4)
    
    fig.subplots_adjust(right=0.9)
    # put colorbar at desire position
    cbar_ax = fig.add_axes([0.93, 0.15, 0.04, 0.7])
    #cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.02])
    cb = fig.colorbar(im, cax=cbar_ax, ticks=np.arange(-4000, 0.1, 500))
    cb.set_label('Depth (m)', fontsize=fz)
      
    
    

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
    
    #debug_here()
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
    
def plot_layer_avg_anom(exp_name, vname, t0, zr, exp_name_alias='', clvls=[], fz=14, 
                        cmap=plt.cm.viridis, lcol='w', bg_col='k', add_clines=True, 
                        save_plots=False, cls=4):
    
    """
    function to plot and compare layer-averaged quantities at two different time instances
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))

    title_str = 'Mean %s' %vname
    
#     if t0<0:
#         # grab most recent 

    tr = [t0, t0]
    
    # load output
    if vname.lower() == 'ssh':

        # load variable at t=t0
        
        #load params
        params, _ = hf.get_exp_params_py(exp_name, vname='PHIHYD')
        
        vble_t0 = hf.readIters(vname='PHIHYD', tr=tr, params=params)
        vble_t0 = vble_t0/params['gravity']
        theta = hf.readIters(vname='THETA', tr=tr, params=params)
        vble_t0 = np.ma.masked_where(theta.mask, vble_t0)
        vble_t0_zavg = vble_t0[...,0] # not really a zavg. Just trying to keep the names consistent
        
    elif vname.startswith('SPD'):
        
        if vname.endswith('inst'):
            ext = '_inst'
        else:
            ext = ''
                
        cbar_ext = 'max' 
        
        params, _ = hf.get_exp_params_py(exp_name, vname='UVEL'+ext)
        
        vble_t0_a = hf.readIters(vname='UVEL'+ext, tr=tr, params=params, printStatus=False)
        vble_t0_b = hf.readIters(vname='VVEL'+ext, tr=tr, params=params, printStatus=False)
        vble_t0 = np.ma.sqrt(vble_t0_a**2 + vble_t0_b**2)
        vble_t0_zavg = hf.get_layer_mean(vble_t0, zr, vname, params)
        
        title_str = 'Mean %s-%sm %s' %(zr[0], zr[1], vname)

    elif vname == 'BT_psi':
        
        params, _ = hf.get_exp_params_py(exp_name, vname='THETA')
        
        theta = hf.readIters(vname='THETA', tr=tr, params=params)
        vble_t0_zavg = hf.calc_BT_streamfunction(exp_name, tr)
        vble_t0_zavg = np.ma.masked_where(theta[:, :, 0].mask, vble_t0_zavg[1:, 1:]/1e6)

    else:

        params, _ = hf.get_exp_params_py(exp_name, vname=vname)
        
        vble_t0 = hf.readIters(vname=vname, tr=tr, params=params)
        vble_t0_zavg = hf.get_layer_mean(vble_t0, zr, vname, params)

        title_str = 'Mean %s-%sm %s' %(zr[0], zr[1], vname)


    #vble_zavg_ano = vble_t1_zavg - vble_t0_zavg

    # plot vble at t0, t1, and their difference

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

def plot_zonal_sect_anom(exp_name, vname, t0, xr, exp_name_alias='', clvls=[], add_clines=True, 
                         fz=14, cmap=plt.cm.viridis, plot_MLD=False):
    
    """
    function to plot zonally-averaged quantities
    """
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    #load control params
    params, forcing = hf.get_exp_params_py(exp_name, vname=vname)
   
    tr = [t0, t0]
    
    # load output
    if vname.lower() == 'ssh':
        raise ValueError("%s not supported for section plots." %vname)
    else:
        vble = hf.readIters(vname=vname, tr=tr, params=params)

    vble_xavg = hf.get_zonal_mean(vble, xr, vname, params)
    
    
    # compute MLD
    if plot_MLD:
        
        print("Warning: MLD function uses a temperature threshold. Need to update!")
        
        if vname.endswith('_inst'):
            suffix = '_inst'
        else:
            suffix = ''
            
        print(suffix)
        
        temp = hf.readIters(vname='THETA'+suffix, tr=tr, params=params)
        mld_xavg = hf.compute_MLD_fast(temp, params, xr=xr, dT=0.2, get_xavg=True)
        
        
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
        
    #clines = [0, 1, 2, 3, 4, 5]
    
    
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
    
    #plt.text(0.93*YY.max(), 0.1*ZZ.min(), "(A)", fontsize=fz+4, color='k')
    #plt.text(0.02*YY.max(), 0.95*ZZ.min(), "(A)", fontsize=fz+4, color='w')


#     axes[1].set_title("Mean %s between $%i<X<%i$ km for years %s-%s\n%s" 
#                   %(vname, xr[0]/1e3, xr[-1]/1e3, tr[0], tr[-1], exp_names_str[1]), fontsize=fz+2)
#     axes[2].set_title("B - A", fontsize=fz+2)
#     plt.subplots_adjust(hspace=0.4)
    plt.sca(ax)
    # add labels
    plt.xlabel("Y (km)", fontsize=fz)
    plt.ylabel("Depth (m)", fontsize=fz)

    # set label size
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)

    # change BG color
    ax.set_facecolor('k') 
    
    
def plot_streamfunction(exp_name, t0, zr, exp_name_alias='', fz=14, clvls=[], 
                           clvls_ano=[], cmap=plt.cm.viridis, save_plots=False):

    #load control params
    #exp_name_ctrl = exp_names[0]
    params, _ = hf.get_exp_params_py(exp_name)
    
    tr = [t0, t0]
    
    #load perturbation params 
#     exp_name_prtb = exp_names[1]
#     params_prtb, _ = hf.get_exp_params_py(exp_name_prtb)
    
    # compute barotropic streamfunction
    psi = hf.calc_streamfunction(exp_name, tr, zr)
    #psi_prtb = hf.calc_BT_streamfunction(exp_name_prtb, tr)
    
    theta = hf.readIters(vname='THETA', tr=tr, params=params, printStatus=False)
    psi = np.ma.masked_where(theta[:, :, 0].mask, psi[1:, 1:]/1e6)
    
    # plot ctrl, perturb and difference
    #fig, axes = plt.subplots(3, 1, figsize=(12, 15))
    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    
    YY, XX = np.meshgrid(params['yy']/1000, params['xx']/1000)
    
    if len(clvls)==0:
        clvls = np.linspace(-0.9*psi.max(), 0.9*psi.max(), 25)
    
    if len(clvls_ano)==0:
        clvls_ano = np.linspace(-0.9*psi.max(), 0.9*psi.max(), 25)
        
    if len(exp_name_alias)>0:
        exp_name_str = exp_name_alias
    else:
        exp_name_str = exp_name

    zr = np.array(zr)
    if np.abs(zr).max()>=4000:    
        title_str = 'Mean barotropic streamfunction'
    else:
        title_str = 'Mean %s-%s m  $\Psi$' %(zr[0], zr[1])
    #clines = np.arange(-2, 13, 1)
    clines = clvls[::2]
    lcol = 'k'
    
    plt.sca(ax)
    im = plt.contourf(XX, YY, psi, clvls, cmap=cmap, extend='both')
    cs = plt.contour(XX, YY, psi, clines, colors=lcol, linewidths=0.5)
    plt.gca().clabel(cs, inline=1, fontsize=fz-2, fmt='%.1f')
    cb = plt.colorbar(im, extend='both', ticks=clines[::2])
    cb.set_label('$\Psi$ ($10^6$ m$^3$/s)', fontsize=fz)
    plt.title("%s for years %s-%s\n%s" %(title_str, tr[0], tr[-1], exp_name_str), fontsize=fz+2)
    #plt.text(0.95*XX.min(), 0.9*YY.max(), "(A)", fontsize=fz+4, color='k')
    

#     plt.sca(axes[1])
#     im = plt.contourf(XX, YY, psi_prtb, clvls, cmap=cmap, extend='both')
#     plt.contour(XX, YY, psi_prtb, clines, colors=lcol, linewidths=0.5)
#     cb = plt.colorbar(im, extend='both', ticks=clines[::2])
#     cb.set_label('$\Psi_{BT}$ (m$2$/s)', fontsize=fz)
#     plt.text(0.95*XX.min(), 0.9*YY.max(), "(B)", fontsize=fz+4, color='k')

#     plt.sca(axes[2])
#     im = plt.contourf(XX, YY, psi_prtb-psi_ctrl, clvls_ano, cmap=plt.cm.bwr, extend='both')
#     cb = plt.colorbar(im, extend='both')
#     cb.set_label('$\Psi_{BT}$ (m$2$/s)', fontsize=fz)
#     plt.text(0.95*XX.min(), 0.9*YY.max(), "(C)", fontsize=16, color='k')

    
#     axes[1].set_title("%s for years %s-%s\n%s" %(title_str1, tr[0], tr[-1], exp_names_str[1]), fontsize=fz+2)
#     axes[2].set_title("B - A", fontsize=fz+2)
    plt.xlabel("X (km)", fontsize=fz)
    plt.ylabel("Y (km)", fontsize=fz)
    plt.subplots_adjust(hspace=0.4)
    # set label size
    plt.gca().tick_params(axis='both', which='major', labelsize=fz)
    cb.ax.tick_params(labelsize=fz) 

    # change BG color
    ax.set_facecolor('k')
    
#     if save_plots:
#         fname = 'BT_Psi_yrs%s-%s.pdf' %(tr[0], tr[1])
#         plt.savefig(os.path.join(plot_dir, fname), bbox_inches='tight')
        