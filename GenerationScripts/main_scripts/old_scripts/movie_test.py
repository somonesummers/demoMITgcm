import numpy as np
import matplotlib.pylab as plt
import imp
import xarray as xr
import os


import quick_plot_funcs as qpf
import overturning_funcs as of
import helper_functions as hf
import plot_functions as pf
import plot_functions_prtb as pfp

import imp

imp.reload(qpf)
imp.reload(of)
imp.reload(hf)
imp.reload(pf)
imp.reload(pfp)

import xgcm
from xmitgcm import open_mdsdataset
from xmovie import Movie

exp_names = ['gyre_run_SS_WAP_DP_AABW_mode2_tauW_v10_1yr', 'gyre_run_SS_WAP_DP_NWRx2h_AABW_mode2_tauW_v10_1yr']
exp_names_alias = ['SS+DP', 'SS+DP+NWRx2h']
exp_dir = '/central/groups/AndyThompsonGroup/earlew/MITgcm_PG/experiments/' 
fz = 14  # fontsize


params1,_ = hf.get_exp_params(exp_names[0])
params2,_ = hf.get_exp_params(exp_names[1])

ds1 = open_mdsdataset(params1['results_path'], prefix=['UVEL_inst', 'VVEL_inst'], geometry='cartesian')
ds2 = open_mdsdataset(params2['results_path'], prefix=['UVEL_inst', 'VVEL_inst'], geometry='cartesian')
ds = [ds1, ds2]

grid1 = xgcm.Grid(ds1, periodic=['X'], boundary='extend')
grid2 = xgcm.Grid(ds1, periodic=['X'], boundary='extend')

grid = [grid1, grid2]

u_transport = []
v_transport = []

for dd in [ds1, ds2]:
    u_transport.append(dd.UVEL * dd.dyG * dd.hFacW * dd.drF)
    v_transport.append(dd.VVEL * dd.dxG * dd.hFacS * dd.drF)

div_uv = []
for ii in range(len(ds)):
    xx = (grid[ii].diff(u_transport[ii], 'X') + grid[ii].diff(v_transport[ii], 'Y')) / ds[ii].rA
    div_uv.append(xx)

zeta = []
for ii in range(2):
    zeta.append((-grid[ii].diff(ds[ii].UVEL * ds[ii].dxC, 'Y') + grid[ii].diff(ds[ii].VVEL * ds[ii].dyC, 'X'))/ds[ii].rAz)

zeta_bt = []
for ii in range(len(zeta)):
    zeta_bt.append((zeta[ii] * ds[ii].drF).sum(dim='Z'))
    
    


secsInDay = 3600*24


xx_km = ds[0]['XC']/1e3
yy_km = ds[0]['YC']/1e3
fz = 14

def save_zeta(outdir='zeta_test2', testing=True, clvls=np.arange(-3, 3.1, 0.2), cmap=plt.cm.RdBu_r, dpi=100):
    movie_dir = 'movies/%s/' %outdir
    os.makedirs(movie_dir, exist_ok=True)
    
    tvec = (zeta[0]['time'] - zeta[0]['time'][0])*params1['deltaT']/secsInDay 
    
    vmin = clvls.min()
    vmax = clvls.max()
    
    for tt in range(len(tvec)): #range(len(tvec)):
        fig, axes = plt.subplots(2, 1, figsize=(1200/dpi, 800/dpi), dpi=dpi)
        if tt%5==0:
            print("t = %s/%s" %(tt, len(tvec)))
                  
        for ii in range(2):
            plt.sca(axes[ii])
            zeta_ma = 1e5*np.ma.masked_where(zeta[ii][tt, 0].values==0, zeta[ii][tt, 0].values)
            #im = plt.pcolormesh(xx_km, yy_km, zeta_ma, vmin=-4, vmax=4, cmap=cmap)
            im = plt.contourf(xx_km, yy_km, zeta_ma, clvls, vmin=vmin, vmax=vmax, cmap=cmap)
            plt.xlabel("X (km)", fontsize=fz)
            plt.ylabel("Y (km)", fontsize=fz)
            plt.gca().set_facecolor('k') 
            plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        
        plt.suptitle("surface vorticity (day %i)" %(tvec[tt]+1), fontsize=fz+4, y=0.92)
        plt.subplots_adjust(hspace=0.3)
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        cb = plt.colorbar(im, cax=cbar_ax)
        cb.set_label('$\zeta$ $10^{-5}$ (ms)$^{-1}$', fontsize=fz)
        plt.sca(axes[ii])
        
        fpath = os.path.join(movie_dir, 'frame_%05d.png' %tt)
        plt.savefig(fpath, dpi=dpi)
        
        plt.close('all')
        
        if testing and tt>=5:
            break
            
            
def save_spd(testing=True, clvls=np.arange(0, 0.6, 0.02), dpi=100, cmap=plt.cm.inferno):
    
    movie_dir = 'movies/spd_test1/'
    os.makedirs(movie_dir, exist_ok=True)
    
    spd = []
    for dd in [ds1, ds2]:
        uvel_s = dd.UVEL[:,:3,:,:].mean(dim='Z').values
        vvel_s = dd.VVEL[:,:3,:,:].mean(dim='Z').values
    
        spd.append(np.sqrt(uvel_s**2 + vvel_s**2))
    
    #print(spd[0].shape)
    YY,XX = np.meshgrid(yy_km, xx_km)
    
    
    tvec = (ds[0]['time'] - ds[0]['time'][0])*params1['deltaT']/secsInDay 
    vmin = clvls.min()
    vmax = clvls.max()
    
    for tt in range(len(tvec)): #range(len(tvec)):
        fig, axes = plt.subplots(2, 1, figsize=(1200/dpi, 800/dpi), dpi=dpi)
        if tt%5==0:
            print("t = %s/%s" %(tt, len(tvec)))
                  
        for ii in range(2):
            plt.sca(axes[ii])
            vble_ma = np.ma.masked_where(spd[ii][tt, :, :]==0, spd[ii][tt, :, :])
            im = plt.contourf(xx_km, yy_km, vble_ma, clvls, vmin=vmin, vmax=vmax, cmap=cmap, extend='max')
            plt.xlabel("X (km)", fontsize=fz)
            plt.ylabel("Y (km)", fontsize=fz)
            plt.gca().set_facecolor('0.6') 
            plt.gca().tick_params(axis='both', which='major', labelsize=fz)
        
        plt.suptitle("Near-surface speed (day %i)" %(tvec[tt]+1), fontsize=fz+4, y=0.92)
        plt.subplots_adjust(hspace=0.3)
        
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
        cb = plt.colorbar(im, cax=cbar_ax, extend='max')
        cb.set_label('Speed (m/s)', fontsize=fz)
        plt.sca(axes[ii])
        
        fpath = os.path.join(movie_dir, 'frame_%05d.png' %tt)
        plt.savefig(fpath, dpi=dpi)
        
        plt.close('all')
        
        if testing and tt>=5:
            break
                 
            
#x=ds[0]['XC'].values/1e3, y=ds[0]['YC'].values/1e3,     
#mov = Movie(zeta[0][:,0,:,:], vmin=-5e-5, vmax=5e-5, cmap=plt.cm.bwr)

#mov.save('movies/movie.mp4')    
