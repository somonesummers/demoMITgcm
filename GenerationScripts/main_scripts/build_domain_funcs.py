import numpy as np
import matplotlib.pylab as plt


def add_bathy(run_config, domain_params, grid_params, fz=12, show_plot=True, 
              y_slices_m=[], x_slices_m=[], cmap='cividis'): # y_slices: 850, 1500
    
    # TODO: 
    # + add doc string
    # + make variable names more consistent
    
    # zonal grid
    dx = domain_params['Lx']/grid_params['Nx'] # x-grid res
    xx = np.arange(grid_params['Nx'] )*dx
    xx = xx-np.mean(xx) # offset so that x=0 is in the center of the domain

    # meridional grid   
    dy = domain_params['Ly']/grid_params['Ny']   # y-grid res
    yy = np.cumsum(dy*np.ones(grid_params['Ny'] )) - dy/2 
    # Y increases with increasing latitude. Also, following Andrew and offseting Y by dy/2)

    # Plotting mesh
    YY, XX = np.meshgrid(yy, xx)
    
    # create flat-bottomed ocean
    H = domain_params['H']
    h = -H
    
    # add southern shelf
    print("Adding southern shelf...")
    if domain_params['add_SS']:
        h += domain_params['SS_H']*(1-np.tanh((YY-domain_params['SS_slope_Y'])/
                                             domain_params['SS_slope_hW']))/2
    
    if domain_params['add_SC']:
        print("Adding Shelf Canyon...")
        # create the canyon
        cc = domain_params['SC_slope_cc']
        X0 = domain_params['SC_X0']
        dX = domain_params['SC_dX']
        H_MR = domain_params['SC_MR']
        Y0 = domain_params['SC_Y0']
        Yend = domain_params['SC_Yend']
        ss = domain_params['SC_ss']
        
#        SC_cut = add_meridional_channel(XX, YY, x0, dX, cc, H_sill, H, MR_y=None)
#        h_SC, SC_y = add_meridional_canyon_old(XX, YY, X0, dX, cc, H_MR, Y0, dY)
        h_MR, SC_y = add_meridional_canyon(XX, YY, X0, dX, cc, H_MR, h, Y0, Yend, ss)
                
        # add to bathymetry
        h += h_MR
   
    if domain_params['add_DP']:
        print("Adding meridional continent...")
        # create meridional continent that spans the domain
        dd = domain_params['slope_hW']/domain_params['DP_dX'] # ratio of slope length-scale to continental width 
        cc = domain_params['slope_hW']
        X0 = domain_params['DP_X']
        dX = domain_params['DP_dX']
        
        try: 
            MR_HoS = domain_params['MR_HoS']
        except KeyError:
            MR_HoS = 0
        
        H_MR = H+MR_HoS # max height of meridional ridge
        
        h_MR, MR_x = add_meridional_ridge(XX, YY, X0, dX, cc, H_MR)

        
        #cut out drake passage 
        print("cutting out Drake Passage...")
        Y0 = domain_params['DP_Y_mid']
        dY = domain_params['DP_dY']
        H_sill = domain_params['DP_sill_H']
        
        DP_cut = add_zonal_channel(XX, YY, Y0, dY, cc, H_sill, H_MR, MR_x=MR_x)
           
        # add to bathymetry
        h += h_MR +DP_cut
        
        
        
    if domain_params['add_ZR']:
        print("adding zonal ridge...")
        Y0 = domain_params['ZR_Y_mid']
        dY = domain_params['ZR_dY']
        cc = domain_params['slope_hW']
        ZR_H = domain_params['ZR_H']
        
        X0 = domain_params['ZR_X_min']+domain_params['ZR_dX']/2 
        dX = domain_params['ZR_dX']
        
        h_ZR, ZR_y = add_zonal_ridge(XX, YY, Y0, dY, cc, ZR_H, X0=X0, dX=dX)
        
        # add to bathymetry
        h += h_ZR
        
        
        if domain_params['add_SAC']:
            print("adding zonal ridge...")
            for cut_X0 in domain_params['SAC_X0']:
                SA_cut = add_meridional_channel(XX, YY, X0=cut_X0, 
                                                dX=domain_params['SAC_dX'], cc=cc, 
                                                H_sill=domain_params['SAC_sill_H'], 
                                                H=ZR_H, MR_y=ZR_y)
                
                h += SA_cut
                
                
    ## Enforce minimum ocean depth and add walls along northern and southern boundaries
    h[h>-100] = 0 # EW NOTE: Not entirely why this done
    h[:, 0] = 0 
    h[:, -1] = 0
                
        
    #print(XX.shape)
    if show_plot:
       
        
        # plot plan view
        plt.figure(figsize=(10, 5))
        pp = plt.pcolormesh(XX/1e3, YY/1e3, h, vmin=h.min(), vmax=0, cmap=cmap)
        plt.contour(XX/1e3, YY/1e3, h, [-5], colors='k')
        plt.axis('equal') # added by ED
        ax0 = plt.gca()
        plt.ylabel('Y (m)', fontsize=fz)
        plt.xlabel('X (m)', fontsize=fz)
        
        cb= plt.colorbar(pp)
        cb.set_label('Depth', fontsize=fz)
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        
        # plot meridonal slices
        plt.sca(axes[0])
        
        for x_slice in x_slices_m:
            x_idx = np.argmin(np.abs(xx-x_slice*1e3))
            plt.plot(yy/1e3, h[x_idx, :], '.-', label="X=%s m" %(xx[x_idx]/1000))
            ax0.axvline(xx[x_idx]/1000, linestyle='--', color='w')
        
        
        plt.xlabel('Y (m)', fontsize=fz)
        plt.ylabel('H (m)', fontsize=fz)
        plt.title('Meridional slices', fontsize=fz)
        
        major_ticks_x=np.linspace(0,600,11)
        major_ticks_z=np.linspace(-4000,0,9)
        minor_ticks=np.linspace(0,600,21)

        axes[0].set_xticks(major_ticks_x)
        axes[0].set_yticks(major_ticks_z)
        axes[0].set_xticks(minor_ticks,minor=True)
        axes[0].grid(which="major",alpha=0.6)
        axes[0].grid(which="minor",alpha=0.3)

        #plt.grid()
        #plt.legend()
        
        # plot zonal slices
        plt.sca(axes[1])
        
        for y_slice in y_slices_m:
            y_idx = np.argmin(np.abs(yy-y_slice*1e3))
            plt.plot(xx/1e3, h[:, y_idx], '.-', label="Y=%s m" %(yy[y_idx]/1000))
            ax0.axhline(yy[y_idx]/1000, linestyle='--', color='w')
        
        plt.xlabel('X (m)', fontsize=fz)
        plt.ylabel('H (m)', fontsize=fz)
        plt.title('Zonal slices', fontsize=fz)
        plt.grid(True)
        plt.legend()

        major_ticks_x=np.linspace(-500,500,11)
        major_ticks_z=np.linspace(-4000,0,9)
        minor_ticks=np.linspace(-500,500,21)
        
        plt.subplots_adjust(wspace=0.3)
        
    return h
        
    
def add_meridional_ridge(XX, YY, X0, dX, cc, H_MR, Y0=np.inf, dY=np.inf):
    
    """
    function creates an elongated ridge, centered at X0, with longitudinal extent of dX and max height of H_MR.
    The steepness of the slopes (in the x-direction) is controlled by cc. The latitudinal extent is infinite 
    by default, but can be constrained by Y0 and dY.
   
    """
    
    dd = cc/dX
    # dd: measures the prominence of the plateau. If the plateau is narrow, use a gaussian function.
    # otherwise, use tanh function. This is necessary since the tanh functions (as implemented here) can 
    # grossly underestimate the height of narrow features
    
    if dd>0.2:
        print("using gaussian funcs...")
        # create ridge using gaussian function spanning the length of the domain
        MR_x = np.exp(-(XX-X0)**2/(0.5*dX**2)) 
        h_MR = H_MR*MR_x
        
    else:
        print("using tanh funcs...")
        # combine two tanh functions to get meridional ridge
        MR_x = (1 - np.tanh((XX-X0-dX/2)/cc))/2 + (1 + np.tanh((XX-X0+dX/2)/cc))/2
        MR_x -= 1 # need to offset by H to get the right height
        h_MR = H_MR*MR_x
               
    # if ridge is finite, limit meridional extent
    if YY.min()<Y0<YY.max() and dY<(YY.max()-YY.min()):
        print("setting y bounds...")
        # wrap in meridional function 
        dd2 = cc/dY
        
        if dd2>0.2:
            MR_y = np.exp(-(YY-Y0)**2/(0.5*dY**2))
        else:
            MR_y = (1 - np.tanh((YY-Y0-dY/2)/cc))/2 + (1 + np.tanh((YY-Y0+dY/2)/cc))/2-1
            
        h_MR = MR_y*h_MR
        
    return h_MR, MR_x
             
            
def add_zonal_ridge(XX, YY, Y0, dY, cc, H_ZR, X0=np.inf, dX=np.inf):   
    
    #print(H)

    h_ZR, ZR_y = add_meridional_ridge(YY, XX, Y0=X0, dY=dX, cc=cc, H_MR=H_ZR, X0=Y0, dX=dY)
            
    return h_ZR, ZR_y   
            
            
def add_zonal_channel(XX, YY, Y0, dY, cc, H_sill, H, MR_x=None):
    
    """
    Function to carve out channel/outlet out of a meridional ridge/continent permit zonal flow. 
    
    XX, YY: arrays defining grid
    Y0: latitudinal center of the channel
    DY: latitudinal width of the channel
    cc: steepness of the channel slopes (half-width)
    H_sill: sill depth of the channel (relative to seafloor)
    H: max seafloor depth
    MR_x: zonal profile of ridge
    
    """
    
    dd = cc/dY
    # dd: measures the meridional extent of the sill. If the sill is narrow, the cut is implemented
    # using a gaussian function. Otherwise, a set of tanh functions is used.
    
    #print(dd)
#     print("H: %s" %H)
#     print("H_sill: %s" %H_sill)
    
    if dd>0.2:
        print("using gaussian funcs...")
        # create ridge using gaussian function spanning the length of the domain
        h_cut = (H_sill-H)*np.exp(-(YY-Y0)**2/(0.5*dY**2)) 
        
    else:
        print("using tanh funcs...")
        # cut out channel using combination of tanh functions
        cut1 = (H-H_sill)*(1 + np.tanh((YY-Y0-dY/2)/cc))/2
        cut2 = (H-H_sill)*(1 - np.tanh((YY-Y0+dY/2)/cc))/2
        h_cut = cut1+cut2+H_sill-H # spans full meridional extent

    # limit cut in x-direction by wrapping in function used to define meridional continent
    if MR_x is not None:
        h_cut = MR_x*h_cut

    return h_cut


def add_meridional_channel(XX, YY, X0, dX, cc, H_sill, H, MR_y=None):
    
    """
    Function to carve out channel/outlet out of a meridional ridge/continent permit zonal flow. 
    This is a wrapper for add_zonal_channel()
    
    XX, YY: arrays defining grid
    X0: longitudinal center of the channel
    DX: longitudinal width of the channel
    cc: steepness of the channel slopes (half-width)
    H_sill: sill depth of the channel (relative to seafloor)
    H: max seafloor depth
    MR_y: meridional profile of ridge
    
    """

    # flip coordinates and use function to create zonal ridge
    h_cut = add_zonal_channel(XX=YY, YY=XX, Y0=X0, dY=dX, cc=cc, H_sill=H_sill, H=H, MR_x=MR_y)


    return h_cut


def add_meridional_canyon_old(XX, YY, X0, dX, cc, H_MR, Y0, dY):
    
    """
    function creates a canyon in the shelf, centered at X0, with longitudinal extent of dX and max depth of H_MR.
    The steepness of the slopes (in the x-direction) is controlled by cc. The latitudinal extent is infinite 
    by default, but can be constrained by Y0 and dY.
   
    """
    
    dd = cc/dX
    # dd: measures the prominence of the plateau. If the plateau is narrow, use a gaussian function.
    # otherwise, use tanh function. This is necessary since the tanh functions (as implemented here) can 
    # grossly underestimate the height of narrow features
    
    if dd>0.2:
        print("using gaussian funcs...")
        # create ridge using gaussian function spanning the length of the domain
        MR_x = np.exp(-(XX-X0)**2/(0.5*dX**2)) 
        h_MR = H_MR*MR_x
        
    else:
        print("using tanh funcs...")
        # combine two tanh functions to get meridional ridge
        MR_x = (1 - np.tanh((XX-X0-dX/2)/cc))/2 + (1 + np.tanh((XX-X0+dX/2)/cc))/2
        MR_x -= 1 # need to offset by H to get the right height
        h_MR = H_MR*MR_x
        
        
    print("test...")
    print(YY.min())
    print(YY.max())
    print(Y0)
    print(dY)
   
    # if ridge is finite, limit meridional extent
    if YY.min()<=Y0<YY.max() and dY<(YY.max()-YY.min()):
        
        # wrap in meridional function 
        dd2 = cc/dY
        
        if dd2>0.2:
            MR_y = np.exp(-(YY-Y0)**2/(0.5*dY**2))
        else:
            MR_y = (1 - np.tanh((YY-Y0-dY/2)/cc))/2 + (1 + np.tanh((YY-Y0+dY/2)/cc))/2-1
            
        h_MR = MR_y*h_MR

        
def add_meridional_canyon(XX, YY, X0, dX, cc, H_MR, h, Y0, Yend, ss):
    
    """
    function creates a canyon in the shelf, centered at X0, with longitudinal extent of dX and max depth of H_MR.
    The steepness of the slopes (in the x-direction) is controlled by cc. The latitudinal extent is infinite 
    by default, but can be constrained by Y0 and dY.
   
   
    H_MR is the max depth of the canyon
    """
    
    dd = cc/dX
    # dd: measures the prominence of the plateau. If the plateau is narrow, use a gaussian function.
    # otherwise, use tanh function. This is necessary since the tanh functions (as implemented here) can 
    # grossly underestimate the height of narrow features
    
    if dd>0.2:
        print("using gaussian funcs...")
        # create ridge using gaussian function spanning the length of the domain
        MR_x = np.exp(-(XX-X0)**2/(0.5*dX**2)) 
        h_MR = H_MR*MR_x
        
    else:
        print("using tanh funcs...")
        # combine two tanh functions to get meridional canyon
        MR_x = (1 - np.tanh((XX-X0-dX/2)/cc))/2 + (1 + np.tanh((XX-X0+dX/2)/cc))/2
        MR_x -= 1 # need to offset by H to get the right height
        h_MR = H_MR*MR_x
        
    
    # print(h_MR.shape)
    
    # Y0 = 50 #m
    # Yend = 300 #m
    # ss = 10*1000
    scale = -2 + (1+np.tanh((YY-Y0*1000)/ss)/2) + (1-np.tanh((YY-Yend*1000)/ss)/2)
       
    h_MR = scale*h_MR 
        
    return h_MR, MR_x
















            
            
            
            