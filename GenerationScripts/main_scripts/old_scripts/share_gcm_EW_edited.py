'''
This file contains some global variables and functions that can be accessed from outside.
'''

def main():
    print("Don't run this directly!")

if __name__ == "__main__":
    main()

# Modules imported
import numpy as np
import pythonstartup as ps
import pickle
import glob
from mitgcm import mds    # mds.get_dim(dir, varn) returns a dictionary that contains the grid information
import MITgcm_Fortran as mfortran
import sys
from scipy import interp
from grid import *
from numba import jit


# @jit(forceobj=True)
# def calc_heat_loss(flux, lat, dx2d, dyz):
#     '''
#     Here we calculate the heat loss to the south of 60S
#     '''
#     FF   = 0
#     for i in range(0, ny):
#         if lat[i] > -40:
#             continue
#         for j in range(0, nx):
#             FF += flux[i,j] * dx2d[i,j] * dyz
#     return FF

def obtain_files(rootdir, var):
    # This will return a list of files pertain to the variable
    # specified by $var. This will also call mds to return the grid
    # information
    list_file    = glob.glob(rootdir + var + "*.data")
    list_file.sort()
    grid         = mds.get_dim(rootdir, var)
    dim          = grid["dim"]
    return list_file,dim

class MITgcm:
    '''
    This is a post-processing module that will return a bunch of variables:
    (1) Global volume-averaged temperature
    (2) Global volume-averaged salinity
    (3) Eulerian-mean MOC for each output step
    (4) Mean variables for the last N files
    (5) Isopycnal overturning circulation calculated over the last 
    M files (M can be equal to N)
    '''

    def __init__(self, rootdir, num=-1):
        '''MITgcm(rootdir, num=?).  num indicates the number of the last
        step. default -1 means that it will be the real last one.

        '''
        self.rootdir    = rootdir + "/"
        self.Series     = {}
        self.Mean       = {}
        self.MOC        = {}
        self.ROC        = {}
        self.num        = num
        # load the grid information
        # atlantic        = np.zeros((ny, nx))
        # atlantic[21:, :30] = 1
        # pacific         = np.zeros((ny, nx))
        # pacific[21:, 30:] = 1
        with open("Basin_Mask.bin", "rb") as fid:
            D = pickle.load(fid)

        self.BMask      = D
        self.ATL        = D["ATL"]
        self.PAC        = D["IND"] + D["PAC"]
        self.lat        = y
        self.z0         = np.insert(np.cumsum(dr), 0, 0)
        dx1d            = R0*np.cos(y*deg2rad)*dx*deg2rad
        self.dx2d       = np.repeat(dx1d.reshape(ny, 1), nx, axis=1)
        self.dyz        = R0*dy*deg2rad
        
    
    def calc_time_series(self, start=-1):
        # This calculate the time series of global volume averaged data
        # Two variables: T & S
        list_T,dim_T    = obtain_files(self.rootdir, "Ttave")
        list_S,dim_S    = obtain_files(self.rootdir, "Stave")
        list_V,dim_V    = obtain_files(self.rootdir, "vVeltave")
        list_vh,dim_vh  = obtain_files(self.rootdir, "layers_VH-tave")
        list_hs,dim_hs  = obtain_files(self.rootdir, "layers_Hs-tave")
        list_flux,dim_flux = obtain_files(self.rootdir, "tFluxtave")
        list_eta,dim_eta= obtain_files(self.rootdir, "ETAtave")
        list_U,dim_U    = obtain_files(self.rootdir, "uVeltave")
        [nz,ny,nx]      = dim_T
        nT    = len(list_T)
        nS    = len(list_S)
        nV    = len(list_V)
        if nT == nS and nT == nV:
            print("We are going to start calculating the global-volume average of Temperature and Salinity!!")
        else:
            print("Warning: number of files not consistent between T and S!! Some of the file not accounted.")
        
        dx2d   = self.dx2d
        if self.num > 0:
            m = min(nT, nS, self.num)
        else:
            m     = min(nT, nS)
        T     = np.zeros(m)
        S     = np.zeros(m)
        mmoc  = np.zeros(m)    # maximum AMOC streamfunction - eulerian mean
        zmoc  = np.zeros(m)    # depth of moc - eulerian mean
        amoc30= np.zeros(m)
        gmoc30= np.zeros(m)
        pmoc30= np.zeros(m)    # pmoc at the max amoc isopycnal
        k30   = ps.find_closest(self.lat, -30) + 1
        nd    = dim_vh[0]
        amoct30 = np.zeros((nd+1, m))
        gmoct30 = np.zeros((nd+1, m))
        pmoct30 = np.zeros((nd+1, m))
        za30    = np.zeros((nd+1, m))
        zp30    = np.zeros((nd+1, m))
        zg30    = np.zeros((nd+1, m))
        zatl    = np.zeros(m)
        zpac    = np.zeros(m)
        zmoc30  = np.zeros(m)    # depth of the maximum at 30 S
        eta     = np.zeros((m, ny, nx))
        ITF     = np.zeros(m)
        TFlux   = {}
        for  basin in ["ATL", "IND", "PAC", "SO"]:
            TFlux[basin] = np.zeros(m)

        
        #FLUX    = np.zeros(m)    # here is the heat loss in the Southern Ocean. 
        for i in range(max(0, start), m):
            print("Processing No.%d"%i)
            TEMP   = ps.read_bin(list_T[i], dtype=">f").reshape(dim_T)
            SALT   = ps.read_bin(list_S[i], dtype=">f").reshape(dim_S)
            VVEL   = ps.read_bin(list_V[i], dtype=">f").reshape(dim_V)
            UVEL   = ps.read_bin(list_U[i], dtype=">f").reshape(dim_U)
            vh     = ps.read_bin(list_vh[i],dtype=">f").reshape(dim_vh)
            hs     = ps.read_bin(list_hs[i],dtype=">f").reshape(dim_hs)
            tflux  = ps.read_bin(list_flux[i], dtype=">f").reshape(dim_flux)
            eta[i] = ps.read_bin(list_eta[i], dtype=">f").reshape(dim_eta)
            # Let's do the calculation by calling another external
            # function written in fortran
            T[i],S[i]   = mfortran.calc_volume_average(TEMP,SALT,dx2d, self.dyz,
                                                       dz, nz,ny,nx)
            
            mmoc[i],mocz= mfortran.calc_eulerian_amoc(VVEL, self.ATL, 
                                                      self.lat, dx2d, 
                                                      dz, nz, ny, nx)
            # calculate the MOC depth
            # 1: find max and min location
            kmax        = ps.find_last(mocz == max(mocz))
            kmin        = ps.find_last(mocz == min(mocz[kmax:]))
            if kmax >= kmin:
                zmoc[i] = np.nan
                print("AMOC not well defined!")
            else:
                zmoc[i]     = interp(0, mocz[kmin:kmax:-1], 
                                     self.z0[kmin:kmax:-1])
            # Now calculate the MOC strength at 30S
            amoct30[:,i],gmoct30[:,i],pmoct30[:,i],amoc30[i],gmoc30[i], \
                pmoc30[i],za30[:,i],zg30[:,i],zp30[:,i] = \
                    mfortran.calc_moc_30(vh, hs, self.ATL, self.PAC, \
                                         k30, dx2d, nd, ny, nx)

            # calculate the heat loss
            #FLUX[i]  = calc_heat_loss(tflux, y, dx2d, self.dyz)

            # Calculate the depth of the isotherm
            # if i == 0:
            #     nmax = amoct30[:,0].argmax() # here is the level where the AMOC reaches maximum
            # zmax = np.sum(hs[:nmax, :, :], axis=0)
            # zatl[i] = np.sum(zmax*self.ATL*dx2d)/np.sum(self.ATL*dx2d)
            # zpac[i] = np.sum(zmax*self.PAC*dx2d)/np.sum(self.PAC*dx2d)

            # mmax      = amoct30[:,i].argmax()
            # zmax1     = np.sum(hs[:mmax, k30, :], axis=0)
            # zmoc30[i] = np.sum(zmax1*self.ATL[k30,:])/np.sum(self.ATL[k30,:])
            # calculate ITF
            for ky1 in range(70, 74):
                ITF[i]   += np.sum(UVEL[:, ky1, 120]*self.dyz*dz)

            # Calculate flux!!!
            for basin in ["ATL", "IND", "PAC", "SO"]:
                TFlux[basin][i] = np.sum(tflux*self.BMask["TAREA"]*self.BMask[basin])   # this is in W
            

        self.Series["TEMP"] = T
        self.Series["SALT"] = S
        self.Series["mmoc"] = mmoc
        self.Series["zmoc"] = zmoc
        self.Series["amoc30"] = amoc30
        self.Series["gmoc30"] = gmoc30
        self.Series["pmoc30"] = pmoc30
        self.Series["amoct30"]= amoct30/1e6
        self.Series["pmoct30"]= pmoct30/1e6
        self.Series["gmoct30"]= gmoct30/1e6
        self.Series["za30"]   = za30
        self.Series["zg30"]   = zg30
        self.Series["zp30"]   = zp30
        #self.Series["heat_loss"]= FLUX
        self.Series["zatl"] = zatl
        self.Series["zpac"] = zpac
        self.Series["zmoc30"]=zmoc30
        self.Series["eta"]  = eta
        self.Series["ITF"] = -ITF/1e6
        self.Series["HeatFlux"]=TFlux

    def calc_average(self, lastn=10, center=None, window=None):      
        # Calculate the last $lastn files
        # Default is 10
        varn    = ["tFluxtave","Ttave", "ETAtave", "layers_VH-tave", "Stave",
                   "layers_Hs-tave", "layers_Hw-tave", "layers_UH-tave", "vVeltave", "wVeltave",
                   "Convtave", "uVeltave"]  #, "PTRtave01", "Tdiftave"]
        for var in varn:
            print("Calculating %s"%var)
            list_file, dim  = obtain_files(self.rootdir, var)
            vv              = np.zeros(dim)
            nn              = len(list_file)
            if center is None:
                if self.num < 0:
                    for i in range(-min(lastn, nn), 0):
                        vt    = ps.read_bin(list_file[i], dtype=">f").reshape(dim)
                        if var == "Stave":
                            vt[vt<0] = 0
                        vv   += vt/lastn
                else:
                    for i in range(self.num-lastn, self.num):
                        vt    = ps.read_bin(list_file[i], dtype=">f").reshape(dim)
                        if var == "Stave":
                            vt[vt<0]  = 0     # This is weird but it happens, probably some 
                        vv   += vt/lastn
            else:
                dn  = int(window/2)
                for i in range(center-dn, center+dn+1):
                    vt    = ps.read_bin(list_file[i], dtype=">f").reshape(dim)
                    if var == "Stave":
                        vt[vt<0] = 0
                    vv   += vt/window

            if var in ["layers_VH-tave","layers_Hs-tave", "layers_Hw-tave", "layers_UH-tave"]:
                #self.Mean[var] = vv[::-1, :, :]
                self.Mean[var] = vv
            else:
                self.Mean[var] = vv

    def calc_moc(self):
        # This will return the Eulerian-mean overturning circulation using the calculated mean data
        if "vVeltave" in self.Mean.keys():
            pass
        else:
            print("The mean value should be calculated first before calculating the MOC! Aborting...")
            sys.exit()
        
        V     = self.Mean["vVeltave"]
        dx2d  = self.dx2d
        nz    = len(dz) + 1
        dz2   = np.insert(np.array(dz), 0, 0)
        # Grid for MOC
        zz    = np.cumsum(dz2)
        lat   = y
        # MOC
        amoc  = np.zeros((nz, ny))
        gmoc  = np.zeros((nz, ny))


        ATL   = self.ATL
        # Calculating GMOC
        for i in range(1, nz):
            gmoc[i, :] = gmoc[i-1,:] + dz[i-1] * np.sum(V[i-1,:, :]*dx2d, axis=1)
            amoc[i, :] = amoc[i-1,:] + dz[i-1] * np.sum(V[i-1,:, :]*dx2d*ATL, axis=1)
        
        self.MOC["amoc"] = amoc/1e6
        self.MOC["gmoc"] = gmoc/1e6
        self.MOC["lat"]  = lat
        self.MOC["zz"]   = zz
        
    # def calc_roc(self, process_years=None):
    #     if "layers_VH-tave" in self.Mean.keys() or process_years is not None:
    #         pass
    #     else:
    #         print("The mean value should be calculated first before calculating the MOC! Aborting...")
    #         sys.exit()
        
    #     dx2d  = self.dx2d
    #     nz    = len(dz)
        
    #     # load pd
    #     pd_r  = ps.read_bin(self.rootdir + "layers1TH.data", dtype=">f")
    #     pd    = pd_r[::-1]
    #     lat   = y
    #     nd    = len(pd)

    #     # Now we calculate ROC  --- This is much simpler because of the regular grid.
    #     if process_years is None:
    #         vh    = self.Mean["layers_VH-tave"]
    #         Hs    = self.Mean["layers_Hs-tave"]
    #     else:
    #         # In this case, I shall calculate vh and Hs from the original data
    #         nyear        = len(process_years)
    #         list_hs,dim0 = obtain_files(self.rootdir, "layers_Hs-tave")
    #         list_vh,dim1 = obtain_files(self.rootdir, "layers_VH-tave")
    #         vh1 = np.zeros(dim1)
    #         Hs1 = np.zeros(dim0)
    #         for k in process_years:
    #             vh1 += ps.read_bin(list_vh[k], dtype=">f").reshape(dim1)/nyear
    #             Hs1 += ps.read_bin(list_hs[k], dtype=">f").reshape(dim0)/nyear

    #         vh = vh1[::-1,:,:]
    #         Hs = Hs1[::-1,:,:]
            
    #     gmoc  = np.zeros((nd, ny))
    #     amoc  = np.zeros((nd, ny))
    #     pmoc  = np.zeros((nd, ny))
    #     vg    = np.zeros((nd, ny))   # Volume above 
    #     va    = np.zeros((nd, ny))   # Volume to calculate the depth of isopycnals
    #     vp    = np.zeros((nd, ny))
    #     ATL   = self.ATL
    #     PAC   = self.PAC
    #     # Let's integrate downward!!!!
    #     for i in range(1,nd):
    #         gmoc[i,:] = gmoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d, axis=1)
    #         amoc[i,:] = amoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d * ATL, axis=1)
    #         pmoc[i,:] = pmoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d * PAC, axis=1)
    #         vg[i,:]   = vg[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d, axis=1)
    #         va[i,:]   = va[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d * ATL, axis=1)
    #         vp[i,:]   = vp[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d * PAC, axis=1)
            
    #     # Now we calculate the depth by interpolation
    #     # Prepare the data
    #     vgp  = np.zeros((nz+1, ny))
    #     vap  = np.zeros((nz+1, ny))
    #     vpp  = np.zeros((nz+1, ny))
    #     zz   = np.cumsum(np.insert(dz, 0, 0))
    #     if process_years is None:
    #         TT    = self.Mean["Ttave"]    # used to decide if it is ocean or land
    #     else:
    #         list_T,dim = obtain_files(self.rootdir, "Ttave")
    #         TT    = ps.read_bin(list_T[0], dtype=">f").reshape(dim)
            
    #     mask = np.invert(TT==0)          # mask for ocean points
    #     dz3 = np.repeat(np.repeat(np.array(dz).reshape(nz, 1, 1), ny,
    #                               axis=1), nx, axis=2) * mask
    #     for i in range(0, nz):
    #         vgp[i+1, :] = vgp[i,:] + np.sum(dz3[i,:,:] * dx2d, axis=1)
    #         vap[i+1, :] = vap[i,:] + np.sum(dz3[i,:,:] * dx2d * ATL, axis=1)
    #         vpp[i+1, :] = vpp[i,:] + np.sum(dz3[i,:,:] * dx2d * PAC, axis=1)
        
    #     # This is the depth of isopycnals here
    #     za   = np.zeros((nd, ny))
    #     zg   = np.zeros((nd, ny))
    #     zp   = np.zeros((nd, ny))
    #     for j in range(0, ny):
    #         if vap[1,j] > 0:
    #             za[:,j] = interp(va[:,j], vap[:,j], zz)
    #         if vgp[1,j] > 0:
    #             zg[:,j] = interp(vg[:,j], vgp[:,j], zz)
    #         if vpp[1,j] > 0:
    #             zp[:,j] = interp(vp[:,j], vpp[:,j], zz)
            
    #     lat2   = np.repeat(lat.reshape(1,ny), nd, axis=0)
    #     self.ROC["amoc"]  = amoc/1e6
    #     self.ROC["gmoc"]  = gmoc/1e6
    #     self.ROC["pmoc"]  = pmoc/1e6
    #     self.ROC["pd"]    = pd
    #     self.ROC["lat"]   = lat
    #     self.ROC["za"]    = za
    #     self.ROC["zg"]    = zg
    #     self.ROC["zp"]    = zp
    #     self.ROC["lat2"]  = lat2

    #     #return self.ROC
    def calc_roc(self):
        if "layers_VH-tave" in self.Mean.keys():
            pass
        else:
            print("The mean value should be calculated first before calculating the MOC! Aborting...")
            sys.exit()
        dx2d  = self.dx2d
        nz    = len(dz)

        # load pd
        pd    = ps.read_bin(self.rootdir + "layers1RHO.data", dtype=">f")
        lat   = self.lat
        nd    = len(pd)
        print(f"number of layers: {nd}")

        # Now we calculate ROC  --- This is much simpler because of the regular grid.
        vh    = self.Mean["layers_VH-tave"]
        Hs    = self.Mean["layers_Hs-tave"]

        gmoc  = np.zeros((nd, ny))
        amoc  = np.zeros((nd, ny))
        pmoc  = np.zeros((nd, ny))
        vg    = np.zeros((nd, ny))   # Volume above
        va    = np.zeros((nd, ny))   # Volume to calculate the depth of isopycnals
        vp    = np.zeros((nd, ny))
        ATL   = self.ATL
        PAC   = self.PAC
        # Let's integrate downward!!!!
        for i in range(1,nd):
            gmoc[i,:] = gmoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d, axis=1)
            amoc[i,:] = amoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d * ATL, axis=1)
            pmoc[i,:] = pmoc[i-1,:] + np.sum(vh[i-1,:,:] * dx2d * PAC, axis=1)
            vg[i,:]   = vg[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d, axis=1)
            va[i,:]   = va[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d * ATL, axis=1)
            vp[i,:]   = vp[i-1,:]   + np.sum(Hs[i-1,:,:] * dx2d * PAC, axis=1)

        # Now we calculate the depth by interpolation
        # Prepare the data
        vgp  = np.zeros((nz+1, ny))
        vap  = np.zeros((nz+1, ny))
        vpp  = np.zeros((nz+1, ny))
        zz   = np.cumsum(np.insert(dz, 0, 0))
        V    = self.Mean["vVeltave"]    # used to decide if it is ocean or land
        mask = np.invert(V==0)          # mask for ocean points
        dz3  = np.repeat(np.repeat(dz.reshape(nz, 1, 1), ny, axis=1), nx, axis=2) * mask
        for i in range(0, nz):
            vgp[i+1, :] = vgp[i,:] + np.sum(dz3[i,:,:] * dx2d, axis=1)
            vap[i+1, :] = vap[i,:] + np.sum(dz3[i,:,:] * dx2d * ATL, axis=1)
            vpp[i+1, :] = vpp[i,:] + np.sum(dz3[i,:,:] * dx2d * PAC, axis=1)

        # This is the depth of isopycnals here
        za   = np.zeros((nd, ny))
        zg   = np.zeros((nd, ny))
        zp   = np.zeros((nd, ny))
        for j in range(0, ny):
            if vap[1,j] > 0:
                za[:,j] = interp(va[:,j], vap[:,j], zz)
            if vgp[1,j] > 0:
                zg[:,j] = interp(vg[:,j], vgp[:,j], zz)
            if vpp[1,j] > 0:
                zp[:,j] = interp(vp[:,j], vpp[:,j], zz)
        lat2   = np.repeat(lat.reshape(1,ny), nd, axis=0)
        self.ROC["amoc"]  = amoc/1e6
        self.ROC["gmoc"]  = gmoc/1e6
        self.ROC["pmoc"]  = pmoc/1e6
        self.ROC["pd"]    = pd
        self.ROC["lat"]   = lat
        self.ROC["za"]    = za
        self.ROC["zg"]    = zg
        self.ROC["zp"]    = zp
        self.ROC["lat2"]  = lat2            
