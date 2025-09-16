#!/usr/bin/env python
# coding: utf-8
from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics in 3d view')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last from only, double to show plot(s)')
parser.add_argument('-x','--xCrossSection', nargs=1, type=float,default = None,
                    help='optional slice location [m]')
parser.add_argument('-y','--yCrossSection', nargs=1, type=float,default = None,
                    help='optional slice location [m]')
parser.add_argument('-z','--zDepth', nargs=1, type=float,default = None,
                    help='optional slice location [m]')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='option specification of views to plot [default = all]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-s','--shadow', action='count', default=0,
                    help='option of shadow for mélange [default (off)]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 5000
xCrossSection = 55000
zDepth = -50
plotDPI = 100
cleanPNGs = True
usePcolor = False
makeMovie = False
brightBergs = False

if(brightBergs):
    bergAlpha = 0.25
else:
    bergAlpha = 0.1

if(os.path.isfile('input/plotHelperLocal.py')):
    sys.path.append('input')
    from plotHelperLocal import *
    print('Found experiment plotting settings')
elif(os.path.isfile('../plotHelper.py')):
    sys.path.append('../')
    print('no custom plotting settings, using local default')
    from plotHelper import *
else:  
    print('no defaults found')
print('Plot DPI:',plotDPI,'; clean PNGs?',cleanPNGs)

if(args.xCrossSection != None):
    print('** Manual xSlice detected **')
    xCrossSection = args.xCrossSection
if(args.yCrossSection != None):
    print('** Manual ySlice detected **')
    yCrossSection = args.yCrossSection
if(args.zDepth != None):
    print('** Manual zDepth detected **')
    zDepth = args.zDepth

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

#Find Diagnostic file, iterate through them
maxStep = 0
sizeStep = 1e10
startStep = 1e10

for file in os.listdir('results'):
    # print(file)
    if "dynDiag.0" in file:
        words = file.split(".")
        # print(words[1])  
        if int(words[1]) > maxStep:
            maxStep = int(words[1])
        if int(words[1]) < startStep and int(words[1]) > 0:
            startStep = int(words[1])
        if abs(int(words[1]) - startStep) < sizeStep and abs(int(words[1]) - startStep) > 0:
            sizeStep = abs(int(words[1]) - startStep)

if(args.timeRange != None):
    startStep = args.timeRange[0] * 86400 / dt
    if(args.timeRange[1] != 0):
        maxStep = args.timeRange[1] * 86400 / dt

if((maxStep-startStep)/sizeStep > 60):   #if more than 60 frames, downscale to be less than 60
    dwnScale = int(np.ceil(((maxStep-startStep)/sizeStep)/60))
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale
print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Clean up old gifs and pngs
os.system('rm -f figs/cross_*.png')
os.system('rm -f figs/autoCross_*.gif')
os.system('rm -f figs/autoCross_*.m*')

#Import grid
x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = mds.rdmds("results/RC")

#Import bathymetry
if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))
   

# Print actual cross section values
xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
ySlice = np.argmin(np.abs(y[:,0] - yCrossSection))
zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
print('cross section is x =', x[0,xSlice], 'index', xSlice)
print('cross section is y =', y[ySlice,0], 'index', ySlice)
print('depth is z =', z[zSlice,0,0], 'index', zSlice)

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','BRGFlx','ptraceDiag','ptraceDiag']
    name = ["Temp", "Sal", "U", "W", "V", "BRGmltRt",'TracePlume','TraceBerg']
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]","[Vol Frac]","[Vol Frac]"]
else:
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

#In KM for x,y for better axis labeling
x = x/1000
y = y/1000

if(usePcolor):
    print('pcolor not supported for this function yet, using contourf')

ghostAlpha = .1

#NOTE matplotlib x and y and MITgcm x,y are FLIPPED below. Be careful.
if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False
if(args.kValues == None):
    kList = range(len(name))
else:
    kList = args.kValues
for k in kList:
    print('\t',name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111, projection='3d',computed_zorder=False)
        data = mds.rdmds("results/%s"%(dynName[k]), i)
        if(args.shadow > 0): #enable berg shadows here
            dataBergs = mds.rdmds("results/BRGFlx",i)
            if(dataBergs.shape[0] < 6):
                dataBergPlot = dataBergs[0,:,:,:]
                dataBergPlot[dataBergPlot != 0] = .2
            else:
                dataBergPlot = dataBergs[5,:,:,:]
            dataBergPlot[dataBergPlot == 0] = np.nan
        kk = k
        if k == 0:
            lvl = tempRange
            cm = tempCmap
        elif k == 1:
            lvl = saltRange
            cm = saltCmap
        elif k == 2:
            lvl = uRange
            cm = uCmap
        elif k == 3:
            lvl = wRange
            cm = wCmap
        elif k == 4:
            lvl = vRange
            cm = vCmap
        elif k == 5:
            lvl = meltRange
            cm = meltCmap
            kk = 2
        elif k == 6:
            lvl = plumeTracerRange
            cm = plumeTracerCmap
            kk = 0
        elif k == 7:
            lvl = bergTracerRange
            cm = bergTracerCmap
            kk = 1

        #Long profile
        XX,ZZ = np.meshgrid(np.squeeze(x[ySlice,:]),np.squeeze(z))
        cp = ax.contourf(
        np.squeeze(data[kk, :, ySlice, :]),
        XX,
        ZZ,
        levels=lvl,
        extend="both",
        cmap=cm,
        alpha= .9,
        zdir='x',offset=y[ySlice,0],zorder=2
        )
        if(args.shadow > 0):
            cp2 = plt.contourf(
                np.squeeze(dataBergPlot[:, ySlice, :]),
                XX,
                ZZ,
                [.4,.6,.8,.9,.95],
                extend="both",
                alpha=bergAlpha,
                cmap='cmo.gray',
                zdir='x',offset=y[ySlice,0],zorder=2)

        #viewers close width half profile
        YY,ZZ = np.meshgrid(np.squeeze(y[:,xSlice]),np.squeeze(z))    
        cp = ax.contourf(
            YY[:,0:ySlice+1],
            np.squeeze(data[kk, :, 0:ySlice+1, xSlice]),
            ZZ[:,0:ySlice+1],
            levels=lvl,
            extend="both",
            cmap=cm,
            alpha= .9,
            zdir='y',offset=x[0,xSlice],zorder=3
        )

        if(args.shadow > 0):
            cp2 = plt.contourf(
                YY[:,0:ySlice+1],
                np.squeeze(dataBergPlot[:, 0:ySlice+1, xSlice]),
                ZZ[:,0:ySlice+1],
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=bergAlpha,
                cmap='cmo.gray',
                zdir='y',offset=x[0,xSlice],zorder=3)
       
        #viewer's far width half profile
        YY,ZZ = np.meshgrid(np.squeeze(y[:,xSlice]),np.squeeze(z))    
        cp = ax.contourf(
            YY[:,ySlice:],
            np.squeeze(data[kk, :, ySlice:, xSlice]),
            ZZ[:,ySlice:],
            levels=lvl,
            extend="both",
            cmap=cm,
            alpha= .9,
            zdir='y',offset=x[0,xSlice],zorder=0
        )
        if(args.shadow > 0):
            cp2 = plt.contourf(
                YY[:,ySlice:],
                np.squeeze(dataBergPlot[:, ySlice:, xSlice]),
                ZZ[:,ySlice:],
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=bergAlpha,
                cmap='cmo.gray',
                zdir='y',offset=x[0,xSlice],zorder=0)
        
        #map view ghost
        #Front
        ghost = ax.contourf(
            np.squeeze(y[0:ySlice+1,0:xSlice+1]),
            np.squeeze(x[0:ySlice+1,0:xSlice+1]),
            np.ones(np.shape(x[0:ySlice+1,0:xSlice+1])),
            levels=[0, 1, 2],
            alpha=ghostAlpha,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=4
        )
        #Back
        ghost = ax.contourf(
            np.squeeze(y[ySlice:,xSlice:]),
            np.squeeze(x[ySlice:,xSlice:]),
            np.ones(np.shape(x[ySlice:,xSlice:])),
            levels=[0, 1, 2],
            alpha=ghostAlpha,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=-1
        )
        #Left
        if(xSlice > 1): #skip if slice is forward most slice
            ghost = ax.contourf(
                np.squeeze(y[ySlice:,0:xSlice]),
                np.squeeze(x[ySlice:,0:xSlice]),
                np.ones(np.shape(x[ySlice:,0:xSlice])),
                levels=[0, 1, 2],
                alpha=ghostAlpha,
                cmap='cmo.gray',
                zdir='z',offset=z[zSlice,0,0],zorder=1
            )
        #Right
        ghost = ax.contourf(
            np.squeeze(y[0:ySlice+1,xSlice:]),
            np.squeeze(x[0:ySlice+1,xSlice:]),
            np.ones(np.shape(x[0:ySlice+1:,xSlice:])),
            levels=[0, 1, 2],
            alpha=ghostAlpha,
            cmap='cmo.gray',
            zdir='z',offset=z[zSlice,0,0],zorder=2
        )
        #Map view, projected onto bottom of frame
        cp = ax.contourf(
            np.squeeze(y),
            np.squeeze(x),
            np.squeeze(data[kk, zSlice, :, :]),
            levels=lvl,
            extend="both",
            cmap=cm,
            alpha= .9,
            zdir='z',offset=z.min(),zorder=-1
        )
        if(args.shadow > 0):
            cp2 = plt.contourf(
                np.squeeze(y),
                np.squeeze(x),
                np.squeeze(dataBergPlot[zSlice, :, :]),
                [.05,.1,.2,.4,.8],
                extend="min",
                alpha=bergAlpha,
                cmap='cmo.gray',
                zdir='z',offset=z.min(),zorder=-1)

        cbar = fig.colorbar(cp,pad=0.0,fraction=.04,orientation='horizontal')
        cbar.set_label(cbarLabel[k])
        
        #Topography along fjord
        ax.plot(x[ySlice,:xSlice],topo[ySlice,:xSlice],color='black',zdir='x',zs=y[ySlice,0],zorder=2)        
        ax.plot(x[ySlice,xSlice:],topo[ySlice,xSlice:],color='black',zdir='x',zs=y[ySlice,0],zorder=2)
        
        #Topography across fjord
        ax.plot(y[ySlice:,xSlice],topo[ySlice:,xSlice],color='black',zdir='y',zs=x[0,xSlice],zorder=0)
        ax.plot(y[:ySlice+1,xSlice],topo[:ySlice+1,xSlice],color='black',zdir='y',zs=x[0,xSlice],zorder=3)
        ax.invert_xaxis()
        plt.title("%s at x,y,z (%i,%i,%i) at %.02f days" % (name[k], x[0,xSlice]*1000,y[ySlice,0]*1000,z[zSlice,0,0], i/86400.0*dt))
        ax.set_xlabel('Width [km]')
        ax.set_ylabel('Along [km]')
        ax.set_zlabel('Depth [m]')
        
        ax.axes.set_xlim3d(left=y.max(), right=y.min())
        ax.axes.set_ylim3d(bottom=x.min(), top=x.max()) 
        ax.axes.set_zlim3d(bottom=z.min(), top=z.max()) 
        ax.set_box_aspect([2,6,1]) # y,x,z
        ax.view_init(elev=25., azim=-40)
        j = i/sizeStep
        plt.tight_layout()
        str = "figs/cross_%s%05i.png" % (name[k],j)
        plt.tight_layout()
        plt.savefig(str, format='png', dpi=plotDPI)
        if(args.quick > 1):
            plt.show()
        plt.close()
    if(args.quick == 0):
        os.system('magick -delay %f figs/cross_%s*.png -colors 256 -depth 256 figs/autoCross_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))
        if(makeMovie):
            os.system('ffmpeg -r %f -i figs/cross_%s%%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p figs/autoCross_%s.mp4' %(80/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/cross_*.png')



