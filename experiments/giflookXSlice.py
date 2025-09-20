from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics at xCrossSection')
parser.add_argument('-x','--xCrossSection', nargs=1, default=None, type=float,
                    help='optional x location [m]')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last from only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='option specification of views to plot [default = all]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
showDensity = True
usePcolor = False

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
    print('** Manual xCrossSection detected **')
    xCrossSection = args.xCrossSection

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

maxStep = 0
sizeStep = 1e10
startStep = 1e10

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

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
    maxStep = args.timeRange[1] * 86400 / dt

if((maxStep-startStep)/sizeStep > args.numFrames):   #if more than numFrames, downscale to be less than numFrames
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/args.numFrames)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

os.system('rm -f figs/sideX*.png')
# os.system('rm -f figs/autosideX*.gif')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(isBerg):
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
    bergMaskNums = bergMaskNums.reshape(np.shape(x))
    bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
    bergsPerCell = bergsPerCell.reshape(np.shape(x))
    maxDepth = np.zeros(np.shape(x))

    ## deepest contour
    # for i in range(np.shape(x)[1]):
    #     for j in range(np.shape(x)[0]):
    #         bergCount = int(bergsPerCell[j,i])
    #         if(bergMask[j,i] == 1 and bergCount > 0):  #only go in if bergs here
    #             depthFile = 'input/iceberg_depth_%05i.txt' % int(bergMaskNums[j,i])
    #             depths = np.zeros(bergCount)
    #             with open(depthFile,'r') as readFile:
    #                 ii = 0
    #                 for line in readFile:
    #                     if ii >= bergCount:
    #                         print('berg count mismatch in depth')
    #                         break
    #                     depths[ii] = float(line)
    #                     ii += 1
    #             readFile.close()
    #             maxDepth[j,i] = np.max(depths)


    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1


xSlice = np.argmin(np.abs(x[0,:] - xCrossSection))
print('cross section is x =', x[0,xSlice],'index', xSlice)

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'BRGFlx','ptraceDiag','ptraceDiag']
    name = ["Temp", "Sal", "U", "W", "V", "BRGmltRt",'TracePlume','TraceBerg']
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]","[Vol Frac]","[Vol Frac]"]
else:
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V"]
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]"]

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
        if(isBerg and os.path.isfile('results/BRGFlx.%010i.001.001.data' % i)):
            localBergs = True
        else:
            localBergs = False
        data = mds.rdmds("results/%s"%(dynName[k]), i)
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
            kk = k - 3
        elif k == 6:
            lvl = plumeTracerRange
            cm = plumeTracerCmap
            kk = 0
        elif k == 7:
            lvl = bergTracerRange
            cm = bergTracerCmap
            kk = 1
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(data[kk, :, :, xSlice]),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(data[kk, :, :, xSlice]),
                lvl,
                extend="both",
                cmap=cm,
            )
        if(showDensity and (dynName[k] == 'dynDiag')):
            salt = np.squeeze(data[1,:,:,xSlice])
            if(i == startStep): #only calc pressure once
                pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
            CT = gsw.CT_from_t(salt, data[0,:,:,xSlice], 0)
            density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
            densityLevels = np.linspace(25,28,16)
            cc = plt.contour(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(density),
                densityLevels,
                colors='black',
                linewidths=0.5,
                alpha=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        plt.plot(y[:,xSlice],topo[:,xSlice],color='black')
        # plt.plot(y[:,xSlice],ice[:,xSlice],color='gray')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        if(localBergs):
            # plt.plot(y[:,xSlice],-np.max(maxDepth,axis=1),color='gray',linestyle='dotted')
            cp2 = plt.contourf(
                np.squeeze(y[:,xSlice]),
                np.squeeze(z),
                np.squeeze(openFrac[:, :, xSlice]),
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')
        plt.xlabel('Across Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[kk, :, :, xSlice]),np.nanmax(data[kk, :, :, xSlice]),np.max(np.isnan(data[kk, :, :, xSlice]))))
        plt.ylabel('Depth [m]')
        plt.title("%s x = %i at %.02f days" % (name[k], x[0,xSlice], i/86400.0*dt))
        j = i/sizeStep + startStep
        
        str = "figs/sideX%s%05i.png" % (name[k],j)
        
        plt.savefig(str, format='png', dpi=plotDPI)
        if(args.quick > 1):
            plt.show()
        plt.close()

    if(args.quick == 0):
        locStr = ''
        timeStr = ''
        if(args.xCrossSection != None):
            locStr = f'{int(np.abs(args.xCrossSection[0]))}'
        if(args.timeRange != None):
            timeStr = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
        os.system(f'magick -delay {500/((maxStep-startStep)/sizeStep)} figs/sideX{name[k]}*.png -colors 256 -depth 256 figs/autosideX_{locStr}{name[k]}{timeStr}.gif') 

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/sideX*.png')
