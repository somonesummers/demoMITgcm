from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import gsw
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics width averaged')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last from only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='option specification of views to plot [default = all]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('-s','--shadow', action='count', default=0,
                    help='option of shadow for mélange [default (off)]')
parser.add_argument('-xl','--xlimit', nargs='?', type=int, default = None,
                    help='optional max x [default = entire frame]')
parser.add_argument('-dpi','--dpi', nargs='?', type=int, default = 125,
                    help='dpi spec [default = 125]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
cleanPNGs = True
showQuiver = False
showZeros = True
showDensity = True
makeMovie = False   

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

# if(args.yCrossSection != None):
#     print('** Manual ySlice detected **')
#     yCrossSection = args.yCrossSection
# #Overwrite local settings here if desired
# # usePcolor = True

print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

dt = 0.0   
for line in fileinput.input('input/data'):
        if "deltaT=" in line:
            dt = float(line[8:-2])
print('dt is loaded as', dt)

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

## Custom timeshots available here
# startStep = int(((maxStep - startStep)/startStep) * 1 / 2 )* startStep 
# maxStep = maxStep - 10 * sizeStep

if(args.timeRange != None):
    startStep = args.timeRange[0] * 86400 / dt
    maxStep = args.timeRange[1] * 86400 / dt

if((maxStep-startStep)/sizeStep > args.numFrames):   #if more than numFrames, downscale to be less than numFrames
    dwnScale = np.ceil(((maxStep-startStep)/sizeStep)/args.numFrames)
    print('Reducing time resolution by', dwnScale)
    sizeStep = sizeStep * dwnScale

print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)

#Decide if iceBerg data files exist
if(os.path.isfile('input/bergMask.bin')):
    isBerg = True
    print('Found icebergs for this run')
else:
    isBerg = False

# os.system('rm -f figs/side_*.png')
# os.system('rm -f figs/autoside_*.gif')
# os.system('rm -f figs/autoside_*.mov')

x = mds.rdmds("results/XC")
y = mds.rdmds("results/YC")
z = np.squeeze(mds.rdmds("results/RC"))

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

oldAveraging = True
if(isBerg):
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    bergMaskNums = np.fromfile('input/bergMaskNums.bin', dtype='>f8')
    bergMaskNums = bergMaskNums.reshape(np.shape(x))
    bergsPerCell = np.fromfile('input/numBergsPerCell.bin', dtype='>f8')
    bergsPerCell = bergsPerCell.reshape(np.shape(x))
    bergContaingingCells = int(np.sum(bergMask))
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
                if topo[j,i] == 0 and i > 1:
                    openFrac[:,j,i] = 0 #zero weight non-ocean cell
                    
if(np.min(openFrac) < 0 or np.max(openFrac) > 1):
    print(f'Error in openFrac min/max: {np.min(openFrac):0.2f}/{np.max(openFrac):0.2f}')
print('averaging over all cross sections')

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
    print("\t" + name[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(showQuiver):
            dataQuiv = mds.rdmds("results/dynDiag", i)
        if(isBerg):
            dataBergs = mds.rdmds("results/BRGFlx",i)
            if(dataBergs.shape[0] == 6):
                # print('ttest')
                openFrac[:,:,:] = dataBergs[5,:,:,:] #directly saved for these runs
                for jj in range(np.shape(x)[0]): #clean up non-berg parts of this mask
                    for ii in range(np.shape(x)[1]):
                        if bergMask[jj,ii] == 0:
                            openFrac[:,jj,ii] = 1
                        if topo[jj,ii] == 0 and ii > 1:
                            openFrac[:,jj,ii] = 0 #zero weight non-ocean cell
                # print(np.sum(openFrac[1,:,:],axis=0))
                # print(np.sum(openFrac2[1,:,:],axis=0))
                # oldAveraging = False
        if(args.shadow > 0): #enable berg shadows here
            if(dataBergs.shape[0] < 6):
                dataBergPlot = dataBergs[0,:,:,:]
                dataBergPlot[dataBergPlot != 0] = .2
            else:
                dataBergPlot = dataBergs[5,:,:,:]
            dataBergPlot[dataBergPlot == 0] = 1 #np.nan
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
        plt.figure(figsize=(12, 5))
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(np.average(data[kk, :, 1:-1, :], weights=openFrac[:,:,:],axis=1)),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(np.average(data[kk, :,:, :], weights=openFrac[:,:,:],axis=1)),
                lvl,
                extend="both",
                cmap=cm,
            )
        # if(k == 5): #special bounds to highlight refreezing areas
        #     cp.cmap.set_under('r')
        plt.plot(x[0,:],topo[int(np.shape(x)[0]/2),:],color='black')
        if(args.shadow > 0):
            cp2 = plt.contourf(
                x[0,:],
                np.squeeze(z),
                np.squeeze(np.nanmean(dataBergPlot[:, :, :],axis=1)),
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=.2,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')
        cbar = plt.colorbar(cp)
        cbar.set_label(cbarLabel[k])
        if(showDensity and (dynName[k] == 'dynDiag')):
            salt = np.squeeze(np.average(data[1,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1))
            if(i == startStep): #only calc pressure once
                pressure = -1 * np.ones(salt.shape) * 1020 * 9.81 * np.repeat(np.expand_dims(z,1), salt.shape[1], axis=1) /10e3
            CT = gsw.CT_from_t(salt, np.average(data[0,:,1:-1,:],weights=openFrac[:,1:-1,:],axis=1), pressure)
            density = gsw.rho(salt, CT, 0) - 1000 #in-stu density less 1000
            densityLevels = np.linspace(25,28,16)
            cc = plt.contour(
                np.squeeze(x[0,:]),
                np.squeeze(z),
                np.squeeze(density),
                densityLevels,
                colors='black',
                linewidths=0.5,
                alpha=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        if(showZeros):
            if( k == 2 or k == 6 or k == 7):
                cc = plt.contour(
                    np.squeeze(x[0,:]),
                    np.squeeze(z),
                    np.squeeze(np.average(data[kk, :, 1:-1, :], weights=openFrac[:,1:-1,:],axis=1)),
                    [0],
                    colors='gray',
                    linewidths=0.5,
                    alpha=0.5
                )
                plt.clabel(cc, inline=3, fontsize=8)
        if(args.xlimit != None):
            plt.xlim([0, args.xlimit])
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(data[kk, :, 1:-1, :]),np.nanmax(data[kk, :, 1:-1, :]),np.max(np.isnan(data[kk, :, 1:-1, :]))))
        plt.ylabel('Depth [m]')
        plt.title("%s Width Averaged at %.02f days" % (name[k], i/86400.0*dt))
        j = i/sizeStep
        
        str = "figs/sideAvg_%s%05i.png" % (name[k],j)
        plt.tight_layout()
        plt.savefig(str, format='png', dpi=args.dpi)
        if(args.quick > 1):
            plt.show()
        plt.close()

    
    if(args.quick == 0):
        os.system('magick -delay %f figs/sideAvg_%s*.png -colors 256 -depth 256 figs/sideAvg_%s.gif' %(500/((maxStep-startStep)/sizeStep), name[k], name[k]))
        if(makeMovie):
            os.system('ffmpeg -r %f -i figs/sideAvg_%s%%05d.png -c:v libx264 -r 30 -pix_fmt yuv420p figs/sideAvg_%s.mov' %(80/((maxStep-startStep)/sizeStep), name[k], name[k]))

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/sideAvg_*.png')