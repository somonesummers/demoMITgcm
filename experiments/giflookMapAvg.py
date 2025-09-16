from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics depth averaged')
parser.add_argument('-z','--zDepth', nargs='?', type=float,default = None,
                    help='optional depth limit [m]. Pos values average above, neg avg below')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last frame only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=int, default = None,
                    help='optional specification of views to plot [default = all] [Temp, Sal, U, W, V, BRGmltRt,TracePlume,TraceBerg,SPD]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('--quiver', nargs="?", type=int, default = 1,
                    help='specify if quiver arrows are shown [Default = True]')
args = parser.parse_args()

# Pick cross section to view from file or default
zDepth = -50
plotDPI = 125
cleanPNGs = True
usePcolor = True
showQuiver = args.quiver
showZeros = True

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
print('Plot DPI:',plotDPI,'; clean PNGs:',cleanPNGs, '; usePcolor:', usePcolor)

fullDepth = True
if(args.zDepth != None):
    print('** Manual zDepth detected **')
    zDepth = args.zDepth
    fullDepth = False

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

#Clean up old gifs and pngs
os.system('rm -f figs/mapAvg*.png')
# os.system('rm -f figs/autoMap*.gif')

y = mds.rdmds("results/YC")
# x = mds.rdmds("results/XC")-8000
x = mds.rdmds("results/XC")
z = mds.rdmds("results/RC")

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))


if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(isBerg):
    # contourf plot
    openFrac = np.fromfile('input/openFrac.bin', dtype='>f8')
    openFrac = openFrac.reshape((np.shape(z)[0], np.shape(x)[0], np.shape(x)[1]))
    bergMask = np.fromfile('input/bergMask.bin', dtype='>f8')
    bergMask = bergMask.reshape(np.shape(x))
    for j in range(np.shape(x)[0]): #clean up non-berg parts of this mask
            for i in range(np.shape(x)[1]):
                if bergMask[j,i] == 0:
                    openFrac[:,j,i] = 1

zSlice = None
if(fullDepth == False):
    avgAbove = True
    if(zDepth < 0):
        avgAbove = False
    zDepth = np.abs(zDepth)
    zSlice = np.argmin(np.abs(z[:,0,0] + zDepth))
    print('depth is z =', z[zSlice,0,0], 'index', zSlice)

if(isBerg):
    dynName = ['dynDiag', 'dynDiag', 'dynDiag', 'dynDiag', 'dynDiag','BRGFlx','ptraceDiag','ptraceDiag','dynDiag']
    name = ["Temp", "Sal", "U", "W", "V", "BRGmltRt",'TracePlume','TraceBerg','SPD']
    cbarLabel = ["[C]", "[ppt]", "[m/s]", "[m/s]", "[m/s]", "[m/d]","[Vol Frac]","[Vol Frac]","[m/s]"]
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
        if(isBerg and os.path.isfile('results/BRGFlx.%010i.001.001.data' % i)):
            localBergs = False
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
        elif k == 8:
            lvl = np.linspace(0,np.max(uRange),31)
            cm = 'cmo.speed'

        if(name[k] == "SPD"):
            dataTmp = np.sqrt(data[2, :, :, :]**2 + data[3, :, :, :]**2 + data[4, :, :, :]**2)
        else:
            dataTmp = data[kk, :, :, :]
        if(fullDepth):
            dataPlot = np.nanmean(dataTmp,axis=0)
        elif(avgAbove):
            dataPlot = np.nanmean(dataTmp[:zSlice,:,:],axis=0)
        else:
            dataPlot = np.nanmean(dataTmp[zSlice:,:,:],axis=0)
        plt.figure(figsize=(10, 4))
        if(usePcolor):
            cp = plt.pcolormesh(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                cmap=cm,
                vmin=np.min(lvl),
                vmax=np.max(lvl),
            )
        else:
            cp = plt.contourf(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                lvl,
                extend="both",
                cmap=cm,
            )
        cbar = plt.colorbar(cp,orientation="horizontal",fraction=0.06,format='%.2f')
        cbar.set_label(cbarLabel[k])
        ax1 = plt.gca()
        ax1.set_aspect('equal')
        plt.xlabel('Along Fjord [m] %.3f %.3f nan: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.max(np.isnan(dataPlot))))
        plt.ylabel('Across Fjord [m]')
        if(fullDepth):
            plt.title("%s depth averaged at %.02f days" % (name[k] ,i/86400.0*dt))
        elif(avgAbove):
            plt.title("%s depth averaged above %.02f at %.02f days" % (name[k], z[zSlice,0,0],i/86400.0*dt))
        else:
            plt.title("%s depth averaged below %.02f at %.02f days" % (name[k], z[zSlice,0,0],i/86400.0*dt))
        plt.contour(x,
                    y,
                    topo,
                    [zDepth],
                    colors = 'black')
        if(showQuiver):
            downSampleX = 5
            downSampleY = 3
            if(fullDepth):
                u = np.nanmean((dataQuiv[2, :, :, :]),axis=0)[::downSampleY,::downSampleX]
                v = np.nanmean((dataQuiv[4, :, :, :]),axis=0)[::downSampleY,::downSampleX]
            elif(avgAbove):
                u = np.nanmean((dataQuiv[2, :zSlice, :, :]),axis=0)[::downSampleY,::downSampleX]
                v = np.nanmean((dataQuiv[4, :zSlice, :, :]),axis=0)[::downSampleY,::downSampleX]
            else:
                u = np.nanmean((dataQuiv[2, zSlice:, :, :]),axis=0)[::downSampleY,::downSampleX]
                v = np.nanmean((dataQuiv[4, zSlice:, :, :]),axis=0)[::downSampleY,::downSampleX]
            plt.quiver(
                np.squeeze(x)[::downSampleY,::downSampleX],
                np.squeeze(y)[::downSampleY,::downSampleX],
                u/np.sqrt(u**2 + v**2 + 1e-12),
                v/np.sqrt(u**2 + v**2 + 1e-12),
                alpha=.2,
                width = .001, #width of line
                scale = 120
                )
        if(showZeros and (k == 2 or k == 3 or k == 4)):
            cc = plt.contour(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                [0],
                colors='gray',
                linewidths=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        if(localBergs):
            cp2 = plt.contourf(np.squeeze(x),
                np.squeeze(y),
                np.squeeze(openFrac[zSlice, :, :]),
                [.4,.6,.8,.9,.95],
                extend="min",
                alpha=.1,
                cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')

        # plt.xlim([-8000, 25000]) # if zooming into a specific region
        j = i/sizeStep + startStep
        str = "figs/mapAvg%s%05i.png" % (name[k],j)
        # plt.xlim([0,1000])        
        plt.tight_layout()
        plt.savefig(str, format='png', dpi=plotDPI)
        if(args.quick > 1):
            plt.show()
        plt.close()
    if(args.quick == 0):
        depthStr = ''
        timeStr = ''
        if(args.zDepth != None):
            depthStr = f'{int(np.abs(args.zDepth))}'
        if(args.timeRange != None):
            timeStr = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
        os.system(f'magick -delay {500/((maxStep-startStep)/sizeStep)} figs/mapAvg{name[k]}*.png -colors 256 -depth 256 figs/autoMapAvg{depthStr}{name[k]}{timeStr}.gif')

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/map*.png')
