from MITgcmutils import mds
from matplotlib import pyplot as plt
import numpy as np
import os
import sys
import cmocean
import fileinput
import argparse

parser = argparse.ArgumentParser(description='Plot dynamics at zDepth')
parser.add_argument('-z','--zDepth', nargs=1, type=float,default = None,
                    help='optional depth location [m]')
parser.add_argument('-q','--quick', action='count', default=0,
                    help='quick option for last frame only, double to show plot(s)')
parser.add_argument('-k','--kValues', nargs='*', type=str, default = ['Eta','U','V'],
                    help='optional specification of views to plot [default = Eta,U,V]')
parser.add_argument('-n','--numFrames', nargs='?', type=int, default = 60,
                    help='optional specification of numFrames [default = 60]')
parser.add_argument('-t','--timeRange', nargs=2, type=int, default = None,
                    help='optional specification of start and endtime in DAYS [default = Full Range]')
parser.add_argument('--quiver', nargs="?", type=int, default = 1,
                    help='specify if quiver arrows are shown [Default = True]')
parser.add_argument('--plotX', nargs="?", type=float,default = None,
                    help='optional max X value for plotting [m]')
args = parser.parse_args()

# Pick cross section to view from file or default
yCrossSection = 1000
xCrossSection = 5000
zDepth = -50
plotDPI = 100
cleanPNGs = True
usePcolor = False
showQuiver = args.quiver
showZeros = True
flatZ = False #plotting is different if only 1 slab ie 2D

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

if(args.zDepth != None):
    print('** Manual zDepth detected **')
    zDepth = args.zDepth

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
    if f'{args.kValues[0]}.0' in file and not('D' in file): #U has issue with other files, so D is not allowed. No states have D in them it seems.
        words = file.split(".")
        # print(words)
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




#Clean up old gifs and pngs
os.system('rm -f figs/map*.png')
# os.system('rm -f figs/autoMap*.gif')

y = mds.rdmds("results/YC")/1e3
# x = mds.rdmds("results/XC")-8000
x = mds.rdmds("results/XC")/1e3
z = mds.rdmds("results/RC")

if(np.shape(z) == (1,1)):
    flatZ = True

if(os.path.isfile('input/bathymetry.bin')):
    topo = np.fromfile('input/bathymetry.bin', dtype='>f8')
    topo = topo.reshape(np.shape(x))
else:
    topo = np.zeros(np.shape(x))

if(not flatZ):
    zSlice = np.argmin(np.abs(z[:,0,0]- zDepth))
    print('depth is z =', z[zSlice,0,0], 'index', zSlice)
else:
    print(f'depth is single layer z = {-2*z[0,0]}')


if(args.quick > 0):
    startStep = maxStep
    cleanPNGs = False
print('startStep,sizeStep,maxStep:',startStep,sizeStep,maxStep)
kList = range(len(args.kValues))
for k in kList:
    print("\t" + args.kValues[k])
    for i in np.arange(startStep, maxStep + 1, sizeStep):
        if(args.kValues[k] =='SPD'):
            dataU = mds.rdmds(f"results/U", i)
            dataV = mds.rdmds(f"results/V", i)
            data = np.sqrt(dataU**2 + dataV**2)
        else:
            data = mds.rdmds(f"results/{args.kValues[k]}", i)
        if(showQuiver):
            dataU = mds.rdmds(f"results/U", i)
            dataV = mds.rdmds(f"results/V", i)
        if args.kValues[k] =='Eta':
            lvl = np.linspace(-.1,.1,31)
            cm = 'cmo.curl'
            cbarLabel = "Surf Anom [m]"
        elif args.kValues[k] =='U':
            lvl = np.linspace(-0.2,0.2,31)
            cm = 'cmo.balance'
            cbarLabel = "Speed [m/s]"
        elif args.kValues[k] =='V':
            lvl = np.linspace(-0.2,0.2,31)
            cm = 'cmo.balance'
            cbarLabel = "Speed [m/s]"
        elif args.kValues[k] =='W':
            lvl = np.linspace(-0.001,0.001,31)
            cm = 'cmo.curl'
            cbarLabel = "Vertical Speed [m/s]"
        elif args.kValues[k] =='T':
            lvl = np.linspace(0,25,31)
            cm = 'cmo.thermal'
            cbarLabel = "Temperature [C]"
        elif args.kValues[k] =='S':
            lvl = np.linspace(20,35,31)
            cm = 'cmo.haline'
            cbarLabel = "Salt [PSU]"
        elif args.kValues[k] =='SPD':
            lvl = np.linspace(0,0.2,31)
            cm = 'cmo.speed'
            cbarLabel = "Speed [m/s]"
        # if(name[k] == "SPD"):
        #     dataPlot = np.sqrt(data[2, zSlice, :, :]**2 + data[3, zSlice, :, :]**2 + data[4, zSlice, :, :]**2)
        # else:
        if(args.kValues[k] =='Eta' or flatZ):
            dataPlot = data[:, :]
        else:
            dataPlot = data[zSlice, :, :]
        if(showQuiver):
            if(flatZ):
                dataUPlot = dataU
                dataVPlot = dataV
            else:
                dataUPlot = dataU[zSlice, :, :]
                dataVPlot = dataV[zSlice, :, :]
        plt.figure(figsize=(5, 4))
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
        cbar = plt.colorbar(cp,orientation="vertical",fraction=0.06,format='%.2g')
        cbar.set_label(cbarLabel)
        ax1 = plt.gca()
        ax1.set_aspect('equal')
        plt.xlabel('Easting [km] %.3f %.3f nan: %i' %(np.nanmin(dataPlot),np.nanmax(dataPlot),np.max(np.isnan(dataPlot))))
        plt.ylabel('Northing [km]')
        if(flatZ):
            plt.title("%s full depth at %.02f days" % (args.kValues[k] ,i/86400.0*dt))
        else:
            plt.title("%s depth %.2f m at %.02f days" % (args.kValues[k], z[zSlice,0,0] ,i/86400.0*dt))
        plt.contour(x,
                    y,
                    topo,
                    [zDepth],
                    colors = 'black')
        if(showQuiver):
            downSampleX = 5
            downSampleY = 3
            u = np.squeeze(dataUPlot[::downSampleY,::downSampleX])
            v = np.squeeze(dataVPlot[::downSampleY,::downSampleX])
            plt.quiver(
                np.squeeze(x)[::downSampleY,::downSampleX],
                np.squeeze(y)[::downSampleY,::downSampleX],
                u/np.sqrt(u**2 + v**2 + 1e-12),
                v/np.sqrt(u**2 + v**2 + 1e-12),
                alpha=.2,
                # width = .1, #width of line
                # scale = 120
                )
        if(showZeros):
            cc = plt.contour(
                np.squeeze(x),
                np.squeeze(y),
                np.squeeze(dataPlot),
                [0],
                colors='gray',
                linewidths=0.5
            )
            plt.clabel(cc, inline=3, fontsize=8)
        # if(localBergs):
        #     cp2 = plt.contourf(np.squeeze(x),
        #         np.squeeze(y),
        #         np.squeeze(openFrac[zSlice, :, :]),
        #         [.4,.6,.8,.9,.95],
        #         extend="min",
        #         alpha=.1,
        #         cmap='cmo.gray')
            #cbar2 = plt.colorbar(cp2)
            #cbar2.set_label('Ocean Fraction')

        # plt.xlim([-8000, 25000]) # if zooming into a specific region
        j = i/sizeStep + startStep
        if(startStep == maxStep):
            snapShotNum = 1
        else:
            snapShotNum = (maxStep-startStep)/sizeStep
        str = "figs/map%s%05i.png" % (args.kValues[k],j)
        if(args.plotX != None):
            plt.xlim([0,args.plotX])        
        plt.tight_layout()
        plt.savefig(str, format='png', dpi=plotDPI)
        if(args.quick > 1):
            plt.show()
        plt.close()
    if(args.quick == 0):
        depthStr = ''
        timeStr = ''
        if(args.zDepth != None):
            depthStr = f'{int(np.abs(args.zDepth[0]))}'
        if(args.timeRange != None):
            timeStr = f'_T_{args.timeRange[0]}_{args.timeRange[1]}'
        os.system(f'magick -delay {500/(snapShotNum)} figs/map{args.kValues[k]}*.png -colors 256 -depth 256 figs/autoMap{depthStr}{args.kValues[k]}{timeStr}.gif')

#Clean up intermediate pngs
    if(cleanPNGs):
        os.system('rm -f figs/map*.png')
