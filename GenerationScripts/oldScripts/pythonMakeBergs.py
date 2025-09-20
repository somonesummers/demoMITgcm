#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
from matplotlib import pyplot as plt
import cmocean
from bisect import bisect_left

writeFiles = False #set true if you'd like files

def find_closest_indices(sorted_A, sorted_B):
    closest_indices = []
    for a in sorted_A:
        pos = bisect_left(sorted_B, a)  # Find position in B where a would fit
        # Compare neighbors to find the closest
        if pos == 0:
            closest_indices.append(0)
        elif pos == len(sorted_B):
            closest_indices.append(len(sorted_B) - 1)
        else:
            before = pos - 1
            after = pos
            closest_indices.append(before if abs(sorted_B[before] - a) <= abs(sorted_B[after] - a) else after)
    return closest_indices

def setUpPrint(msg):
    print(msg)

def write_bin(fname, data):
    setUpPrint(fname + " " + str(np.shape(data)))
    if(writeFiles):
        data.astype(">f8").tofile(run_config['run_dir']+'/input/'+fname)
    else:
        setUpPrint('Not saving')

#=======================================================================================
# Make Bergs, now all in python
setUpPrint('====== Making mélange =====')

hfacThreshold = .95
# Iceberg configuration =========================
iceBergDepth = 150 # max iceberg depth [meters], used for ICEBERG package
iceExtent = 18000 # [meters] of extent of ice
iceCoverage = 80 # % of ice cover in melange, stay under 90% ideally
doMelt = 1 # do we melt (0/1 = no/yes)
doBlock = 1 # do we block and drag (0/1 = no/yes)


# Cell resolution
deltaX = 500
deltaY = 500
deltaZ = 20

# Dimensions of grid
nx = 80 # 15 km long
ny = 12 # 10 km wide (plus 500 m walls)
nz = 50 # 1000 m deep

dz = np.ones(nz)*deltaZ

#makeMasks
bergMask = np.zeros([ny,nx])
driftMask = np.zeros([ny,nx])
meltMask = np.zeros([ny,nx])
barrierMask = np.zeros([ny,nx])
bergConc = np.zeros([ny,nx])
bergMaskNums = np.zeros([ny,nx])
numBergsPerCell = np.zeros([ny,nx],dtype=np.int64)

# Berg parameters
bergType = 1 # 1 = block 2 = cone (not implemented)
alpha = 1.8 # slope of inverse power law size frequency distribution
scaling = 2 # 1 = Sulak 2017 2 = Barker 2004
maxBergDepth = iceBergDepth # (m) - set to zero if 'prescribing' max iceberg width, set at top here
minBergDepth= 20 # (m)
maxBergWidth = 0 # (m) - set to zero if 'prescribing' max iceberg depth
minBergWidth = 20 # (m)

iceExtentIndex = int(np.round(iceExtent/deltaY))

# Iceberg mask
bergMask[1:-1,1:iceExtentIndex] = 1 # icebergs in inner 5 km, all oriented east-west

# Drift mask, No drift for Melange experiments, but can toggle on here if you want
# driftMask[1:-1,1:iceExtentIndex] = 1 # calculate effect of iceberg drift on melt rates 

# Melt mask, only let bergs melt in this region (make melt water, these don't change size)
meltMask[1:-1,1:iceExtentIndex] = doMelt # Allow focus on blocking effect only

# Barrier mask
barrierMask[1:-1,1:iceExtentIndex] = doBlock # make icebergs a physical barrier to water flow
# barrierMask[plume_loc,icefront] = 0 #if using ICEPLUME, do not block in locations with an active plume

# Iceberg concentration (# of each surface cell that is filled in plan view)
bergConc[1:-1,1:iceExtentIndex] = np.linspace(iceCoverage,40,(iceExtentIndex-1)) # iceberg concentration set at top
# bergConc[1:-1,1:iceExtentIndex] = iceCoverage # iceberg concentration set at top

# print(bergConc[1:-1,1:iceExtentIndex])

desiredBergArea = np.sum(bergConc/100.0*deltaX*deltaY)
bergMaskArea = np.sum(bergMask*deltaX*deltaY)
setUpPrint('Area where bergs live: ' + str(bergMaskArea) + ' m^2')
setUpPrint('Desired berg area: ' + str(desiredBergArea) + ' m^2')
setUpPrint('Ratio: ' + str(desiredBergArea/bergMaskArea*100) + '%')

if(scaling == 1): # then use Sulak17 volume-area scaling volume = 6.0*area^1.30
    # assumes volume = L*W*D and W = L/1.62 (Dowdeswell et al 1992)
    if(maxBergWidth==0):
        maxBergWidth = 0.0642449*maxBergDepth**(5/3)
        # minBergWidth = 0.0642449*minBergDepth**(5/3)
    elif(maxBergDepth==0):
        maxBergDepth = 5.19155*maxBergWidth**(5/3)
        minBergDepth = 5.19155*minBergWidth**(5/3)
elif(scaling == 2): # Then use Barker04 width-depth relationship
    # Depth = 2.91*Width^0.71
    if(maxBergWidth==0):
        maxBergWidth = (100*10**(58/71)*maxBergDepth**(100/71)) / (291*291**(29/71))
        #minBergWidth = (100*10**(58/71)*minBergDepth**(100/71)) / (291*291**(29/71))        
    elif(maxBergDepth==0):
        maxBergDepth = 2.91*maxBergWidth^0.71
        minBergDepth = 2.91*minBergWidth^0.71

numberOfBergs = 1500 #low start, immediately doubled by scheme below, so guess low, high guesses (300%+) can cause to fail
bergTopArea = 0
areaResidual = 1
# Generate the Inverse Power Law cumulative distribution function
# over the range minBergWidth-maxBergWidth with a slope of alpha.
setUpPrint('Making bergs, this can take a few loops...')
loop_count = 1

np.random.seed(2)
setUpPrint('random seed set, not really random anymore')

while(np.abs(areaResidual) > .005 ): # Create random power dist of bergs, ensure correct surface area
    numberOfBergs = round(numberOfBergs * (1 + areaResidual))  
    setUpPrint('\tnumberOfBergs: ' + str(numberOfBergs))
    x_width = np.arange(minBergWidth, maxBergWidth, (maxBergWidth-minBergWidth)/(numberOfBergs*1e2))
    x_depth = np.arange(minBergDepth, maxBergDepth, (maxBergDepth-minBergDepth)/(numberOfBergs*1e2))
    inversePowerLawPDF_width = ((alpha-1) / minBergWidth) * (x_width/minBergWidth) ** (-alpha)
    inversePowerLawPDF_depth = ((alpha-1) / minBergDepth) * (x_depth/minBergDepth) ** (-alpha)
        # Get the CDF numerically
    inversePowerLawCDF_width = np.cumsum(inversePowerLawPDF_width)
    inversePowerLawCDF_depth = np.cumsum(inversePowerLawPDF_depth)
        # Normalize
    inversePowerLawCDF_width = inversePowerLawCDF_width / inversePowerLawCDF_width[-1]
    inversePowerLawCDF_depth = inversePowerLawCDF_depth / inversePowerLawCDF_depth[-1]
        
        # Generate number_of_bergs uniformly distributed random numbers.
    uniformlyDistributedRandomNumbers = np.random.uniform(0,1,numberOfBergs)
    
    inversePowerLawDistNumbers_width = np.zeros(uniformlyDistributedRandomNumbers.size);
    inversePowerLawDistNumbers_depth = np.zeros(uniformlyDistributedRandomNumbers.size);
    nearestIndex_width = [0] * uniformlyDistributedRandomNumbers.size
    nearestIndex_depth = [0] * uniformlyDistributedRandomNumbers.size
    
    # for i in range(uniformlyDistributedRandomNumbers.size):  #this is pretty slow 
    #     nearestIndex_width[i] = np.abs(uniformlyDistributedRandomNumbers[i]-inversePowerLawCDF_width).argmin();
    #     nearestIndex_depth[i] = np.abs(uniformlyDistributedRandomNumbers[i]-inversePowerLawCDF_depth).argmin();
    # This works by leaveraging that inversPowerLaw is sorted
    nearestIndex_width = find_closest_indices(uniformlyDistributedRandomNumbers,inversePowerLawCDF_width)
    nearestIndex_depth = find_closest_indices(uniformlyDistributedRandomNumbers,inversePowerLawCDF_depth)


    inversePowerLawDistNumbers_width = x_width[nearestIndex_width];
    inversePowerLawDistNumbers_length = inversePowerLawDistNumbers_width/1.62 # Widths are bigger 
    tooWide = np.count_nonzero(inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold))  #disallow completely full cells
    tooLong = np.count_nonzero(inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold))
    inversePowerLawDistNumbers_width[inversePowerLawDistNumbers_width > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(hfacThreshold) # Max width is grid cell (assumed square)
    inversePowerLawDistNumbers_length[inversePowerLawDistNumbers_length > deltaX * np.sqrt(hfacThreshold)] = deltaX * np.sqrt(hfacThreshold) # Max length is grid cell (assumed square)
    if(tooLong + tooWide > 0):
        setUpPrint('\t\tBergs clipped: %i for width, %i for length' % (tooWide, tooLong))
    
    inversePowerLawDistNumbers_depth = x_depth[nearestIndex_depth]; #depths don't get clipped
    
    bergTopArea = sum(inversePowerLawDistNumbers_width*inversePowerLawDistNumbers_length)
    areaResidual = (desiredBergArea - bergTopArea)/desiredBergArea
    setUpPrint('\t\t%.2f %% Bergs' % (bergTopArea/bergMaskArea*100))
    setUpPrint('\t\tareaResidual %.2f %%' % (areaResidual * 100))
    loop_count += 1
setUpPrint('====== Success! Found our bergs =====')
setUpPrint('Width min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_width),np.mean(inversePowerLawDistNumbers_width),np.max(inversePowerLawDistNumbers_width)))
setUpPrint('Depth min/mean/max: %f/%f/%f [m]' % (np.min(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth),np.max(inversePowerLawDistNumbers_depth)))
setUpPrint('Total Berg Area %f' % bergTopArea)
setUpPrint('Total Berg fract: %.2f %%' % (bergTopArea/bergMaskArea*100))

# Now we sort these berg into cell, randomly 
bergMaski = 0  #Bad name, but this is the count of cells that will recieve bergs
bergDict = {}

for j in range(ny):
    for i in range(nx):
        if(bergMask[j,i] == 1):
            # print('i,j, bergmask',i,j,bergMask[j,i])
            bergMaski = 1 + bergMaski #Needs to start at 1, as non-bergs will be 0
            bergMaskNums[j,i] = bergMaski #Assign Mask Nums, not random as we'll randomly place bergs in cells
            bergDict[bergMaski] = [j,i] #This lets us do 1-D loops for the whole grid
setUpPrint('%i cells with bergs' % bergMaski)

# Sort my bergs
sorted_indices = np.argsort(-inversePowerLawDistNumbers_depth) # Sort backwards to get descending from big to small bergs 
sorted_depth = inversePowerLawDistNumbers_depth[sorted_indices]
sorted_width = inversePowerLawDistNumbers_width[sorted_indices]
sorted_length = inversePowerLawDistNumbers_length[sorted_indices]
assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) # In this script, every berg has a home

# Array for bergs
bergsPerCellLimit = 500
icebergs_depths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_widths = np.zeros([bergMaski,bergsPerCellLimit])
icebergs_length = np.zeros([bergMaski,bergsPerCellLimit])  #careful, not plural as to length match

np.random.seed(2)
assignedCell = np.random.randint(0,bergMaski,[numberOfBergs]) #every Berg has a spot

icebergs_per_cell = np.zeros([bergMaski],dtype=np.int16)
icebergs_area_per_cell = np.zeros([bergMaski])

for i in range(numberOfBergs): 
    j = assignedCell[i]
    # print('looking at mask number',j,'at berg',i)
    # print('Berg number', icebergs_per_cell[j],'in this cell')
    bergArea = sorted_width[i] * sorted_length[i]
    loopLimiter = 0
    while(bergArea > (deltaX * deltaY  * (bergConc[bergDict[j+1][0],bergDict[j+1][1]]/100) - icebergs_area_per_cell[j])): #if above 'full', pick random new cell, accept with decreasing probability
        j_old = j
        j = np.random.randint(0,bergMaski)
        loopLimiter += 1
        if((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY) < hfacThreshold - .01): #only consider accepting if under 95
            odds = np.abs(np.random.normal(0,.5,1))  #randomly accepts those that are big in overfull cells, but at decreasing frequency
            overFull = ((bergArea + icebergs_area_per_cell[j])/(deltaX * deltaY)*100 - bergConc[bergDict[j+1][0],bergDict[j+1][1]])
            if(odds > overFull):
                 # print('accepting overfull')
                 assignedCell[i] = j  #if we it a shuffling critera
                 break
        if(loopLimiter > bergMaski*20): #eventually we have to force some in
            indexesAllowed = np.where((deltaX * deltaY * hfacThreshold - icebergs_area_per_cell)  > bergArea)
            randi = np.random.randint(0,len(indexesAllowed[0]))
            j = indexesAllowed[0][randi]
            assignedCell[i] = j  #if we it a shuffling critera, must line up for calculation below
            # setUpPrint('\t Randomly missed, will force into cell with room: %i' % j)
            if((np.min(icebergs_area_per_cell) + bergArea)/(deltaX * deltaY) > hfacThreshold):
                setUpPrint('WARNING cell very full: %.2f%%' %((np.min(icebergs_area_per_cell) + bergArea)*100/(deltaX * deltaY)))
            break
     
    icebergs_depths[j,icebergs_per_cell[j]] = sorted_depth[i]
    icebergs_widths[j,icebergs_per_cell[j]] = sorted_width[i]
    icebergs_length[j,icebergs_per_cell[j]] = sorted_length[i]
    icebergs_per_cell[j] += 1
    # icebergs_area_per_cell[j] = np.sum(icebergs_widths[j,:]*icebergs_length[j,:])
    icebergs_area_per_cell[j] += bergArea
setUpPrint('Bergs per cell and filled faction at surface for spot check')     
setUpPrint(icebergs_per_cell)
setUpPrint(np.round(icebergs_area_per_cell/(deltaX*deltaY),2))
setUpPrint('Max fill is: %.2f%%' % (np.nanmax(icebergs_area_per_cell/(deltaX*deltaY))*100))

# All bergs now sorted 
openFrac = np.zeros([nz,ny,nx])
SA = np.zeros([nz,ny,nx])
SA[:,:,:] = np.nan
cellVolume = deltaX*deltaY*deltaZ

#This loop knows where all bergs are already, different from searching for all bergs across entire grid
for i in range(bergMaski):
    bergCount = icebergs_per_cell[i]
    numBergsPerCell[bergDict[i+1][0],bergDict[i+1][1]] = bergCount
    if(bergCount > 0):
        lengths = icebergs_length[i,icebergs_length[i,:] > 0] #return only non-zeros
        widths = icebergs_widths[i,icebergs_widths[i,:] > 0] #return only non-zeros
        depths = icebergs_depths[i,icebergs_depths[i,:] > 0] #return only non-zeros
        for k in range(nz):
            d_bot = k*deltaZ + deltaZ #bottom of depth bin
            d_top = k*deltaZ
            volume1 = deltaZ * lengths[depths > d_bot] * widths[depths > d_bot]
            SA1 = deltaZ*2*(lengths[depths > d_bot] + widths[depths > d_bot])
            partialFill = (depths < d_bot) & (depths > d_top)
            #partial fill
            volume2 = (depths[partialFill] - d_top) * lengths[partialFill] * widths[partialFill]
            #partial sides
            SA2 = (depths[partialFill] - d_top)*2*(lengths[partialFill] + widths[partialFill]) 
            #bottom
            SA3 = lengths[partialFill] * widths[partialFill]
            #print(np.sum(volume1), np.sum(volume2))
            openFrac[k,bergDict[i+1][0],bergDict[i+1][1]] = 1-((np.sum(volume1) + np.sum(volume2))/cellVolume)
            SA[k,bergDict[i+1][0],bergDict[i+1][1]] = np.sum(SA1) + np.sum(SA2) + np.sum(SA3)
    elif(bergCount == 0):
        openFrac[:,bergDict[i+1][0],bergDict[i+1][1]] = 1
        SA[:,bergDict[i+1][0],bergDict[i+1][1]] = 0

# Plots for reference on whats happening berg-wise
fig = plt.figure()
plt.subplot(2,2,1)
for i in range(bergMaski):
    plt.plot(openFrac[:,bergDict[i+1][0],bergDict[i+1][1]],-np.cumsum(dz),alpha=.5,color='xkcd:gray',linewidth=.5)
plt.plot(np.mean(openFrac[:,bergMask==1],1),-np.cumsum(dz),alpha=1,color='xkcd:black',linewidth=1,linestyle='--',label='Average Bergs')
plt.plot([0,1],[-maxBergDepth,-maxBergDepth],color = 'xkcd:red',linestyle=':', label='Target Max Depth')
plt.plot([1-np.max(bergConc)/100,1-np.max(bergConc)/100],[-nz*deltaZ,0],color = 'xkcd:gray',linestyle=':',label='Target Max Berg Conc')
plt.xlabel('Open Fraction of Cells')
plt.ylabel('Depth [m]')
# plt.legend()  #not quite room so off for now

plt.subplot(2,2,2)
plt.hist(inversePowerLawDistNumbers_depth,bins = 50)
plt.ylabel('Count')
plt.xlabel('Depth [m]')

plt.subplot(2,2,3)
plt.hist(inversePowerLawDistNumbers_width,bins = 50)
plt.ylabel('Count')
plt.xlabel('Width [m]')

plt.subplot(2,2,4)
plt.hist(inversePowerLawDistNumbers_length,bins = 50)
plt.ylabel('Count')
plt.xlabel('Length [m]')
fig.tight_layout()
if(writeFiles):
    plt.savefig(run_config['run_dir']+'/input/bergStatistics.png', format='png', dpi=200)
plt.show()

fig = plt.figure()
pltHelper = 1-openFrac[0,:,:]
pltHelper[bergMask == 0] = np.nan
plt.subplot(211)
pc = plt.pcolormesh(pltHelper,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Iceberg Cover and Residual')
plt.ylabel('Cell across fjord')
cbar.set_label('iceberg cover')

plt.subplot(212)
pc = plt.pcolormesh(pltHelper - bergConc/100,cmap='cmo.ice')
cbar = plt.colorbar(pc)
pc_min = np.nanmin(pltHelper) * 100
pc_max = np.nanmax(pltHelper) * 100
plt.xlabel("Cell along fjord, Ice coverage is %.2f%% - %.2f%%" % (pc_min, pc_max))
plt.ylabel('Cell across fjord')
cbar.set_label('cover resid')
fig.tight_layout()
if(writeFiles):
    plt.savefig(run_config['run_dir']+'/input/bergMap.png', format='png', dpi=200)
plt.show()

fig = plt.figure()
pc = plt.pcolormesh(meltMask,cmap='cmo.ice_r')
cbar = plt.colorbar(pc)
plt.suptitle('Melt Mask')
plt.ylabel('Cell across fjord')
plt.xlabel("Cell along fjord")
fig.tight_layout()
if(writeFiles):
    plt.savefig(run_config['run_dir']+'/input/meltMask.png', format='png', dpi=200)
plt.show()

## write iceberg txt files
# setUpPrint('Saving text files for bergs...')
# if(writeFiles):
#     for i in range(bergMaski):
#         with open(run_config['run_dir'] + '/input/iceberg_depth_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_depths[i,icebergs_depths[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))
#         with open(run_config['run_dir']+'/input/iceberg_width_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_widths[i,icebergs_widths[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))
#         with open(run_config['run_dir'] + '/input/iceberg_length_%05i.txt' % (i+1) , 'w') as file_handler:
#             for item in icebergs_length[i,icebergs_length[i,:] > 0]:
#                 file_handler.write("{}\n".format(item))

icebergs_depths2D = np.zeros([bergsPerCellLimit,ny,nx])
icebergs_widths2D = np.zeros([bergsPerCellLimit,ny,nx])
icebergs_length2D = np.zeros([bergsPerCellLimit,ny,nx])

for k in range(bergMaski):
    j = bergDict[k+1][0]
    i = bergDict[k+1][1]
    icebergs_depths2D[:,j,i] = icebergs_depths[k,:]
    icebergs_widths2D[:,j,i] = icebergs_widths[k,:]
    icebergs_length2D[:,j,i] = icebergs_length[k,:]

# write global files
write_bin('bergMask.bin',bergMask)
write_bin('bergMaskNums.bin',bergMaskNums)
write_bin('numBergsPerCell.bin',numBergsPerCell)
write_bin('openFrac.bin',openFrac)
write_bin('totalBergArea.bin',SA)
write_bin('meltMask.bin',meltMask)
write_bin('driftMask.bin',driftMask)
write_bin('barrierMask.bin',barrierMask)
write_bin('icebergs_depths.bin',icebergs_depths2D)
write_bin('icebergs_widths.bin',icebergs_widths2D)
write_bin('icebergs_length.bin',icebergs_length2D)



setUpPrint('Berg setup is done.')