# -*- coding: utf-8 -*-

"""

Created on Mon Feb 06 07:43:50 2017

@author: Jesse Trinity, Jason Coleman

Last used: 9-9-17

"""



import open_calcium_data_fast as ocdfast

import StepCodeFile as scf

import numpy as np

from matplotlib import pyplot as plt

import calcium_imaging_data_fast as cidfast

from copy import deepcopy

from collections import OrderedDict

import matplotlib.image as mpimg

import readroi as roizip

import dill    


#Gets Step code file
#fileDirectoryTS = 'F:/coleman lab/jasonc/thygcamp6s_test2/'
fileDirectoryTS = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_LT4/'
#fileDirectory = 'F:/coleman lab/jasonc/thygcamp6s_test2/'
fileDirectory = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_LT4/'

#fileBin = 'mThy6s2_alldrift_D2_001_12_data.bin'
#fileBin = 'mThy6s_001_D3_6_data.bin'
fileBin = 'mThy6s_001_D3Z1t2_13_data.bin'

#imgname = 'STD_mThy6s2_alldrift_D2_001.tif'
#imgname = 'STD_mThy6s2_alldrift_D3_001Z1.tif'
imgname = 'STD_mThy6s2_alldrift_D3_001Z1t2.tif'

#roizipname = 'mThy6s2_alldrift_D2_001_ROI.zip'
#roizipname = 'mThy6s2_alldrift_D3_001Z1_ROI.zip'
roizipname = 'mThy6s2_alldrift_D3_001Z1t2_ROI.zip'

filenamepkl = roizipname.replace('ROI.zip','SESSION.pkl')

daqRecordingFreq = 1000.0 #sampling rate - Hz
stimList = [0,45,90,135,180,225,270,315]
csvfiletype = 3 # 1= ; 2= ; 3=FIJI/ImageJ csv intensity file

datatype = 'raw' # use raw IntDen/Area data, for normalizing to visual stim
# datatype = 'filtered' # use EWMA smoothing, etc
# Parameters for Konnerth lab dF/F; EWMA
#for a 30Hz (Fs) imaging systems the following parameter setup is
#    recommended (empirical note on Nature paper): 
#    t_0= 0.2;
#    t_1=0.75;
#    t_2=3;
t_0 = 0.2 #0.2 for 30fps; 0.8 for 7.6fps
t_1 = 0.75
t_2 = 3
mpmSamplingFreq = 30 #~fps for 2p

# gray time (in between stims)
gray_offset= 7.0 #seconds



# Begin decoding and data extraction

code_file = scf.StepCodeFile(fileDirectoryTS+fileBin,8,8,1000)

stims = code_file.get_stim_angles(stimList)

Bframes = [ts[0] for ts in code_file.timestamps[2]]
"""
recording_freq = 1000 # samples/sec

numberOfFrames = 3250

frame_rate= 7.6

Bframes = [ts[0] for ts in code_file.timestamps[2]]

Bframe_start = 23289

Bframe_length = (numberOfFrames/frame_rate) * recording_freq

timestamp = code_file.data[0][Bframe_start:]

step_data = code_file.data[1][Bframe_start:]

#Bframes = np.arange(0, Bframe_length, recording_freq*(1.0/frame_rate))
#make Bframe
"""

#get intensity file

cellsraw, cellsmean, areas, xycoords, Data = ocdfast.openData(csvfiletype)

# calculate DeltaF/F

# Run Konnerth lab deltaF/F and EWMA filter - USAGE: (cells = raw intensity data)                      

if datatype == 'filtered':
    dff, dff_ewma, dff_offset = cidfast.run_deltaf_ewma(cellsmean, t_0, t_1, t_2, mpmSamplingFreq)
    cells = deepcopy(dff_ewma)
elif datatype == 'raw':
    cells = deepcopy(cellsmean)
    
handler = cidfast.FileHandler()

#makes dictionary of [cellnumber][orientation][block]

gray_offset *= daqRecordingFreq

response_data = OrderedDict()

grays = OrderedDict()

response_indices = OrderedDict()

grays_indices = OrderedDict()

for cell in range(len(cells)):

    response_data[cell] = OrderedDict()

    grays[cell] = OrderedDict()

    response_indices[cell] = OrderedDict()

    grays_indices[cell] = OrderedDict()

    for stim in stims:

        response_data[cell][stim] = list()

        grays[cell][stim] = list()

        response_indices[cell][stim] = list()

        grays_indices[cell][stim] = list()

        for ts in stims[stim]:

            begin = float(ts[0])

            end = float(ts[1])

            begin_frame_time = handler.get_nearest(begin, Bframes)

            begin_gray_time = handler.get_nearest(begin - gray_offset, Bframes)

            

            end_frame_time = handler.get_nearest(end, Bframes)

            end_gray_time = handler.get_nearest(begin, Bframes)

            

            begin_frame = list(Bframes).index(begin_frame_time)

            begin_gray = list(Bframes).index(begin_gray_time)

            

            end_frame = list(Bframes).index(end_frame_time)

            end_gray = list(Bframes).index(end_gray_time)

            

            chunk = cells[cell][int(begin_frame):int(end_frame)]

            gray_chunk = cells[cell][int(begin_gray):int(end_gray)]

            

            response_data[cell][stim].append(chunk)

            grays[cell][stim].append(gray_chunk)

            

            response_indices[cell][stim].append([int(begin_frame),int(end_frame)])

            grays_indices[cell][stim].append([int(begin_gray),int(end_gray)])



#example: plots all 5 block of degree 45 orientation for cell 0

#cell_0_45 = response_data[0][45]

#for block in cell_0_45:

#    plt.plot(block)



#Gets response averages

response_avgs = OrderedDict()

for cell in response_data:
    response_avgs[cell] = OrderedDict()

    for orientation in response_data[cell]:

        cell_ori = response_data[cell][orientation]
        
        #trim signals down to shortest length for averaging
        min_chunk_length = min([len(x) for x in cell_ori])
        cell_ori = [x[0:min_chunk_length] for x in cell_ori]

        A = np.array(cell_ori)

        B = np.mean(A, axis = 0)

        response_avgs[cell][orientation] = B


#example: plots all 45 degree responses

#for cell in response_avgs:

#    plt.plot(response_avgs[cell][45])

#example: plots all orientations for cell 10

#for ori in response_avgs[10]:

#    plt.plot(response_avgs[0][ori])
    
# stim_onsetframes, grays_onsetframes = getOnsetIndices(response_indices, grays_indices)
def getOnsetIndices(indices1, indices2):
    
    tempstims = indices1[0]
    
    tempgrays = indices2[0]
    
    indices1_onsetframes = list()
    
    indices2_onsetframes = list()
    
    for jj in tempstims:
    
        for i in range(len(tempstims[jj])):
    
            indices1_onsetframes.append([tempstims[jj][i][0], tempstims[jj][i][1]])
    
    for kk in tempgrays:
    
        for i in range(len(tempgrays[kk])):
    
            indices2_onsetframes.append([tempgrays[kk][i][0], tempgrays[kk][i][1]])
    
    return indices1_onsetframes, indices2_onsetframes



# Check DeltaF/F
# Plot all traces (whole)
# tracePlot(dff_ewma,)
def tracePlot(data):
    
    stim_onsetframes, grays_onsetframes = getOnsetIndices(response_indices, grays_indices)
    
    plt.subplots()

    rangeX = len(data)
    

    for i in range(rangeX):

        plt.plot(data[i]+i)
        
        
    # plot estimated stim onsets
    for i in range(len(stim_onsetframes)):
    
        plt.axvline(stim_onsetframes[i][0], color = 'g')
    
        plt.axvline(stim_onsetframes[i][1], color = 'k')
    
        #plt.axvline(grays_onsetframes[i][0], color = 'b')
    
        #plt.axvline(grays_onsetframes[i][1], color = 'r')


# Plot all orientation response means or pre-stim gray, etc. on one plot
def plotROIavgs(data,roisOfInterest,ylimit):

    for i in range(len(roisOfInterest)):
        roinumber = roisOfInterest[i]
        plt.subplots()
        
        for ori in data[roinumber]:
            if ori == 315:
                plt.plot(data[roinumber][ori], '#ff8c00')
            elif ori != 315:
                plt.plot(data[roinumber][ori])
            
        plt.ylim([0, ylimit])
        
        plt.title(filenamepkl.replace('_SESSION.pkl','')+' roi: ' + str(roinumber))
        plt.legend(['0', '45', '90', '135', '180', '225', '270', '315'])


# Just trying to get a picture of possible "selectivity" for each cell - idea is to have 5 trials plotted for each ori, 
# then overlay bold, black avg trace
def multiPlot(data,avgdata,nrows, ncols, celllist, yaxlimit):

    fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all')  

    cols = ['0', '45', '90', '135', '180', '225', '270', '315']

    rows = ['{}'.format(row) for row in celllist]


    for ax, col in zip(axes[0], cols):

        ax.set_title(col)


    for ax, row in zip(axes[:,0], rows):

        ax.set_ylabel(row, rotation=0, size='large')  


    for cellnum in range(nrows):

        temp=data[celllist[cellnum]]

        tempavg=avgdata[celllist[cellnum]]

        templist = list()

        templistavg = list()


        for ori in temp:

            templist.append(temp[ori])

            templistavg.append(tempavg[ori])

        for j in range(ncols):

            #for k in range(len(templist[j])):

            #    axes[cellnum,j].plot(templist[j][k]) 

            axes[cellnum,j].plot(templistavg[j], linewidth=1.5, color='k')

            axes[cellnum,j].set_ylim([0, yaxlimit])
            


def plotStack(imgdir, imgname, roizipname):

    img = mpimg.imread(imgdir+imgname)

    plt.subplots()

    imgplot = plt.imshow(img)

    imgplot = plt.imshow(img, cmap="gray")#, origin='lower')

    #read ROI zip file from Fiji/ImageJ (magic wand objects)

    #you will need readroi.py file (on GitHub)

    a=roizip.read_roi_zip(imgdir + roizipname)
    

    for j in range(len(a)):
        
        #print a[j]
        #print type(a[j])
        #mac
        #ylist = [a[j][1][i][0] for i in range(len(a[j][1]))]
        
        #windows
        ylist = [a[j][i][0] for i in range(len(a[j]))]
        
        ylist.append(ylist[0])
        
        #mac
        #xlist = [a[j][1][i][1] for i in range(len(a[j][1]))]
        
        #windows
        xlist = [a[j][i][1] for i in range(len(a[j]))]

        xlist.append(xlist[0])

        plt.plot(xlist, ylist, linestyle = '-', linewidth=0.5)

        plt.annotate(str(j+1), xy=(1, 1), xytext=(xlist[i], ylist[i]), color='limegreen', fontsize=8)
        

def get_colormaps(delta_f, keys, minimum, maximum):
    delta_f = np.array(delta_f).astype(np.float)
    fig, ax = plt.subplots()
    heatmap = ax.pcolor(delta_f, vmin = minimum, vmax = maximum)
    ax.set_yticks(np.arange(len(delta_f))+0.5, minor=False)
    ax.set_yticklabels(keys, minor=False, fontsize = 8)
    ax.invert_yaxis()
    plt.colorbar(heatmap)
    plt.show()
    
            
def plotHeatAvgs(data, datakeys, stimonset, minheat, maxheat, plottitle):
    get_colormaps(data, datakeys, minheat, maxheat)
    plt.title(plottitle)
    plt.axvline(x = stimonset, color = 'green', ls = 'dashed')
    plt.axis('tight')  
    
  
# Run some functions
tracePlot(cells)
plotStack(fileDirectory, imgname, roizipname)
multiPlot(response_data,response_avgs,32, 8,range(0,32),2.0)



"""
# Save all variables from the session; to load --> >>> dill.load_sessio(fileDirectory+filenamepkl) 
dill.dump_session(fileDirectory+filenamepkl)


#%% Work out how to divide response_avgs by 1s of gray        
#Gets response averages

# plot gray from [roi6][ori0][session1]
#plt.plot(grays[6][0][1])

gray_avgs = OrderedDict()

for cell in grays:
    gray_avgs[cell] = OrderedDict()

    for pregray in grays[cell]:

        cell_gray = grays[cell][pregray]
        
        #trim signals down to shortest length for averaging (1sec or 30frames for gray)
        min_chunk_length = min([len(x) for x in cell_gray])
        cell_gray = [x[0:min_chunk_length] for x in cell_gray]

        C = np.array(cell_gray)

        D = np.mean(C, axis = 0)

        gray_avgs[cell][pregray] = D
        
#%% Plot each ROI orientation means (dff_ewma)

plotROIavgs(response_avgs,[6,7,9,10,11,13,14,15,20,25,26],1000.0)
#plotROIavgs(gray_avgs,[6,7,9,10,11,13,14,15,20,25,26],2.0)

"""