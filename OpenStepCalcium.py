# -*- coding: utf-8 -*-

"""

Created on Mon Feb 06 07:43:50 2017

@author: Jesse Trinity, Jason Coleman

Last used: 9-9-17

"""

import numpy as np

from scipy import signal

from scipy.signal import butter, lfilter, freqz

import scipy.stats as stats

from matplotlib import pyplot as plt

import matplotlib.image as mpimg

import readroi as roizip

from collections import OrderedDict

import pickle


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


def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    
    return y
    

def get_ori_responses(data):
    """
    Calulate mean of each trace avg for stats
    data[cell][orientation]
    """        
    a = OrderedDict()        
    for cell in data:
        a[cell] = OrderedDict()
        for ori in data[cell]:
            a[cell][ori] = np.mean(data[cell][ori])
            
    return a
   
    
# Check DeltaF/F
# Plot all traces (whole)
# tracePlot(dff_ewma,)
def tracePlot(data, response_indices, grays_indices):
    
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
def plotROIavgs(data,filename,roisOfInterest,ylimit):

    for i in range(len(roisOfInterest)):
        roinumber = roisOfInterest[i]
        plt.subplots()
        
        for ori in data[roinumber]:
            if ori == 315:
                plt.plot(data[roinumber][ori], '#ff8c00')
            elif ori != 315:
                plt.plot(data[roinumber][ori])
            
        plt.ylim(ylimit)
        
        plt.title(filename +' roi: ' + str(roinumber))
        plt.legend(data[0].keys())
        #plt.legend(['0', '45', '90', '135', '180', '225', '270', '315'])
        plt.axvline(x = 30, color = 'black', ls = 'dashed')


# Just trying to get a picture of possible "selectivity" for each cell - idea is to have 5 trials plotted for each ori, 
# then overlay bold, black avg trace
def multiPlot(avgdata,nrows, ncols, celllist, yaxlimit, window, responseindex):

    fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all')  

    cols = avgdata[0].keys()

    rows = ['{}'.format(row) for row in celllist]


    for ax, col in zip(axes[0], cols):

        ax.set_title(col)


    for ax, row in zip(axes[:,0], rows):

        ax.set_ylabel(row, rotation=0, size='large')  


    for cellnum in range(nrows):

        tempavg=avgdata[cellnum]
        
        axisID = -1 # count which column
            
        for ori in tempavg:
            
            axisID += 1

            axes[cellnum,axisID].plot(tempavg.values()[axisID], linewidth=1.5, color='b')
            
            if responseindex[cellnum][ori][0] > 0:
                tempresponsekey = 'co'
            elif responseindex[cellnum][ori][0] < 0:
                tempresponsekey = 'ro'
            elif responseindex[cellnum][ori][0] == 0:
                tempresponsekey = 'wo'
                          
            axes[cellnum,axisID].plot(110.0, 1.0, tempresponsekey, markersize = 10)

            axes[cellnum,axisID].set_ylim(yaxlimit)
            
            axes[cellnum,axisID].axvline(x = window, color = 'limegreen')
            #axes[cellnum,ori].axvline(x = window, color = 'black', ls = 'dashed')
            


def plotStack(imgdir, imgname, roizipname, responses_i):

    img = mpimg.imread(imgdir+imgname)

    plt.subplots()

    imgplot = plt.imshow(img)

    imgplot = plt.imshow(img, cmap="gray")#, origin='lower')

    #read ROI zip file from Fiji/ImageJ (magic wand objects)

    #you will need readroi.py file (on GitHub)

    a=roizip.read_roi_zip(imgdir + roizipname)
    
    # find repsonsive cell (binary)
    c=list()
    for cell in responses_i:
        b = list()
        for ori in responses_i[cell]:
            if responses_i[cell][ori][0] == 1:
                b.append(1)
            elif responses_i[cell][ori][0] != 1:
                b.append(0)
        if sum(b) > 0:
            print cell
            c.append(1)
        elif sum(b) == 0 or sum(b) < 0:
            c.append(0)
    

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
        
        if c[j] > 0: # if cell is responsive (see loop above)

            plt.annotate(str(j), xy=(1, 1), xytext=(xlist[i], ylist[i]), color='limegreen', fontsize=8)
        
        elif c[j] <= 0: # if cell is NOT responsive (see loop above)

            plt.annotate(str(j), xy=(1, 1), xytext=(xlist[i], ylist[i]), color='red', fontsize=8)
        

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
    
    
def get_responseClass(response_avgs, pregray1s_response_avgs, pre_response_post_avgs, stimwindow, pthresh, togglePlot):
    response_indices = dict()
    for cell in response_avgs:
        tempstim = response_avgs[cell]
        tempgray = pregray1s_response_avgs[cell]
        response_indices[cell] = dict()
        
        #test by orientation
        #ori = 0
        for ori in tempstim:
            response_indices[cell][ori] = list()            
            
            fstat,pval = stats.mannwhitneyu(tempstim[ori][0:stimwindow], tempgray[ori])
            
            if (pval < pthresh and np.mean(tempstim[ori][0:stimwindow]) > 0.1):
                print "Response in cell "+str(cell)+" for "+str(ori)+"deg"+" p="+str(pval)
                response_indices[cell][ori].append(1.0)
            elif (pval < pthresh and np.mean(tempstim[ori][0:stimwindow]) < -0.1):
                print "Depression in cell "+str(cell)+" for "+str(ori)+"deg"+" p="+str(pval)
                response_indices[cell][ori].append(-1.0)
            elif ((pval > pthresh or pval < pthresh) and np.mean(tempstim[ori][0:stimwindow]) <= 0.1):
                print "NO response in cell "+str(cell)+" for "+str(ori)+"deg"
                response_indices[cell][ori].append(0)
            elif ((pval > pthresh or pval < pthresh) and np.mean(tempstim[ori][0:stimwindow]) >= -0.1):
                print "NO response in cell "+str(cell)+" for "+str(ori)+"deg"
                response_indices[cell][ori].append(0)
                
            if togglePlot == 1:
                plt.subplots()
                plt.plot(pre_response_post_avgs[cell][ori])
                #plt.legend(str(ori))
                plt.plot(30, np.mean(tempgray[ori]), 'ko')
                plt.plot(100, np.mean(tempstim[ori][0:stimwindow]), 'bo')
                if (pval < pthresh and np.mean(tempstim[ori][0:stimwindow]) > 0.1):
                    plt.title('YES+|Cell#'+str(cell)+'; ori='+str(ori)+'deg'+'; p='+str(pval))
                elif (pval < pthresh and np.mean(tempstim[ori][0:stimwindow]) < -0.1):
                    plt.title('YES-|Cell#'+str(cell)+'; ori='+str(ori)+'deg'+'; p='+str(pval))
                else:
                    plt.title('NO|Cell#'+str(cell)+'; ori='+str(ori)+'deg'+'; p='+str(pval))
    
    return response_indices
    

def load_responses(filedir, picklefilename):

    # D2 Z1
    with open(filedir+picklefilename) as f:  # Python 3: open(..., 'rb')
        all = pickle.load(f)
    a = all['responses_means']
    b = all['all_response_indices']
    
    return a,b
            
    # now use to classify each cell - 1) responsive? 2) preference?
    #   3) potentiated or depressed?


def load_spontaneous(filedir, picklefilename):
    
    """
    loads data from the first X sec prior to first stimulus onset (gray screen)
    """

    # D2 Z1
    with open(filedir+picklefilename) as f:  # Python 3: open(..., 'rb')
        all = pickle.load(f)
    a = all['spontaneous_raw']
    b = all['spontaneous_dff']
    
    return a,b
    

def plot_time_responses(t1_data, t1_indices, t2_data, t2_indices, plotall):
    
    plt.subplots()
    
    #colors = plt.cm.rainbow(np.linspace(0, 1, len(t1_data)))
    colors = plt.cm.jet(np.linspace(0, 1, len(t1_data)))
    C=-1
    

    t1_output_all = list()
    t2_output_all = list()    
    
    for cell in t1_data:
        
        C += 1
        
        #plt.subplots()
        
        for ori in t1_data[cell]:
            
            response1 = t1_indices[cell][ori]
            response2 = t2_indices[cell][ori]
            
            if plotall == 0:
                
                if response1[0] + response2[0] == 2.0:
                    
                    t1_output_all.append(t1_data[cell][ori])
                    t2_output_all.append(t2_data[cell][ori])                    
                    
                    plt.plot(t1_data[cell][ori], t2_data[cell][ori], marker = 'o', color = colors[C])
                    
                elif response1[0] + response2[0] != 2.0:
                    #plt.plot(t1_data[cell][ori], t2_data[cell][ori], 'kx')
                    # to preserve color colding by all of a cell's orientation responses
                
                    t1_output_all.append('NaN') # Is there a better equivalent to Matlab 'NaN' (ie not strings)?
                    t2_output_all.append('NaN')                
                
                    #plt.plot(0, 0, marker = 'o', color = colors[C])
                
            elif plotall == 1:
                
                t1_output_all.append(t1_data[cell][ori])
                t2_output_all.append(t2_data[cell][ori])

                plt.plot(t1_data[cell][ori], t2_data[cell][ori], marker = 'o', color = colors[C])
            
            plt.title('t1 v t2')
            plt.xlabel('time1')
            plt.ylabel('time2')


    plt.ylim([0, 1.0])
    plt.xlim([0, 1.0])
    
    return t1_output_all, t2_output_all
