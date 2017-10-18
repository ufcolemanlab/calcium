# -*- coding: utf-8 -*-
"""
calcium_imaging_data_fast.py
Library of functions for calcium imaging data

Created on Wed Oct 12 16:20:15 2016

@author: jcoleman, jtrinity
copyright UF Pediatrics
"""

import numpy as np
import matplotlib.pyplot as plt
import tkFileDialog
import csv
from copy import deepcopy
from collections import OrderedDict

class SRPdecoder:
    def __init__(self, *args, **kwargs):

        self.voltageThreshold = 2.5
        self.stimLength = 500
        self.baseline = 25


    #Gives back chunks of the input signalChannel corresponding to flip and flop stim lists
    def GetStimLists(self, signalChannel, numStims, avgLength, flipTimeStamps, flopTimeStamps, frame_timestamps):
        frame_timestamps = np.array(frame_timestamps).astype(np.float)
        flips = OrderedDict()
        flops = OrderedDict()
        grays = OrderedDict()
        sessions = OrderedDict()
        block = 0
        
        last_flop_frame_end = None
                    
        for i in range(len(flopTimeStamps)):

            if self.StimLength(flopTimeStamps[i]) < 0.5*avgLength:
                first_flip_frame = (np.abs(frame_timestamps - flipTimeStamps[i][0])).argmin()
                current_flip_block = OrderedDict()
                current_flop_block = OrderedDict()
                current_gray_block = OrderedDict()

                if last_flop_frame_end is not None:                 
                    #first_flip_frame = (np.abs(frame_timestamps - flipTimeStamps[i][0])).argmin()
                    #current_gray_block[(last_flop_frame_end, first_flip_frame)] = signalChannel[last_flop_frame_end:first_flip_frame]
                    grays[block - 1] = {"ts":(last_flop_frame_end, first_flip_frame), "data":signalChannel[last_flop_frame_end:first_flip_frame]}
                for j in range(numStims):
                    flip_frame = (np.abs(frame_timestamps - flipTimeStamps[i+j][0])).argmin()
                    flip_frame_end = (np.abs(frame_timestamps - (flipTimeStamps[i+j][0] + self.stimLength))).argmin()
                    current_flip_block[(flip_frame, flip_frame_end)] = (signalChannel[flip_frame:flip_frame_end])
                    flips[block] = current_flip_block.copy()
                for j in range(numStims - 1):
                    flop_frame = (np.abs(frame_timestamps - flopTimeStamps[1+i+j][0])).argmin()
                    flop_frame_end = (np.abs(frame_timestamps - (flopTimeStamps[1+i+j][0] + self.stimLength))).argmin()
                    #gray_frame_end = (np.abs(frame_timestamps - (flopTimeStamps[1+i+j][0] + 2 * self.stimLength))).argmin()
                    current_flop_block[(flop_frame, flop_frame_end)] = (signalChannel[flop_frame:flop_frame_end])
                    #current_gray_block[(flop_frame_end, gray_frame_end)] = signalChannel[flop_frame_end:gray_frame_end]
                #split the last long block into flip/flop
                last_flop_frame = (np.abs(frame_timestamps - (flipTimeStamps[i + numStims - 1][0] + self.stimLength))).argmin()
                last_flop_frame_end = (np.abs(frame_timestamps - (flipTimeStamps[i + numStims - 1][0] + 2 * self.stimLength))).argmin()
                #last_gray_frame_end = (np.abs(frame_timestamps - (flipTimeStamps[i + numStims - 1][0] + 3 * self.stimLength))).argmin()
                current_flop_block[(last_flop_frame, last_flop_frame_end)] = (signalChannel[last_flop_frame:last_flop_frame_end])
                #current_gray_block[(last_flop_frame_end, last_gray_frame_end)] = signalChannel[last_flop_frame_end:last_gray_frame_end]
                
                flops[block] = current_flop_block.copy()
                
                sessions[block] = {"ts":(first_flip_frame, last_flop_frame_end),"data":signalChannel[first_flip_frame:last_flop_frame_end]}
                
                block = block + 1
                current_flip_block.clear()
                current_flop_block.clear()
                current_gray_block.clear()
                
        gray_length = grays[0]["ts"][1] - grays[0]["ts"][0]
        current_gray_block[(last_flop_frame_end, last_flop_frame_end + gray_length)] = signalChannel[last_flop_frame_end:last_flop_frame_end + gray_length]
        grays[block - 1] = {"ts":(last_flop_frame_end, last_flop_frame_end + gray_length), "data":signalChannel[last_flop_frame_end:last_flop_frame_end + gray_length]}
                
                
        return flips, flops, grays, sessions
    
    #Averages flips/flops in sequential groups of numStims
    def GetAverages(self, flips, numStims):
        avgs = []
        for i in range(0, len(flips), numStims):
            avg = np.zeros(len(flips[0]))
            for j in range(numStims):
                avg += np.array(flips[i+j]-np.average(flips[i+j][:self.baseline]))
            avg /= numStims
            #avg += i*np.ones(len(flips[0]))
            avgs.append(avg)
            
        return avgs
            
                
    #Detects the number of stims in a session
    def StimsPerSession(self, flipLengths, flopLengths, avgLength):
        start = 0
        end = 0
        while(flopLengths[start] > avgLength*0.5):
            start += 1
        while(flipLengths[end] < avgLength*1.5):
            end += 1
        return end - start + 1
        
    #Make this take and return list of lists for using more than 2 channels (maybe a dictionary)
    def GetStimLengths(self, flipTimeStamps, flopTimeStamps):
        flipLengths = [self.StimLength(flip) for flip in flipTimeStamps]
        flopLengths = [self.StimLength(flop) for flop in flopTimeStamps]
        return flipLengths, flopLengths
        
    #takes a list of encoding channels and gives the code at sample point i
    def Decode(self, i, channels):
        return sum(2**a for (a,j) in enumerate([channel[i] for channel in channels]) if j > self.voltageThreshold)
    
    #Returns a code list corresponding to the sample in the signal
    #passes signal for standardization of channel length
    def GetCodeList(self, sig, channels):
        return[self.Decode(i, channels) for i in range(len(sig))]
            
    #Gives total number of sessions
    def GetTotalSessions(self, flopTimeStamps):
        return sum(self.StimLength(i) < 0.5*self.AvgStimLength([flopTimeStamps]) for i in flopTimeStamps)
    
    #Gives (rise, fall) timestamps corresponding to given code (1 for flop, 2 for flip)
    def GetTimeStamps(self, code, timeCodes):    
        rise = [i for i in range (1, len(timeCodes)) if timeCodes[i] == code and timeCodes[i-1] != code]
        fall = [i for i in range (1, len(timeCodes)) if timeCodes[i] != code and timeCodes[i-1] == code]
        return zip(rise, fall)
        
    #gives average stim length. input must be in a list, i.e. f([flipTimeStamps]) or f([flipTimes, flopTimes]) but not f(flipTimeStamps)
    def AvgStimLength(self, timeStampsLists):
        return sum(sum(y-x for (x,y) in stampsList)/len(stampsList) for stampsList in timeStampsLists)/len(timeStampsLists)
        
    def StimLength(self, stim):
        return stim[1]- stim[0]
        

        
class FileHandler():
    def __init__(self):
        self.strobe_up = list()
        self.strobe_down = list()
        self.onset_up = list()
        self.onset_down = list()
        
        self.stim_up = list()
        self.stim_down = list()
        
        self.Bframe_up = list()
        self.Bframe_down = list()
        
    
    #Opens event channel Data from csv
    def open_event_csv(self):
        files = tkFileDialog.askopenfilenames()
        if files:
            files = list(files)
            with open(files[0], 'rb') as csvfile:
                reader = csv.reader(csvfile, delimiter = ',')
                reader = list(reader)
                r = [row[2:9] for row in reader]
                r = r[8:]
                r = np.array(r).astype(np.float64)
                
        return r
    
    #Opens intensity Data from csv
    def open_intensity_csv(self):
        files = tkFileDialog.askopenfilenames()
        if files:
            files = list(files)
            with open(files[0], 'rU') as csvfile:
                reader = csv.reader(csvfile, delimiter = ',')
                reader = list(reader)
                reader = np.array(reader)
                
        return reader
        
    #Gets indices of upward value swing
    def get_upswing(self, channel, threshold):
        return [index for index in range(1, len(channel)) if (channel[index] > threshold) and (channel[index - 1] < threshold)]
    
    #Gets indices of downward value swing
    def get_downswing(self, channel, threshold):
        return [index for index in range(1, len(channel)) if (channel[index] < threshold) and (channel[index - 1] > threshold)]
    
    #Gets cell data, roi area, and x-y coordinate arrays from intensity csv
    def get_cells(self, c):
        cellsraw = list()
        cellsmean = list()
        areas = list()
        xcoords = list()
        ycoords = list()
        for i in range(len(c[0])):
            if 'RawIntDen' in c[0][i]:
                cellsraw.append(np.array(c[1:, i]).astype(np.float64))
            if 'Mean' in c[0][i]:
                cellsmean.append(np.array(c[1:, i]).astype(np.float64))
            if 'Area' in c[0][i]:
                areas.append(np.array(c[1:, i]).astype(np.float64))
            if 'X' in c[0][i]:
                xcoords.append(np.array(c[1:, i]).astype(np.float64))
            if 'Y' in c[0][i]:
                ycoords.append(np.array(c[1:, i]).astype(np.float64))
        return cellsraw, cellsmean, areas, xcoords, ycoords
    
    def get_cells_from_smoothed(self, c):
        cells = list()
        for i in range(len(c[0])):
            cells.append(np.array(c[1:, i]).astype(np.float64))
        return cells
    
    #Gets channels from event csv 
    def get_channels(self, r):
        channels = [r[0:, i] for i in range(len(r[0]))]
        return channels
    
    #Gets event timestamp data from channels
    def get_event_stamps(self, channels, threshold):
        self.strobe_up = self.get_upswing(channels[1], threshold)
        self.strobe_down = self.get_downswing(channels[1], threshold)
        self.onset_up = self.get_upswing(channels[2], threshold)
        self.onset_down = self.get_downswing(channels[2], threshold)
        
        self.stim_up = self.get_upswing(channels[3], threshold)
        self.stim_down = self.get_downswing(channels[3], threshold)
        
        self.Bframe_up = self.get_upswing(channels[6], threshold)
        self.Bframe_down = self.get_downswing(channels[6], threshold)
    
    #Gets nearest value in array to input
    def get_nearest(self, x, array):
        keylist = np.array(array)
        index = (np.abs(keylist-x)).argmin()
        return keylist[index]
    
    #Gives maximum difference in ms bewteen frame onsets and stim onset
    def max_frame_error(self, frame_timestamps, event_timestamps):
        deviations = list()
        for ts in event_timestamps:
            dev = np.abs(np.array(frame_timestamps) - ts).min()
            deviations.append(dev)
        return deviations
    
    #Gives flip and flop dictionaries
    #key is (start, end) frames for that stimulus
    #value is intensity data
    def get_flip_flops(self, cell, strobe_timestamps, stim_timestamps, frame_timestamps):
        #cell, Data.strobe_down, Data.stim_up, Data.Bframe_up
        duration = strobe_timestamps[1] - strobe_timestamps[0]
        flips = OrderedDict()
        flops = OrderedDict()
        gray = OrderedDict()
        for ts in stim_timestamps:
            closest_flip = self.get_nearest(ts, strobe_timestamps)
            closest_flop = self.get_nearest(ts + duration, strobe_timestamps)
            
            flip_frame_ts = self.get_nearest(closest_flip, frame_timestamps)
            flop_frame_ts = self.get_nearest(closest_flop, frame_timestamps)
            
            flip_frame_ts_end = self.get_nearest(closest_flip + duration, frame_timestamps)
            flop_frame_ts_end = self.get_nearest(closest_flop + duration, frame_timestamps)
            gray_frame_ts_end = self.get_nearest(closest_flop + 2*duration, frame_timestamps)
            
            #frames are counted by how many times Bframe channel pulses up and down
            #target frame number should be the array index for the nearest timestamp
            flip_frame = (np.abs(frame_timestamps - flip_frame_ts)).argmin()
            flop_frame = (np.abs(frame_timestamps - flop_frame_ts)).argmin()
            gray_frame = (np.abs(frame_timestamps - gray_frame_ts_end)).argmin()
            
            flip_frame_end = (np.abs(frame_timestamps - flip_frame_ts_end)).argmin()
            flop_frame_end = (np.abs(frame_timestamps - flop_frame_ts_end)).argmin()
            
            
            flips[(flip_frame, flip_frame_end)] = cell[flip_frame:flip_frame_end]
            flops[(flop_frame, flop_frame_end)] = cell[flop_frame: flop_frame_end]
            gray[(flop_frame_end, gray_frame)] = cell[flop_frame_end: gray_frame]
            
        return flips, flops, gray
    
    def get_stim_block(self, cell, stim_up, stim_down, frame_timestamps):
        stim_up.sort()
        stim_down.sort()
        stamps = zip(stim_up, stim_down)
        blocks = OrderedDict()
        for ts in stamps:
            onset_frame_ts = self.get_nearest(ts[0], frame_timestamps)
            offset_frame_ts = self.get_nearest(ts[1], frame_timestamps)
            
            frame_start = (np.abs(frame_timestamps - onset_frame_ts)).argmin()
            frame_end = (np.abs(frame_timestamps - offset_frame_ts)).argmin()
            
            blocks[(frame_start, frame_end)] = cell[frame_start:frame_end]
        return blocks
            
        
    #gets the average of flips or flops for a cell
    def get_avg(self, flips):
        return np.mean(np.array(flips).astype(np.float64), axis = 0)

def Decode(i, channels, thresh):
    return sum(2**a for (a,j) in enumerate([channel[i] for channel in channels]) if j > thresh)
        
#Returns a code list corresponding to the sample in the signal
def GetCodeList(sig, channels, thresh):
    return[Decode(i, channels, thresh) for i in range(len(sig))]

#Gives (rise, fall) timestamps corresponding to given code (1 for flop, 2 for flip)
def GetTimeStamps(code, timeCodes):    
    rise = [i for i in range (1, len(timeCodes)) if timeCodes[i] == code and timeCodes[i-1] != code]
    fall = [i for i in range (1, len(timeCodes)) if timeCodes[i] != code and timeCodes[i-1] == code]
    return zip(rise, fall)        

def run_deltaf_ewma(data, t_0, t_1, t_2, samplingfreq):         
    import process_function_jc as pf    
    """
    From Konnerth lab Nature Protocols paper, for 30Hz:
    t_0 = 0.2
    t_1 = 0.75
    t_2 = 3
    samplingfreq = 30
    """   
    dff = OrderedDict()
    dff_ewma = OrderedDict()
    dff_offset = []
    for i in range(len(data)):
        dff[i], dff_ewma[i], dff_offset = pf.process_function(data[i], t_0, t_1, t_2, samplingfreq)
    if dff_offset > 0:
        dffnans = np.zeros(dff_offset)
        dffnans[:] = np.NAN
        for j in range(len(dff)):
            dff[j] = np.append(dffnans,dff[j])
            dff_ewma[j] = np.append(dffnans,dff_ewma[j])           
    return dff, dff_ewma, dff_offset
    
def chunkSessions(sessions, grays):
    # Chunk up stim and gray sessions
    stimsession_cell = list()
    for cell in sessions:
        temp_session = np.array([sessions[cell][s]["data"] for s in sessions[cell]])
        stimsession_cell.append(temp_session)        
    graysession_cell = list()
    for cell in grays:
        temp_gray = np.array([grays[cell][s]["data"] for s in grays[cell]])
        graysession_cell.append(temp_gray)        
    return stimsession_cell, graysession_cell
    
def combineGrayStim(stimsession_cell, graysession_cell, flop_i, gray_i):
    # Combine STIM and GRAY session DATA (all sessions)
    sessions_stimgray = list()
    for i in range(len(stimsession_cell)):
        temp_stimgray = np.array([np.concatenate((stimsession_cell[i][s], graysession_cell[i][s]), axis=0) for s in range(len(stimsession_cell[i]))]) # zip?
        sessions_stimgray.append(temp_stimgray)   
    # Combine GRAY and STIM session DATA (sessions 2-6 only)
    sessions_graystim = list()
    for i in range(len(stimsession_cell)):
        temp_graystim = np.array([np.concatenate((graysession_cell[i][s-1], stimsession_cell[i][s]), axis=0) for s in range(1,len(stimsession_cell[i]))]) # zip?
        sessions_graystim.append(temp_graystim)        
    # Indicii for gray onset (for avg sesssion, after last flop, which is the max from flop_i + the time interval for the last flop +1)
    gray_session_onset_i = 1 + (flop_i[-1][1]+(flop_i[-1][1] - flop_i[-2][1]))  
    # Need indicii for STIM onset (for avg sesssion, after gray, which is the max from gray_i + the time interval for the last flop +1)
    stim_session_onset_i = (gray_i[0][1]-gray_i[0][0])-1
    return sessions_stimgray, sessions_graystim, gray_session_onset_i, stim_session_onset_i
    
def meanGrayStimSort(sessions_graystim):
    # Find the min length of graystim sessions for averaging, calc. means and sort
    minlength_graystimsession = list()
    for cell in range(len(sessions_graystim)):
        graystim_sessionlengths = list()
        for i in range(len(sessions_graystim[cell])):
            tempval = len(sessions_graystim[cell][i])
            graystim_sessionlengths.append(tempval)
        minlength_graystimsession.append(min(graystim_sessionlengths))
    minlength_graystimsession=min(minlength_graystimsession)  
    avgs_graystim = list()
    for cell in range(len(sessions_graystim)):
        tempval=list()
        for s in range(len(sessions_graystim[cell])):
            tempval.append(sessions_graystim[cell][s][0:minlength_graystimsession])
        avgs_graystim.append(np.nanmean(tempval, axis = 0))
    # zero avgs
    for i in range(len(avgs_graystim)):
        tempdenom = np.nanmean(avgs_graystim[i])
        avgs_graystim[i] = avgs_graystim[i]-tempdenom       
    # Sort cells by latency to max intensity value        
    sorted_graystimavgs = sorted(avgs_graystim, key=lambda x: x.argmax())
    sorted_graystimavgs_keys = np.argsort(np.argmax(avgs_graystim, axis=1))
    sorted_graystimavgs_axis_labels = [x+1 for x in sorted_graystimavgs_keys]       
    #Normalize all data (0-1)    
    norm_graystimavgs = deepcopy(sorted_graystimavgs) #.copy()
    for key in range(len(norm_graystimavgs)):
        norm_graystimavgs[key] -= min(norm_graystimavgs[key])
        norm_graystimavgs[key] /= max(norm_graystimavgs[key])
        norm_graystimavgs[key] = np.array(norm_graystimavgs[key]).astype(np.float)       
    return minlength_graystimsession, avgs_graystim, sorted_graystimavgs, norm_graystimavgs, sorted_graystimavgs_keys, sorted_graystimavgs_axis_labels
    
    
def plotSessionTraces(nrows, ncols, sessiondata, avgdata, stimonset, plottitle, mode):
    if mode == 'share xy':
        fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all')    
    elif mode == 'share off':
        fig, axes = plt.subplots(nrows, ncols)           
    for i in range(nrows):
        for j in range(ncols):
            try:
                for k in range(len(sessiondata[i*ncols+j])):
                    axes[i,j].plot(sessiondata[i*ncols+j][k])
                    axes[i,j].plot(avgdata[i*ncols+j], linewidth=2, color='k')
                    axes[i,j].text(.5, .8,'Cell '+str(i*ncols+j+1), ha='left', va='center', transform=axes[i,j].transAxes)
                    axes[i,j].axvline(x=stimonset, linewidth=0.5, color='k')
                if i==0 and j==0 and mode == 'share xy':
                    axes[i,j].tick_params(axis=u'both', which=u'both',length=0)
                elif i>0 or j>0 and mode == 'share xy':
                    axes[i,j].axis('off')
            except IndexError:
                pass
    fig.suptitle(plottitle, fontsize=11)
    
          
def subplotCellTraces(nrows, ncols, avgdata, skipval, stimonset, mode):  #skipval = number of traces per plot
    if mode == 'share xy':
        fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all') 
    elif mode == 'share off':
        fig, axes = plt.subplots(nrows, ncols)
    index = 0
    for i in range(nrows):
        for j in range(ncols):            
            for k in range(skipval):
                try:
                    axes[i,j].plot(avgdata[index])
                    axes[i,j].axvline(x=stimonset, linewidth=0.5, color='k')
                    if i==0 and j==0 and mode == 'share xy':
                        axes[i,j].tick_params(axis=u'both', which=u'both',length=0)
                    elif i>0 or j>0 and mode == 'share xy':
                        axes[i,j].axis('off')
                    if mode == 'share off':
                        axes[i,j].tick_params(axis=u'both', which=u'both',length=0)
                except IndexError:
                    pass
                #print(str(index))
                index += 1
                
                
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
        
