# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 20:17:21 2016

@author: jcoleman
"""

import numpy as np
from collections import OrderedDict
import calcium_imaging_data_fast as cidfast

#root = tk.Tk()
#root.update()


def openData(filetype):
    Data = cidfast.FileHandler()
    if filetype == 1 or filetype == 2:
        print (' *** open a rows-columns CSV intensity data file ***')
        int_smooth_red = Data.open_intensity_csv()
        cells = Data.get_cells_from_smoothed(int_smooth_red)        
        return cells
    
    elif filetype == 3:
        print (' *** open a Fiji CSV intensity data file ***')        
        intensity_data = Data.open_intensity_csv()
        cellsraw, cellsmean, areas, xcoords, ycoords = Data.get_cells(intensity_data)
        xycoords = np.array([xcoords, ycoords])  # need to check format  
        return cellsraw, cellsmean, areas, xycoords, Data

def processData(cells, timestamp_file, Data):
    threshold = 2.5        
        
    if timestamp_file == 1:
        print (' *** open MATLAB VisStim Timestamp CSV file ***')
        event_data = Data.open_event_csv()
        channels = Data.get_channels(event_data)
        Data.get_event_stamps(channels, threshold)
    #elif timestamp_file == 2:
        # place new class/functions for extracting "voltage" timestamps
    
    # Run Konnerth lab deltaF/F and EWMA filter - USAGE: (cells = raw intensity data)                      
    t_0 = 0.2
    t_1 = 0.75
    t_2 = 3
    samplingfreq = 30
    
    dff, dff_ewma, dff_offset = cidfast.run_deltaf_ewma(cells, t_0, t_1, t_2, samplingfreq)       
    return channels, dff, dff_ewma, dff_offset

def decodeData(channels, filetype, cells, Data):
    threshold = 2.5
    
    flips = OrderedDict()
    flops = OrderedDict()
    grays = OrderedDict()   
    sessions = OrderedDict()

    timeCodes = cidfast.GetCodeList(channels[2], [channels[2], channels[3]], threshold)

    flopTimeStamps = cidfast.GetTimeStamps(1, timeCodes)
    flipTimeStamps = cidfast.GetTimeStamps(2, timeCodes)
       
    signalChannel = channels[0]
    codeChannels = [channels[2], channels[3]]

    SRP = cidfast.SRPdecoder()
    timeCodes = SRP.GetCodeList(signalChannel, codeChannels)
    flopTimeStamps = SRP.GetTimeStamps(1, timeCodes)
    flipTimeStamps = SRP.GetTimeStamps(2, timeCodes)
    flipLengths, flopLengths = SRP.GetStimLengths(flipTimeStamps, flopTimeStamps)
    avgLength = SRP.AvgStimLength([flipTimeStamps, flopTimeStamps])
    
    SRP.stimLength = avgLength
    
    stimsPerSession = SRP.StimsPerSession(flipLengths, flopLengths, avgLength)

    for i in range(len(cells)):
        flips[i], flops[i], grays[i], sessions[i] = SRP.GetStimLists(cells[i], stimsPerSession, avgLength, flipTimeStamps, flopTimeStamps, Data.Bframe_up)
        
    # get frame indices for flip and flop onsets in first session       
    flip_i = flips[0][0].keys()
    flop_i = flops[0][0].keys()
    # get frame indices for gray onsets in ALL sessions
    gray_i = list()
    for key in grays[0]:
        gray_i.append(grays[0][key]['ts'])
        
    return flips, flops, grays, sessions, flip_i, flop_i, gray_i