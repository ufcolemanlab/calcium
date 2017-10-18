# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 14:04:22 2017

@author: jesse
"""

import numpy as np
import BinFile as BinFile
from matplotlib import pyplot as plt

class StepCodeFile:
    def __init__(self, filename, codes, channels, frequency):
        bf = BinFile.BinFile(filename)
        bf.setTotalChannels(channels)
        bf.open_bin()
        
        self.sampling_frequency = frequency 
        
        self.data = bf.data
        self.timestamps = [list() for i in range(len(bf.data))]
        
        self.timestamp_channel = 0
        self.event_channel = 1
        self.Bframe_channel = 2
        
        self.sampling_window = 50
        self.sampling_window_offset = 200
        
        self.max_voltage = max(self.data[1])
        self.onset_thresh = 0.05 * self.max_voltage
    
        self.level_codes = [i+1 for i in range(codes)]
        self.level_thresholds = [self.max_voltage * (i + 1)/codes for i in range(codes)]
        
        assert len(self.level_codes) == len(self.level_thresholds)
        
        self.get_timestamps()
        
        
        self.levels = list()
                
        self.get_levels()
        
        #self.reduce_frequency()

        self.stimuli = [x[2] for x in self.levels]
        instance_0 = self.stimuli[0]
        for i in range(1, codes + 1):
            assert self.stimuli.count(i) == self.stimuli.count(instance_0)

    #Gives last sample point before upswing
    def get_onset(self, chan):
        return [i for i in range(len(chan)-1) if chan[i+1] > self.onset_thresh and chan[i] < self.onset_thresh]
    
    #Gives last sample point before downswing
    def get_offset(self, chan):
        return [i for i in range(len(chan)-1) if chan[i+1] < self.onset_thresh and chan[i] > self.onset_thresh]
    
    #Creates (Onset, Offset) pairs
    def get_timestamps(self):
        for i in range(len(self.timestamps)):
            self.timestamps[i] = zip(self.get_onset(self.data[i]), self.get_offset(self.data[i]))
    
    #Estimates voltage level by averaging on event stamp channel
    def get_step_value(self, chan, index):
        limit = self.sampling_window / 2
        return np.mean(chan[index - limit : index + limit])
    
    #onset, offset, level_code, level
    def get_levels(self):
        for stamp in self.timestamps[self.timestamp_channel]:
            middle_offset = self.sampling_window_offset + self.sampling_window/2
            level = self.get_step_value(self.data[self.event_channel], stamp[0] + middle_offset)
            
            vals = np.abs(np.array(self.level_thresholds) - level)
            code = self.level_codes[np.argmin(vals)]
            
            self.levels.append((stamp[0], stamp[1], code, level))
    
    def reduce_frequency(self):
        frequency_factor = self.sampling_frequency / 1000.0
        for i in range(len(self.levels)):
            self.levels[i] = tuple((self.levels[i][0] / frequency_factor, self.levels[i][1] / frequency_factor, self.levels[i][2], self.levels[i][3]))
    
    def get_event_timestamps(self):
        pass
    
    #assigns events to given list of angles. returns dictionary of angle (onset, offset) pairs
    def get_stim_angles(self, angles):
        stims = dict()
#        angle_codes = dict()
        assert len(angles) == len(self.level_codes)
        for angle in angles:
            stims[angle] = list()
        for code_tuple in self.levels:
            event_code = code_tuple[2]
            angle = angles[event_code-1]
            stims[angle].append((code_tuple[0],code_tuple[1]))
        return stims
            
            

if __name__ == "__main__":
    code_file = StepCodeFile('C:/Users/jesse/Documents/calcium data/drive-download-20170829T182520Z-001/new data/mThy6s2_alldrift_002_8_data.bin',8,8,1000)
    stims = code_file.get_stim_angles([0,45,90,135,180,225,270,315])
    plt.plot(code_file.data[0])
    plt.plot(code_file.data[1])
    plt.plot(code_file.data[2])