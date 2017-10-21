# -*- coding: utf-8 -*-
"""
Created on Sat Sep 09 00:14:39 2017

@author: jcoleman

Load "*.spydata" file and run script.
"""

from matplotlib import pyplot as plt

import matplotlib.image as mpimg

import readroi as roizip    


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
        
        
plotStack(fileDirectory, imgname, roizipname)
multiPlot(response_data,response_avgs,32, 8,range(0,32),1.0)

