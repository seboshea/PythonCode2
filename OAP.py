# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:54:22 2019

@author: Admin TEMP
"""

import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import math
import bisect
import h5py
import scipy.io as sio
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
#from matplotlib.mlab import bivariate_normal
from matplotlib.dates import DateFormatter
import os
import gc
from GrayscaleFunctions import GetIAT_TimeInS_vector, BinHaloHolo_v2
from MyFunctions import datetime_2_s_sincemidnight
from PICASSO_functions import Avg_CDP_dndDp


#// Constants to define the single particle time stamps
#StrConstant KS_OASIS_PTIME_MATRIX = "ParticleTimesWv"
#Constant K_OASIS_PTIME_COLS = 8
#Constant K_OASIS_PTIME_COL_SECS = 0
#Constant K_OASIS_PTIME_COL_NS = 1
#Constant K_OASIS_PTIME_COL_SLCS = 2
#Constant K_OASIS_PTIME_COL_PNUM = 3
#Constant K_OASIS_PTIME_COL_DOF = 4
#Constant K_OASIS_PTIME_COL_TAS = 5
#Constant K_OASIS_PTIME_COL_BUFFN = 6
#Constant K_OASIS_PTIME_COL_IMOUT = 7
#
#// Constant to define the single particle stats output
#StrConstant KS_OASIS_PSTAT_MATRIX = "ParticleStatsWv"
#Constant K_OASIS_PSTAT_COLS = 20
#Constant K_OASIS_PSTAT_COL_NUM = 0
#Constant K_OASIS_PSTAT_COL_AREA = 1
#Constant K_OASIS_PSTAT_COL_AREA2 = 2
#Constant K_OASIS_PSTAT_COL_XMIN = 3
#Constant K_OASIS_PSTAT_COL_XMAX = 4
#Constant K_OASIS_PSTAT_COL_YMIN = 5
#Constant K_OASIS_PSTAT_COL_YMAX = 6
#Constant K_OASIS_PSTAT_COL_PERIM = 7
#Constant K_OASIS_PSTAT_COL_SLC = 8
#Constant K_OASIS_PSTAT_COL_XMIN_BG = 9
#Constant K_OASIS_PSTAT_COL_XMAX_BG = 10
#Constant K_OASIS_PSTAT_COL_YMIN_BG = 11
#Constant K_OASIS_PSTAT_COL_YMAX_BG = 12
#Constant K_OASIS_PSTAT_COL_N_EDGE0 = 13
#Constant K_OASIS_PSTAT_COL_N_EDGE1 = 14
#Constant K_OASIS_PSTAT_COL_INTVOID = 15
#Constant K_OASIS_PSTAT_COL_N_EDGE0_BG = 16
#Constant K_OASIS_PSTAT_COL_N_EDGE1_BG = 17
#Constant K_OASIS_PSTAT_COL_CENTRE = 18
#Constant K_OASIS_PSTAT_COL_MAXLEN_BG = 19
#
#// Constant to define the single particle stats output
#StrConstant KS_OASIS_PGRAY_MATRIX = "ParticleStatsWv"
#Constant K_OASIS_PGRAY_COLS = 15
#Constant K_OASIS_PGRAY_COL_LVL0 = 20
#Constant K_OASIS_PGRAY_COL_LVL1 = 21
#Constant K_OASIS_PGRAY_COL_LVL2 = 22
#Constant K_OASIS_PGRAY_XMN_LVL0 = 23
#Constant K_OASIS_PGRAY_XMX_LVL0 = 24
#Constant K_OASIS_PGRAY_YMN_LVL0 = 25
#Constant K_OASIS_PGRAY_YMX_LVL0 = 26
#Constant K_OASIS_PGRAY_XMN_LVL1 = 27
#Constant K_OASIS_PGRAY_XMX_LVL1 = 28
#Constant K_OASIS_PGRAY_YMN_LVL1 = 29
#Constant K_OASIS_PGRAY_YMX_LVL1 = 30
#Constant K_OASIS_PGRAY_XMN_LVL2 = 31
#Constant K_OASIS_PGRAY_XMX_LVL2 = 32
#Constant K_OASIS_PGRAY_YMN_LVL2 = 33
#Constant K_OASIS_PGRAY_YMX_LVL2 = 34

#__________________________________________________________________________________
# Loop through OASIS output files in a given folder

def BatchBothChannels():

#    Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/'
#   Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Test/'
#    ExpDate= datetime.datetime(2017, 2, 7, 0, 0, 0)
#    Path2DS = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/oasisoutput/Colocation/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/Clarify/C052/2DS/Output/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/Clarify/C031/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/DropletGun/Dropgun_2ds/Batch1/90um/Run1_20180322_1038_vchannel/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/DropletGun/Dropgun_2ds/Batch1/120um/Run4_20180321_1523_hchannel/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/DropletGun/Dropgun_2ds/Batch1/120um/Run5_20180321_1536_hchannel/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/DropletGun/Dropgun_2ds/Batch1/120um/Run3_20180321_1503_vchannel/'
    Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C097/2DS/Output/'
#    Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C073/2DS/OasisOut/'
    for filena in os.listdir(Path2DS):
        if filena.endswith(".h5") and filena.startswith('base'):
            print(filena)
            FindParticlesOnBothChannels(Path2DS,filena)
            HistMIDPosition2DS(Path2DS,filena)
            #ShatteringDiagnostics(Path2DS,filena)

    

#__________________________________________________________________________________
#
            
#Histogram on mean particle position across the array
def HistMIDPosition2DS(Path2DS,filena):
    PixelSize =10
    MaxAspectRatio = 10
    Streak = 5 
    
    Data_h5 = h5py.File(Path2DS+filena, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    #TmpStats=np.array(Data_h5['PartGrayStatsWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()

    DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    StreakFlag =(np.where(np.logical_and(DiameterBGx == PixelSize, DiameterBGy >= Streak*PixelSize), 1, 0))   
    AspectRatio = DiameterBGy / DiameterBGx

    Channel = ParticleTimesWv[:,4] 
    Streak_Ch0 = StreakFlag[Channel == 0]
    Streak_Ch1 = StreakFlag[Channel == 1]
    AspectRatio_Ch0 = AspectRatio[Channel == 0]
    AspectRatio_Ch1 = AspectRatio[Channel == 1] 

    # Nan edge particles
    shape = np.shape(ParticleStatsWv)
    Edge= np.zeros(shape[0])
    for i in range(shape[0]):
        if (ParticleStatsWv[i,3]==0 or ParticleStatsWv[i,4]==127):
            Edge[i] = 1
    Edge_CH0 = Edge[Channel == 0]
    Edge_CH1 = Edge[Channel == 1]
    
    MIDx = (ParticleStatsWv[:,10] + ParticleStatsWv[:,9])/2
    MIDx_Ch0 = MIDx[Channel == 0]
    MIDx_Ch1 = MIDx[Channel == 1]

    # Nan edge particles
    MIDx_Ch0[Edge_CH0 == 1] = np.nan
    MIDx_Ch1[Edge_CH1 == 1] = np.nan
    
    # Nan streak particles
    MIDx_Ch0[Streak_Ch0 == 1] = np.nan
    MIDx_Ch1[Streak_Ch1 == 1] = np.nan
    MIDx_Ch0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    MIDx_Ch1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan
    
    HistBins = np.linspace(0,128,num=129, endpoint=True)
    HistBinsMid = (HistBins[:-1] + HistBins[1:]) / 2
    MIDx_Ch0_hist, tmp= np.histogram(MIDx_Ch0,bins=HistBins)
    MIDx_Ch1_hist, tmp= np.histogram(MIDx_Ch1,bins=HistBins)
    
    percentile = 1.1*max(np.percentile(MIDx_Ch0_hist,98),np.percentile(MIDx_Ch1_hist,99))
    
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 12})
        plt.plot(HistBinsMid,MIDx_Ch0_hist, label= 'Ch = 0')
        plt.plot(HistBinsMid,MIDx_Ch1_hist, label= 'Ch = 1' )
        plt.xlabel('Array elements')
        plt.ylabel('Counts')
        plt.ylim([0,percentile]) #restrict y to below 90 percentile
        plt.legend()
        plt.savefig(Path2DS+filena[:-3]+'_MIDx_edgeRej.png',dpi=200)
        plt.close(fig)


#__________________________________________________________________________________
#
# Count number of shattered particles per shattering event.

 
def ShatteringDiagnostics(Path2DS,filena):
    ColocationThreshold = 5E-7
    IAT_treshold = 1E-5
    SaveFile =1
    PixelSize = 10 
    MaxAspectRatio = 10
    Streak = 5
    MeanXYFlag = 1
    ZdMax = 2
    Lambda_um = 0.658 

    Data_h5 = h5py.File(Path2DS+filena, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()
    
    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    Channel = ParticleTimesWv[:,4] 
    MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)
    DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    MaxDiameterBG= np.sqrt(DiameterBGx**2 + DiameterBGy**2)
    StreakFlag =(np.where(np.logical_and(DiameterBGx == PixelSize, DiameterBGy >= Streak*PixelSize), 1, 0))   
    AspectRatio = DiameterBGy / DiameterBGx
    
    #Select channels
    MeanXY_BG_CH0 = MeanXY_BG[Channel == 0]
    Seconds_CH0 = ParticleTimeSeconds[Channel == 0]
    #Edge_CH0 = Edge[Channel == 0]
    MeanXY_BG_CH1 = MeanXY_BG[Channel == 1]
    Seconds_CH1 = ParticleTimeSeconds[Channel == 1]
    #Edge_CH1 = Edge[Channel == 1]
    Streak_Ch0 = StreakFlag[Channel == 0]
    Streak_Ch1 = StreakFlag[Channel == 1]
    AspectRatio_Ch0 = AspectRatio[Channel == 0]
    AspectRatio_Ch1 = AspectRatio[Channel == 1] 
    MaxDiameterBG_Ch0 = MaxDiameterBG[Channel == 0]
    MaxDiameterBG_Ch1 = MaxDiameterBG[Channel == 1]
#    MIDx_Ch0 = MIDx[Channel == 0]
#    MIDx_Ch1 = MIDx[Channel == 1]
    
    Streak_Ch0 = StreakFlag[Channel == 0]
    Streak_Ch1 = StreakFlag[Channel == 1]
    AspectRatio_Ch0 = AspectRatio[Channel == 0]
    AspectRatio_Ch1 = AspectRatio[Channel == 1] 
    MaxDiameterBG_Ch0 = MaxDiameterBG[Channel == 0]
    MaxDiameterBG_Ch1 = MaxDiameterBG[Channel == 1]

    for i in range(len(Seconds_CH0)): 
        if (AspectRatio_Ch0[i] >= MaxAspectRatio or Streak_Ch0[i] == 1) : 
            Seconds_CH0[i] = np.nan
            
    for i in range(len(Seconds_CH1)): 
        if (AspectRatio_Ch1[i] >= MaxAspectRatio or Streak_Ch1[i] == 1) : 
            Seconds_CH1[i] = np.nan  
            
    MeanXY_BG_CH0 = MeanXY_BG_CH0[~np.isnan(Seconds_CH0)]
    MaxDiameterBG_Ch0 = MaxDiameterBG_Ch0[~np.isnan(Seconds_CH0)]
    #Edge_CH0 = Edge_CH0[~np.isnan(Seconds_CH0)]
    Seconds_CH0= Seconds_CH0[~np.isnan(Seconds_CH0)]

    MeanXY_BG_CH1 = MeanXY_BG_CH1[~np.isnan(Seconds_CH1)]
    MaxDiameterBG_Ch1 = MaxDiameterBG_Ch1[~np.isnan(Seconds_CH1)]
    #Edge_CH1 = Edge_CH1[~np.isnan(Seconds_CH1)]       
    Seconds_CH1= Seconds_CH1[~np.isnan(Seconds_CH1)]

    #Get IAT times for each channel seprately
    IAT_CH0 = GetIAT_TimeInS_vector(Seconds_CH0)    
    IAT_CH1 = GetIAT_TimeInS_vector(Seconds_CH1)

    # Count number of shattered particles per event
    Counts = 0
    for i in range(len(IAT_CH0)) :
        if IAT_CH0[i] < IAT_treshold : 
            Counts += 1

    IAT_CH0_logic = IAT_CH0 < IAT_treshold
    ShatteredNPerEventCH0 = np.diff(np.nonzero(np.concatenate(([IAT_CH0_logic[0]],IAT_CH0_logic[:-1] != IAT_CH0_logic[1:])))[0])[::2]

    IAT_CH1_logic = IAT_CH1 < IAT_treshold
    ShatteredNPerEventCH1 = np.diff(np.nonzero(np.concatenate(([IAT_CH1_logic[0]],IAT_CH1_logic[:-1] != IAT_CH1_logic[1:])))[0])[::2]

    bins = np.arange(0.5,20.5,1)
    histCH0,tmp = np.histogram(ShatteredNPerEventCH0, bins)
    histCH1,tmp = np.histogram(ShatteredNPerEventCH1, bins)
    MidBins = (bins[1:]+bins[:-1]) / 2
    
    
    fig=plt.figure(figsize=(7,7)) 
    plt.rcParams.update({'font.size': 12})
    plt.plot(MidBins, histCH0, 'o',label = 'Ch0')
    plt.plot(MidBins, histCH1, 'o',label = 'Ch1')
    plt.yscale('log')
    plt.xlabel('Number of particles per shattering event')
    plt.ylabel('Counts')
    plt.xlim([0,20])
    plt.title('IAT = ' +str(IAT_treshold))
    plt.legend()
    plt.xticks(np.arange(0, 22, step=2))
    plt.savefig(Path2DS+filena[:-3]+'_NshatteredParticles.png',dpi=200)
    plt.close(fig)
    
        # save matching particles
    if SaveFile == 1 : 
        SavePath= Path2DS+'IAT_'+filena
        file = h5py.File(SavePath, 'w')
        file.create_dataset('Seconds_CH0', data=Seconds_CH0)
        file.create_dataset('Seconds_CH1', data=Seconds_CH1)
        file.create_dataset('IAT_CH0', data=IAT_CH0)
        file.create_dataset('IAT_CH1', data=IAT_CH1)
        file.create_dataset('MaxDiameterBG_Ch0', data=MaxDiameterBG_Ch0)
        file.create_dataset('MaxDiameterBG_Ch1', data=MaxDiameterBG_Ch1)
        file.close()

#__________________________________________________________________________________
# Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/'
# filena = 'base180207155526.h5'
        
# Select particles with time seperation below ColocationThreshold on the other channel

def FindParticlesOnBothChannels(Path2DS,filena):
    SaveFile =1
    #Set particle filtering thresholds
    ColocationThreshold = 5E-7
    IAT_treshold = 1E-5
    #IAT_treshold = 1E-6
    PixelSize = 10 
    MaxAspectRatio = 10
    Streak = 5
    MeanXYFlag = 0
    ZdMax = 2
    Lambda_um = 0.785 
    #load OASIS stats
    Data_h5 = h5py.File(Path2DS+filena, 'r')              
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'])
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()
    ImageID = ParticleTimesWv[:,3]
    BufferNumber = ParticleTimesWv[:,6]
    #Get particle buffer time in seconds from midnight. 3 decimal places
    ParticleBufferTimeS= np.zeros(len(BufferNumber))*np.nan
    for i in range(len(BufferNumber)):
        if BufferNumber[i] == 4294967295 :
            ParticleBufferTimeS[i] = np.nan
        else:   
            ParticleBufferTimeS[i] = 3600*HeaderMatrixWv[BufferNumber[i],3]+ 60*HeaderMatrixWv[BufferNumber[i],4] + HeaderMatrixWv[BufferNumber[i],5] + (HeaderMatrixWv[BufferNumber[i],6]/1000)
    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    Channel = ParticleTimesWv[:,4] 
    MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)
    DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    MaxDiameterBG= np.sqrt(DiameterBGx**2 + DiameterBGy**2)
    StreakFlag =(np.where(np.logical_and(DiameterBGx == PixelSize, DiameterBGy >= Streak*PixelSize), 1, 0))   
    AspectRatio = DiameterBGy / DiameterBGx
    #MIDx = (ParticleStatsWv[:,10] + ParticleStatsWv[:,9])/2

    # find the Z at ZdMax 
    if MeanXYFlag == 1: 
        Zum_at_ZdMax = (ZdMax * MeanXY_BG**2) / (4 * Lambda_um)
    else : 
        Zum_at_ZdMax = (ZdMax * MaxDiameterBG**2) / (4 * Lambda_um)

    # Nan edge particles
    # Edge = 0 not touching array edge
    # Edge = 1 touching array edge
    # Edge = 2 Zd < 2 and not touching array edge
    shape = np.shape(ParticleStatsWv)
    Edge= np.zeros(shape[0])
    for i in range(shape[0]):       
        MinElement = np.maximum(0, 63-np.round(Zum_at_ZdMax[i]/PixelSize)) # Create new edge where Zd < 2 
        MaxElement = np.minimum(127,63+np.round(Zum_at_ZdMax[i]/PixelSize))
        if (ParticleStatsWv[i,3]==0 or ParticleStatsWv[i,4]==127):
            Edge[i] = 1
        if (ParticleStatsWv[i,3]>MinElement and ParticleStatsWv[i,4]<MaxElement):
            Edge[i] = 2 # if Zd < 2 

    #Select channels
    MeanXY_BG_CH0 = MeanXY_BG[Channel == 0]
    Seconds_CH0 = ParticleTimeSeconds[Channel == 0]
    Edge_CH0 = Edge[Channel == 0]
    MeanXY_BG_CH1 = MeanXY_BG[Channel == 1]
    Seconds_CH1 = ParticleTimeSeconds[Channel == 1]
    Edge_CH1 = Edge[Channel == 1]
    Streak_Ch0 = StreakFlag[Channel == 0]
    Streak_Ch1 = StreakFlag[Channel == 1]
    AspectRatio_Ch0 = AspectRatio[Channel == 0]
    AspectRatio_Ch1 = AspectRatio[Channel == 1] 
    MaxDiameterBG_Ch0 = MaxDiameterBG[Channel == 0]
    MaxDiameterBG_Ch1 = MaxDiameterBG[Channel == 1]
    ImageID_Ch0 = ImageID[Channel == 0]
    ImageID_Ch1 = ImageID[Channel == 1]
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS[Channel == 0]
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS[Channel == 1]

    print('File length '+ str(len(Channel[~np.isnan(Channel)])))
    print('Channel 0 = '+ str(len(Seconds_CH0[~np.isnan(Seconds_CH0)])))
    print('Channel 1 = '+ str(len(Seconds_CH1[~np.isnan(Seconds_CH1)])))

    #remove streaking
    for i in range(len(Seconds_CH0)): 
        if (AspectRatio_Ch0[i] >= MaxAspectRatio or Streak_Ch0[i] == 1) : 
            Seconds_CH0[i] = np.nan
            
    for i in range(len(Seconds_CH1)): 
        if (AspectRatio_Ch1[i] >= MaxAspectRatio or Streak_Ch1[i] == 1) : 
            Seconds_CH1[i] = np.nan     
    MeanXY_BG_CH0 = MeanXY_BG_CH0[~np.isnan(Seconds_CH0)]
    MaxDiameterBG_Ch0 = MaxDiameterBG_Ch0[~np.isnan(Seconds_CH0)]
    Edge_CH0 = Edge_CH0[~np.isnan(Seconds_CH0)]
    MeanXY_BG_CH1 = MeanXY_BG_CH1[~np.isnan(Seconds_CH1)]
    MaxDiameterBG_Ch1 = MaxDiameterBG_Ch1[~np.isnan(Seconds_CH1)]
    Edge_CH1 = Edge_CH1[~np.isnan(Seconds_CH1)]       
    ImageID_Ch0 = ImageID_Ch0[~np.isnan(Seconds_CH0)] 
    ImageID_Ch1 = ImageID_Ch1[~np.isnan(Seconds_CH1)] 
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[~np.isnan(Seconds_CH0)]
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[~np.isnan(Seconds_CH1)] 
    
    Seconds_CH0= Seconds_CH0[~np.isnan(Seconds_CH0)]
    Seconds_CH1= Seconds_CH1[~np.isnan(Seconds_CH1)]

    #Remove shattering, Get IAT times for each channel seprately
    IAT_CH0 = GetIAT_TimeInS_vector(Seconds_CH0)    
    IAT_CH1 = GetIAT_TimeInS_vector(Seconds_CH1)
    
    # Plot IAT for both channels indpendently
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 12})
        HistBins = np.logspace(-8,0,num=100)
        IAT_CH0_hist, tmp= np.histogram(IAT_CH0,bins=HistBins)
        HistBinsMid = (HistBins[:-1] + HistBins[1:]) / 2
        plt.plot(HistBinsMid,IAT_CH0_hist,label='Ch=0')       
        IAT_CH1_hist, tmp= np.histogram(IAT_CH1,bins=HistBins)
        plt.plot(HistBinsMid,IAT_CH1_hist,label='Ch=1')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('IAT, s')
        plt.ylabel('Counts')
        plt.legend()
        plt.title('IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
        plt.savefig(Path2DS+filena[:-3]+'_IAT.png',dpi=200)
        plt.close(fig)
    
    Seconds_CH0[IAT_CH0<IAT_treshold] =np.nan
    Seconds_CH1[IAT_CH1<IAT_treshold] =np.nan
    
    #Remove nans again
    MeanXY_BG_CH0 = MeanXY_BG_CH0[~np.isnan(Seconds_CH0)]
    MaxDiameterBG_Ch0 = MaxDiameterBG_Ch0[~np.isnan(Seconds_CH0)]
    Edge_CH0 = Edge_CH0[~np.isnan(Seconds_CH0)]
    MeanXY_BG_CH1 = MeanXY_BG_CH1[~np.isnan(Seconds_CH1)]
    MaxDiameterBG_Ch1 = MaxDiameterBG_Ch1[~np.isnan(Seconds_CH1)]
    Edge_CH1 = Edge_CH1[~np.isnan(Seconds_CH1)]       
    IAT_CH0= IAT_CH0[~np.isnan(Seconds_CH0)]
    IAT_CH1= IAT_CH1[~np.isnan(Seconds_CH1)]
    ImageID_Ch0 = ImageID_Ch0[~np.isnan(Seconds_CH0)] 
    ImageID_Ch1 = ImageID_Ch1[~np.isnan(Seconds_CH1)] 
    ParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[~np.isnan(Seconds_CH0)]
    ParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[~np.isnan(Seconds_CH1)]
    
    Seconds_CH0= Seconds_CH0[~np.isnan(Seconds_CH0)]
    Seconds_CH1= Seconds_CH1[~np.isnan(Seconds_CH1)]
    
    # Look for colocated particles. 
    ChTimeDelta_high= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX_high= np.zeros(len(Seconds_CH1))*np.nan
    ChTimeDelta_low= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX_low= np.zeros(len(Seconds_CH1))*np.nan
    ChTimeDelta= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX= np.zeros(len(Seconds_CH1))*np.nan
    
#    # Go element by element in channel1. looking for the closest time match in channel0
#    for i in range(len(Seconds_CH1)) : 
#        ChIDX_high[i] = np.searchsorted(Seconds_CH0, Seconds_CH1[i], side="left")
#        if (ChIDX_high[i]<len(Seconds_CH0)):
#            ChTimeDelta_high[i]=np.absolute(Seconds_CH1[i] - Seconds_CH0[int(ChIDX_high[i])])
#    #check whether i or i-1 gives lowest time delta
#    ChIDX_low = ChIDX_high - 1
#    ChIDX_low[ChIDX_low<0] = 0
#    ChTimeDelta_low = np.array([np.absolute(Seconds_CH1[x] - Seconds_CH0[int(ChIDX_low[x])]) for x in range(len(ChIDX_low))])
#    ChIDX = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChIDX_low, ChIDX_high))
#    ChTimeDelta = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChTimeDelta_low, ChTimeDelta_high))

    # Go element by element in channel1. looking for the closest time match in channel0
    ChIDX_high = np.searchsorted(Seconds_CH0, Seconds_CH1, side="left") # side="left" means	Seconds_CH0[i-1] < Seconds_CH1 <= Seconds_CH0[i]
    ChIDX_high[ChIDX_high>=len(Seconds_CH0)] = len(Seconds_CH0)-1 # if outside the array
    #check whether i or i-1 gives lowest TimeDelta
    ChTimeDelta_high=np.absolute(Seconds_CH1 - Seconds_CH0[ChIDX_high])
    ChIDX_low = ChIDX_high - 1
    ChIDX_low[ChIDX_low<0] = 0
    ChTimeDelta_low = np.absolute(Seconds_CH1 - Seconds_CH0[ChIDX_low])
    ChIDX = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChIDX_low, ChIDX_high)) # Index with lowest TimeDelta
    ChTimeDelta = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChTimeDelta_low, ChTimeDelta_high)) # Lowest TimeDelta

    # Plot colocation time histogram
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 10})
        HistBins = np.logspace(-8,2,num=100)
        plt.hist(ChTimeDelta,HistBins)
        plt.xscale('log')
        plt.xlabel('Colocation time, s')
        plt.ylabel('Counts')
        plt.title('zeros = '+ str(len(ChTimeDelta[ChTimeDelta == 0]))+', IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
        plt.savefig(Path2DS+filena[:-3]+'_deltat.png',dpi=200)
        plt.close(fig)

    #Select colocated particle stats
    ColocationIDX = (ChIDX[ChTimeDelta <= ColocationThreshold]).astype(int) # Indexes are for channel 0 (same length as channel 1)
    ColocationDelta = ChTimeDelta[ChTimeDelta < ColocationThreshold] 
    #ColocationParticleTime_CH1 = ParticleTime_CH1[ChTimeDelta < ColocationThreshold]
    ColocationMeanXY_BG_CH1 = MeanXY_BG_CH1[ChTimeDelta < ColocationThreshold]
    ColocationMaxDiameterBG_Ch1 = MaxDiameterBG_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationSecondsCH1 = Seconds_CH1[ChTimeDelta < ColocationThreshold]
    ColocationEdgeCH1 = Edge_CH1[ChTimeDelta < ColocationThreshold]
    ColocationImageID_Ch1 = ImageID_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationParticleBufferTimeS_Ch1 = ParticleBufferTimeS_Ch1[ChTimeDelta < ColocationThreshold]
    
    ColocationSecondsCH0 = Seconds_CH0[ColocationIDX]
    ColocationMeanXY_BG_CH0 = MeanXY_BG_CH0[ColocationIDX]
    ColocationMaxDiameterBG_Ch0 = MaxDiameterBG_Ch0[ColocationIDX]
    ColocationEdgeCH0 = Edge_CH0[ColocationIDX]
    ColocationImageID_Ch0 = ImageID_Ch0[ColocationIDX]
    ColocationParticleBufferTimeS_Ch0 = ParticleBufferTimeS_Ch0[ColocationIDX]
    
    # select inter arrival time for only colocated particles 
    IAT_CH0_colocation = IAT_CH0[ColocationIDX]
    IAT_CH1_colocation = IAT_CH1[ChTimeDelta < ColocationThreshold]
    IAT_colocation = np.minimum(IAT_CH0_colocation,IAT_CH1_colocation)
    
    if 1 == 2:
        FileName = 'Export_'+filena
        for i in range(len(ColocationImageID_Ch0)) :
            if (ColocationEdgeCH0[i] == 2 and ColocationEdgeCH1[i] == 2): 
                if (i % 10 == 0) :
                    #print(i)
                    PlotImages_2Channels(Path2DS,FileName,ColocationParticleBufferTimeS_Ch0[i],ColocationImageID_Ch0[i],ColocationParticleBufferTimeS_Ch1[i],ColocationImageID_Ch1[i], Path2DS+'/ColocationImages/')
                if (i % 1000 == 0) :
                    print(i)

       
    # Plot CH1 vs CH0 diameter
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 12})
        plt.plot(ColocationMaxDiameterBG_Ch0,ColocationMaxDiameterBG_Ch1,'o',color='tab:gray',markersize=1)        
        plt.xlabel('Channel 0 max dimension, μm')
        plt.ylabel('Channel 1 max dimension, μm')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim([5,1000])
        plt.ylim([5,1000])
        plt.title('IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
        plt.savefig(Path2DS+filena[:-3]+'_colocationCH0vsCH1.png',dpi=200)
        plt.close(fig)
    
    # Plot IAT between colocated particles and all other particles
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 12})
        HistBins = np.logspace(-8,-2,num=100)
        IAT_colocation_hist, tmp= np.histogram(IAT_colocation,bins=HistBins)
        HistBinsMid = (HistBins[:-1] + HistBins[1:]) / 2
        plt.plot(HistBinsMid,IAT_colocation_hist)       
        plt.xscale('log')
        plt.xlabel('IAT (Colocated particles), s')
        plt.ylabel('Counts')
        plt.title('IAT = '+str(IAT_treshold)+'s, Colocation = ' + str(ColocationThreshold)+ 's' )
        plt.savefig(Path2DS+filena[:-3]+'_colocationIAT.png',dpi=200)
        plt.close(fig)
    
    # save matching particles
    if SaveFile == 1 : 
        SavePath= Path2DS+'Colocate_'+filena
        file = h5py.File(SavePath, 'w')
        file.create_dataset('ColocationSecondsCH1', data=ColocationSecondsCH1)
        file.create_dataset('ColocationSecondsCH0', data=ColocationSecondsCH0)
        file.create_dataset('ColocationDelta', data=ColocationDelta)
        file.create_dataset('ColocationMeanXY_BG_CH1', data=ColocationMeanXY_BG_CH1)
        file.create_dataset('ColocationMeanXY_BG_CH0', data=ColocationMeanXY_BG_CH0)
        file.create_dataset('ColocationMaxDiameterBG_Ch0', data=ColocationMaxDiameterBG_Ch0)
        file.create_dataset('ColocationMaxDiameterBG_Ch1', data=ColocationMaxDiameterBG_Ch1)
        file.create_dataset('ColocationEdgeCH0', data=ColocationEdgeCH0)
        file.create_dataset('ColocationEdgeCH1', data=ColocationEdgeCH1)
        
        file.close()
        #Create PSDs using all data and just colocated data
        if MeanXYFlag == 1:
            PSD_Colocate_1hz(ZdMax, Seconds_CH0, MeanXY_BG_CH0,Edge_CH0, Seconds_CH1, MeanXY_BG_CH1,Edge_CH1,SaveFile,Path2DS,filena,0)
            PSD_Colocate_1hz(ZdMax, ColocationSecondsCH0, ColocationMeanXY_BG_CH0, ColocationEdgeCH0, ColocationSecondsCH1, ColocationMeanXY_BG_CH1,ColocationEdgeCH1,SaveFile,Path2DS,filena,1)
        else : # use max diameter
            PSD_Colocate_1hz(ZdMax, Seconds_CH0, MaxDiameterBG_Ch0,Edge_CH0, Seconds_CH1, MaxDiameterBG_Ch1,Edge_CH1,SaveFile,Path2DS,filena,0)
            PSD_Colocate_1hz(ZdMax, ColocationSecondsCH0, ColocationMaxDiameterBG_Ch0,ColocationEdgeCH0,ColocationSecondsCH1, ColocationMaxDiameterBG_Ch1,ColocationEdgeCH1,SaveFile,Path2DS,filena,1)


#__________________________________________________________________________________
# Creates a psd from particle by particle data between StartAvg and EndAvg.

      
def ColocatePSD(StartAvg, EndAvg, Path2DS, FileName):
#    Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/' 
#    FileName = 'Colocate_base180207155526.h5'
#    ExpDate= datetime.datetime(2018, 2, 7, 0, 0, 0)
#    PltTitle='C078_1602'
    #StartAvg=datetime.datetime(2018, 2, 7, 16, 2, 0)
    #EndAvg=datetime.datetime(2018, 2, 7, 16, 10, 0)  
    
    Data_h5 = h5py.File(Path2DS + FileName, 'r')              
    ColocationSecondsCH1=np.array(Data_h5['ColocationSecondsCH1'])
    ColocationSecondsCH0=np.array(Data_h5['ColocationSecondsCH0'])
    #ColocationMeanXY_BG_CH1=np.array(Data_h5['ColocationMeanXY_BG_CH1'])
    #ColocationMeanXY_BG_CH0=np.array(Data_h5['ColocationMeanXY_BG_CH0'])
    DiameterCH0 =np.array(Data_h5['ColocationMaxDiameterBG_Ch0'])
    DiameterCH1 =np.array(Data_h5['ColocationMaxDiameterBG_Ch1'])
    
    Data_h5.close()
     
    Startime=int(datetime_2_s_sincemidnight(StartAvg))
    Endtime=int(datetime_2_s_sincemidnight(EndAvg))
    Ns=Endtime-Startime
    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
    
    ProbeRes= 10
    ArrayElements = 128
    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
    PSD_SizeMid=SizeBins-ProbeRes*0.5
    TAS = 100 # m/s
    TAS_mm = TAS * 1000 # mm/s

    SizeCH0_idx = np.digitize(DiameterCH0, SizeBins)
    SizeCH1_idx = np.digitize(DiameterCH1, SizeBins)
    
    # Calculate concentration of individual particle   
    ArrayWidthCH0_mm = ((ArrayElements*ProbeRes) - (DiameterCH0))*0.001 # mm
    SVolCH0_mm3_s = TAS_mm * (ArrayWidthCH0_mm**2)
    SVolCH0_L_s = SVolCH0_mm3_s / 1E6
    SVolCH0_L= SVolCH0_L_s * Ns

    ArrayWidthCH1_mm = ((ArrayElements*ProbeRes) - (DiameterCH1))*0.001 # mm
    SVolCH1_mm3_s = TAS_mm * (ArrayWidthCH1_mm**2)
    SVolCH1_L_s = SVolCH1_mm3_s / 1E6
    SVolCH1_L= SVolCH1_L_s * Ns

    # Filter and bin data
    Counts_PSD_CH0_colocation = np.zeros(len(SizeBins))
    PSD_L_CH0_colocation = np.zeros(len(SizeBins)) 

    Counts_PSD_CH1_colocation = np.zeros(len(SizeBins))
    PSD_L_CH1_colocation = np.zeros(len(SizeBins)) 
        
    # for colocated data. Assumimg ColocationSecondsCH0 == ColocationSecondsCH1
    for i in range(len(ColocationSecondsCH0)):      
        if (ColocationSecondsCH0[i] > Startime) & (ColocationSecondsCH0[i] < Endtime) : 
            Counts_PSD_CH0_colocation[SizeCH0_idx[i]]+=1
            PSD_L_CH0_colocation[SizeCH0_idx[i]]+= 1/(SVolCH0_L[i])
            
            Counts_PSD_CH1_colocation[SizeCH1_idx[i]]+=1
            PSD_L_CH1_colocation[SizeCH1_idx[i]]+= 1/(SVolCH1_L[i])      
    
    dNdD_L_CH0_colocation = PSD_L_CH0_colocation/ProbeRes # normalise using bin width
    dNdD_L_CH0_colocation[0:2] =np.nan # nan first two bins
    dNdD_L_CH0_colocation_error = dNdD_L_CH0_colocation / (2*np.sqrt(Counts_PSD_CH0_colocation)) # /2 since plot plus minus

    dNdD_L_CH1_colocation = PSD_L_CH1_colocation/ProbeRes # normalise using bin width
    dNdD_L_CH1_colocation[0:2] =np.nan # nan first two bins
    dNdD_L_CH1_colocation_error = dNdD_L_CH1_colocation / (2*np.sqrt(Counts_PSD_CH1_colocation)) # /2 since plot plus minus

    return PSD_SizeMid, dNdD_L_CH0_colocation, dNdD_L_CH0_colocation_error, dNdD_L_CH1_colocation, dNdD_L_CH1_colocation_error

#__________________________________________________________________________________

#Diameter um
#Lambda um
# ArmSep mm
# ProbeRes um
#TAS m/s

# Sample volume from Lawson et 2006

def Lawson_sVol_array(Diameter,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS):

    radius=Diameter/2 # um
    DoF=(2*c*radius**2)/Lambda ## 2* because there is a +/- in the equation
    DoF/=1000 # return in mm
    
    DoF = np.minimum(DoF,ArmSep)
        
    ArrayWidth = ((ArrayElements*ProbeRes) - (Diameter))*0.001 # mm
    
    TAS*=1000 # mm/s
    
    sVol = DoF * ArrayWidth * TAS # mm/s
   
    sVol /= 1E6 #L/s
    
    return sVol#, DoF



#__________________________________________________________________________________
#Calculate sample volume for 2ds laser overlap

def Colocated_sVol_array(DiameterCH0, DiameterCH1,Lambda_um,ZdMax,ArmSep,ArrayElements,PixelSize,TAS):
    
         
    Zum_at_ZdMax_CH0 = (ZdMax * DiameterCH0**2) / (4 * Lambda_um)     
    MinElement = np.maximum(0, 63-np.round(Zum_at_ZdMax_CH0/PixelSize))
    MaxElement = np.minimum(127,63+np.round(Zum_at_ZdMax_CH0/PixelSize))
    ArrayWidthCH0_mm = ((MaxElement - MinElement + 1) * PixelSize -  DiameterCH0)*0.001 # mm  
    
    Zum_at_ZdMax_CH1 = (ZdMax * DiameterCH1**2) / (4 * Lambda_um)
    MinElement = np.maximum(0, 63-np.round(Zum_at_ZdMax_CH1/PixelSize))
    MaxElement = np.minimum(127,63+np.round(Zum_at_ZdMax_CH1/PixelSize))
    ArrayWidthCH1_mm = ((MaxElement - MinElement + 1) * PixelSize -  DiameterCH1)*0.001 # mm 
    
    TAS*=1000 # mm/s
    SVol_mm3_s = TAS * ArrayWidthCH0_mm * ArrayWidthCH1_mm
    SVol_L_s = SVol_mm3_s / 1E6
    
    return SVol_L_s

#__________________________________________________________________________________
 

#Create arrays with dNdD_L_CH0 and dNdD_L_CH1. 
#ColocatedFlag =0  Lawson et al sample volume 
#ColocatedFlag =1  colocation sample volume 

def PSD_Colocate_1hz(ZdMax,SecondsCH0, DiameterCH0, Edge_CH0, SecondsCH1, DiameterCH1,Edge_CH1,SaveFile,Path2DS,filena,ColocatedFlag) :
    ProbeRes= 10
    ArrayElements = 128
    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
    PSD_SizeMid=SizeBins-ProbeRes*0.5
    TAS = 100 # m/s
    c=6
    ArmSep = 63
    Lambda = 0.658

    Startime=int(np.minimum(SecondsCH0[0], SecondsCH1[0]))
    Endtime=int(np.maximum(SecondsCH0[-1], SecondsCH1[-1]))
    if Endtime > 172800:
        Endtime = 172800
    
    Ns=Endtime-Startime
    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
    TimeMidBins_s=TimeBins_s-0.5 

    TimeCH0_idx = np.digitize(SecondsCH0, TimeBins_s)
    TimeCH1_idx = np.digitize(SecondsCH1, TimeBins_s)
    SizeCH0_idx = np.digitize(DiameterCH0, SizeBins)
    SizeCH1_idx = np.digitize(DiameterCH1, SizeBins)
  
   
    
    # Calculate concentration of individual particle   
    if ColocatedFlag == 1 : 
#        TimeCH0_idx[Edge_CH0 == 2] = np.nan
#        TimeCH1_idx[Edge_CH1 == 2] = np.nan
#        SizeCH0_idx[Edge_CH0 == 2] = np.nan
#        SizeCH1_idx[Edge_CH1 == 2] = np.nan
        SVolCH0_L_s = Colocated_sVol_array(DiameterCH0, DiameterCH1,Lambda,ZdMax,ArmSep,ArrayElements,ProbeRes,TAS)
        SVolCH1_L_s = SVolCH0_L_s
    else : 
#        TimeCH0_idx[Edge_CH0 == 1] = np.nan
#        TimeCH1_idx[Edge_CH1 == 1] = np.nan
#        SizeCH0_idx[Edge_CH0 == 1] = np.nan
#        SizeCH1_idx[Edge_CH1 == 1] = np.nan
        SVolCH0_L_s = Lawson_sVol_array(DiameterCH0,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS)    
        SVolCH1_L_s = Lawson_sVol_array(DiameterCH1,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS) 

    # Filter and bin data
    Counts_PSD_CH0 = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L_CH0 = np.zeros((len(TimeBins_s), len(SizeBins)))

    Counts_PSD_CH1 = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L_CH1 = np.zeros((len(TimeBins_s), len(SizeBins)))  
    
    # Edge = 0 not touching array edge
    # Edge = 1 touching array edge
    # Edge = 2 Zd < 2
    
    for i in range(len(SecondsCH0)):      
        if (SecondsCH0[i] > Startime) & (SecondsCH0[i] < Endtime) : 
            if SizeCH0_idx[i] < ArrayElements:
                if (ColocatedFlag == 0 and Edge_CH0[i] != 1) or (ColocatedFlag == 1 and Edge_CH0[i] == 2 and Edge_CH1[i] == 2): # filter edge particles 
                    Counts_PSD_CH0[TimeCH0_idx[i],SizeCH0_idx[i]]+=1
                    PSD_L_CH0[TimeCH0_idx[i],SizeCH0_idx[i]]+= 1/(SVolCH0_L_s[i])
            
    for i in range(len(SecondsCH1)):      
        if (SecondsCH1[i] > Startime) & (SecondsCH1[i] < Endtime) : 
            if SizeCH1_idx[i] < ArrayElements:
                if (ColocatedFlag == 0 and Edge_CH1[i] != 1) or (ColocatedFlag == 1 and Edge_CH0[i] == 2 and Edge_CH1[i] == 2):  # filter edge particles            
                    Counts_PSD_CH1[TimeCH1_idx[i],SizeCH1_idx[i]]+=1
                    PSD_L_CH1[TimeCH1_idx[i],SizeCH1_idx[i]]+= 1/(SVolCH1_L_s[i])     
    
    dNdD_L_CH0 = PSD_L_CH0/ProbeRes # normalise using bin width
    dNdD_L_CH1 = PSD_L_CH1/ProbeRes # normalise using bin width
    
#    if ColocatedFlag == 1:
#        dNdD_L_CH0[:,0:2] =np.nan # nan first two bins
#        dNdD_L_CH1[:,0:2] =np.nan # nan first two bins

    if SaveFile == 1 : 
        
        if ColocatedFlag == 1 :
            SavePath= Path2DS+'dNdD_L_Colocate_'+filena
        else:
            SavePath= Path2DS+'dNdD_L_'+filena
        file = h5py.File(SavePath, 'w')
        file.create_dataset('dNdD_L_CH0', data=dNdD_L_CH0)
        file.create_dataset('Counts_PSD_CH0', data=Counts_PSD_CH0)
        file.create_dataset('dNdD_L_CH1', data=dNdD_L_CH1)
        file.create_dataset('Counts_PSD_CH1', data=Counts_PSD_CH1)
        file.create_dataset('TimeMidBins_s', data=TimeMidBins_s)
        file.create_dataset('PSD_SizeMid', data=PSD_SizeMid)
        file.close()


#__________________________________________________________________________________
    
##Create arrays with dNdD_L_CH0_colocation vs time 
#
#def ColocatePSD_1hz(Path2DS, filena):
##    Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/' 
##    FileName = 'Colocate_base180207155526.h5'
##    ExpDate= datetime.datetime(2018, 2, 7, 0, 0, 0)
##    PltTitle='C078_1602'
#    #StartAvg=datetime.datetime(2018, 2, 7, 16, 2, 0)
#    #EndAvg=datetime.datetime(2018, 2, 7, 16, 10, 0)  
#    
#    SaveFile = 1
#    
#    Data_h5 = h5py.File(Path2DS + filena, 'r')              
#    #ColocationSecondsCH1=np.array(Data_h5['ColocationSecondsCH1'])
#    ColocationSecondsCH0=np.array(Data_h5['ColocationSecondsCH0'])
#    ColocationMeanXY_BG_CH1=np.array(Data_h5['ColocationMeanXY_BG_CH1'])
#    ColocationMeanXY_BG_CH0=np.array(Data_h5['ColocationMeanXY_BG_CH0'])
#    Data_h5.close()
#     
#    Startime=int(ColocationSecondsCH0[0])
#    Endtime=int(ColocationSecondsCH0[-1])
#    Ns=Endtime-Startime
#    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
#    TimeMidBins_s=TimeBins_s-0.5
#    
#    ProbeRes= 10
#    ArrayElements = 128
#    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
#    PSD_SizeMid=SizeBins-ProbeRes*0.5
#    TAS = 100 # m/s
#    TAS_mm = TAS * 1000 # mm/s
#
#    Time_idx = np.digitize(ColocationSecondsCH0, TimeBins_s)
#    SizeCH0_idx = np.digitize(ColocationMeanXY_BG_CH0, SizeBins)
#    SizeCH1_idx = np.digitize(ColocationMeanXY_BG_CH1, SizeBins)
#    
#    # Calculate concentration of individual particle   
#    ArrayWidthCH0_mm = ((ArrayElements*ProbeRes) - (ColocationMeanXY_BG_CH0))*0.001 # mm
#    ArrayWidthCH1_mm = ((ArrayElements*ProbeRes) - (ColocationMeanXY_BG_CH1))*0.001 # mm
#    
#    SVolCH0_mm3_s = TAS_mm * ArrayWidthCH0_mm * ArrayWidthCH1_mm
#    SVolCH0_L_s = SVolCH0_mm3_s / 1E6 # 1 second sample volume
#    #SVolCH0_L= SVolCH0_L_s * Ns
#
#    SVolCH1_mm3_s = TAS_mm * ArrayWidthCH0_mm * ArrayWidthCH1_mm
#    SVolCH1_L_s = SVolCH1_mm3_s / 1E6 # 1 second sample volume
#    #SVolCH1_L= SVolCH1_L_s * Ns
#
#    # Filter and bin data
#    Counts_PSD_CH0_colocation = np.zeros((len(TimeBins_s), len(SizeBins)))
#    PSD_L_CH0_colocation = np.zeros((len(TimeBins_s), len(SizeBins)))
#
#    Counts_PSD_CH1_colocation = np.zeros((len(TimeBins_s), len(SizeBins)))
#    PSD_L_CH1_colocation = np.zeros((len(TimeBins_s), len(SizeBins)))
#        
#    # for colocated data. Assumimg ColocationSecondsCH0 == ColocationSecondsCH1
#    for i in range(len(ColocationSecondsCH0)):      
#        if (ColocationSecondsCH0[i] > Startime) & (ColocationSecondsCH0[i] < Endtime) : 
#            Counts_PSD_CH0_colocation[Time_idx[i],SizeCH0_idx[i]]+=1
#            PSD_L_CH0_colocation[Time_idx[i],SizeCH0_idx[i]]+= 1/(SVolCH0_L_s[i])
#            
#            Counts_PSD_CH1_colocation[Time_idx[i],SizeCH1_idx[i]]+=1
#            PSD_L_CH1_colocation[Time_idx[i],SizeCH1_idx[i]]+= 1/(SVolCH1_L_s[i])      
#    
#    dNdD_L_CH0_colocation = PSD_L_CH0_colocation/ProbeRes # normalise using bin width
#    dNdD_L_CH0_colocation[:,0:2] =np.nan # nan first two bins
#    dNdD_L_CH0_colocation_error = dNdD_L_CH0_colocation / (2*np.sqrt(Counts_PSD_CH0_colocation)) # /2 since plot plus minus
#
#    dNdD_L_CH1_colocation = PSD_L_CH1_colocation/ProbeRes # normalise using bin width
#    dNdD_L_CH1_colocation[:,0:2] =np.nan # nan first two bins
#    dNdD_L_CH1_colocation_error = dNdD_L_CH1_colocation / (2*np.sqrt(Counts_PSD_CH1_colocation)) # /2 since plot plus minus
#
#    if SaveFile == 1 : 
#        SavePath= Path2DS+'dNdD_L_'+filena
#        file = h5py.File(SavePath, 'w')
#        file.create_dataset('dNdD_L_CH0_colocation', data=dNdD_L_CH0_colocation)
#        file.create_dataset('dNdD_L_CH0_colocation_error', data=dNdD_L_CH0_colocation_error)
#        file.create_dataset('dNdD_L_CH1_colocation', data=dNdD_L_CH1_colocation)
#        file.create_dataset('dNdD_L_CH1_colocation_error', data=dNdD_L_CH1_colocation_error)
#        file.create_dataset('TimeMidBins_s', data=TimeMidBins_s)
#        file.create_dataset('PSD_SizeMid', data=PSD_SizeMid)
#        file.close()
#__________________________________________________________________________________
# CIRCCREX flight B895 create size distributions for each run. Plot and save data   

     
def ColocateB895():
    
    SavePlot = 1
    SaveData = 1
    
    FlightDate = datetime.datetime(2015, 3, 13, 0, 0, 0)
    RunStart = [datetime.datetime(2015, 3, 13, 10, 37,33),datetime.datetime(2015, 3, 13, 10, 42,31),datetime.datetime(2015, 3, 13, 10, 47,40),datetime.datetime(2015, 3, 13, 10, 51,26),datetime.datetime(2015, 3, 13, 11, 1,2),datetime.datetime(2015, 3, 13, 11, 4,57),datetime.datetime(2015, 3, 13, 11, 22,1)]
    RunEnd = [datetime.datetime(2015, 3, 13, 10, 42,31), datetime.datetime(2015, 3, 13, 10, 47,40), datetime.datetime(2015, 3, 13, 10, 51,26), datetime.datetime(2015, 3, 13, 10, 56,21), datetime.datetime(2015, 3, 13, 11, 4,57), datetime.datetime(2015, 3, 13, 11, 10,16), datetime.datetime(2015, 3, 13, 11, 26,7)]
    RunTemperature = [226,229,232,234,237,239,245]
    Path2DS = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/oasisoutput/Colocation/'
    FileName = 'base150313084741.h5'
    
    Data_h5 = h5py.File(Path2DS + 'dNdD_L_Colocate_'+FileName, 'r')
    dNdD_L_CH0_colocation=np.array(Data_h5['dNdD_L_CH0'])
    dNdD_L_CH1_colocation=np.array(Data_h5['dNdD_L_CH1'])
    Counts_PSD_CH0_colocation=np.array(Data_h5['Counts_PSD_CH0'])
    Counts_PSD_CH1_colocation=np.array(Data_h5['Counts_PSD_CH1'])
    TimeMidBins_s_colocation=np.array(Data_h5['TimeMidBins_s'])
    PSD_SizeMid_colocation=np.array(Data_h5['PSD_SizeMid'])
    Data_h5.close()
    #Load all PSDs
    Data_h5 = h5py.File(Path2DS + 'dNdD_L_'+FileName, 'r')
    dNdD_L_CH0=np.array(Data_h5['dNdD_L_CH0'])
    dNdD_L_CH1=np.array(Data_h5['dNdD_L_CH1'])
    Counts_PSD_CH0=np.array(Data_h5['Counts_PSD_CH0'])
    Counts_PSD_CH1=np.array(Data_h5['Counts_PSD_CH1'])
    TimeMidBins_s=np.array(Data_h5['TimeMidBins_s'])
    PSD_SizeMid=np.array(Data_h5['PSD_SizeMid'])
    Data_h5.close()
    
    #Set up time to average psd over
    #AvgBinsMid = AvgBinsEdge[0:-2:1] + AvgBinsEdge[1:-1:1]
    
    Start_s= [(RunStart[x]- FlightDate).total_seconds() for x in range(len(RunStart))]  
    End_s= [(RunEnd[x]- FlightDate).total_seconds() for x in range(len(RunEnd))]
    
    IDX_Start_colocation = np.searchsorted(TimeMidBins_s_colocation,Start_s)
    IDX_End_colocation = np.searchsorted(TimeMidBins_s_colocation,End_s)
    
    IDX_Start = np.searchsorted(TimeMidBins_s,Start_s)
    IDX_End = np.searchsorted(TimeMidBins_s,End_s)
    
    #Load HALOHolo
    PathHALO = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/HALAholoSizeDistros_v3/'
    HALOName= 'RunsHALOnnd.csv'
    HaloHoloData =pd.read_csv(PathHALO+HALOName)
    HaloholoMidsize = HaloHoloData['MidSize']
    HaloholoMidsize= np.array(HaloholoMidsize[1:].astype(np.float))
    
    CIP100Name = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/CIP100/20150313085414/OASIS_OUTPUT/'
    CIP100Data = 'dNdD_L_CIP100_maxD.h5'
    Data_h5 = h5py.File(CIP100Name+CIP100Data, 'r')
    dNdD_L_CIP100=np.array(Data_h5['dNdD_L'])
    Counts_PSD_CIP100=np.array(Data_h5['Counts_PSD'])
    TimeMidBins_s_CIP100=np.array(Data_h5['TimeMidBins_s'])
    PSD_SizeMid_CIP100=np.array(Data_h5['PSD_SizeMid'])
    Data_h5.close()
    
    IDX_Start_CIP100 = np.searchsorted(TimeMidBins_s_CIP100,Start_s)
    IDX_End_CIP100 = np.searchsorted(TimeMidBins_s_CIP100,End_s)
      
#    OASISPath = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/'
#    OASISName = 'B895_OASIS_2ds.csv'
#    OASISData = pd.read_csv(OASISPath+OASISName)
#    OASISMidsize = OASISData['MidSize']
#    OASISMidsize= np.array(OASISMidsize[1:].astype(np.float))
    
    fig=plt.figure(figsize=(15,15)) 
    plt.rcParams.update({'font.size': 14})
    
    Run_dNdD_L_CH0_colocation_avg= np.zeros([128,len(IDX_Start_colocation)])*np.nan
    Run_dNdD_L_CH1_colocation_avg= np.zeros([128,len(IDX_Start_colocation)])*np.nan
    Run_dNdD_L_CH0_avg= np.zeros([128,len(IDX_Start_colocation)])*np.nan
    Run_dNdD_L_CH1_avg= np.zeros([128,len(IDX_Start_colocation)])*np.nan
    
    for i in range(len(IDX_Start_colocation)):
        
        plt.subplot(3, 3, i+1)
        #fig=plt.figure(figsize=(7,7)) 
        StartIdx = IDX_Start_CIP100[i]
        EndIdx = IDX_End_CIP100[i]
        dNdD_L_CIP100_avg = np.nanmean(dNdD_L_CIP100[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CIP100_avg[0] =np.nan # Ditch first bin
        dNdD_L_CIP100_error = dNdD_L_CIP100_avg / (np.sqrt(np.nansum(Counts_PSD_CIP100[StartIdx:EndIdx][:],axis=0)))       
        
        
        StartIdx = IDX_Start[i]
        EndIdx = IDX_End[i] 
        dNdD_L_CH0_avg = np.nanmean(dNdD_L_CH0[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_avg = np.nanmean(dNdD_L_CH1[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_error = dNdD_L_CH0_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH0[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_error = dNdD_L_CH1_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH1[StartIdx:EndIdx][:],axis=0)))
        #dNdD_L_CH0_avg[dNdD_L_CH0_error == 0] = np.nan
        #dNdD_L_CH1_avg[dNdD_L_CH1_error == 0] = np.nan
        Run_dNdD_L_CH0_avg[:,i]= dNdD_L_CH0_avg
        Run_dNdD_L_CH1_avg[:,i]= dNdD_L_CH1_avg
        
        
        StartIdx = IDX_Start_colocation[i]
        EndIdx = IDX_End_colocation[i]
        dNdD_L_CH0_colocation_avg = np.nanmean(dNdD_L_CH0_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_colocation_avg = np.nanmean(dNdD_L_CH1_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_colocation_error = dNdD_L_CH0_colocation_avg / ( np.sqrt(np.nansum(Counts_PSD_CH0_colocation[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_colocation_error = dNdD_L_CH1_colocation_avg / ( np.sqrt(np.nansum(Counts_PSD_CH1_colocation[StartIdx:EndIdx][:],axis=0)))
        #dNdD_L_CH0_colocation_avg[dNdD_L_CH0_colocation_error==0] = np.nan
        #dNdD_L_CH1_colocation_avg[dNdD_L_CH1_colocation_error==0] = np.nan

        Run_dNdD_L_CH0_colocation_avg[:,i]= dNdD_L_CH0_colocation_avg
        Run_dNdD_L_CH1_colocation_avg[:,i]= dNdD_L_CH1_colocation_avg

        HaloRun = HaloHoloData[str(RunTemperature[i])+'K']
        HaloHoloPSD = np.array(HaloRun[1:].astype(np.float))
        HaloHoloPSD[0:6] = np.nan
        
        L1 = plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH0_avg,yerr=(dNdD_L_CH0_error), color = 'tab:gray', label = '2DS all')
        plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH1_avg,yerr=(dNdD_L_CH1_error), color = 'tab:gray')
        L2 = plt.plot(PSD_SizeMid_CIP100, dNdD_L_CIP100_avg,color='k', label= 'CIP100')
        L3 = plt.errorbar(PSD_SizeMid,dNdD_L_CH0_colocation_avg,yerr=(dNdD_L_CH0_colocation_error), color = 'r', fmt='o',markersize=5, capsize=2, label = 'CH0 colocation')
        L4 = plt.errorbar(PSD_SizeMid,dNdD_L_CH1_colocation_avg,yerr=(dNdD_L_CH1_colocation_error), color = 'c', fmt='o',markersize=5, capsize=2, label = 'CH1 colocation')
        L5 = plt.plot(HaloholoMidsize, HaloHoloPSD,color='b', label='HALOHolo')
        
#        OASISRun = OASISData[str(RunTemperature[i])+'K']
#        OASISPSD = np.array(OASISRun[1:].astype(np.float))
#        plt.plot(OASISMidsize, OASISPSD, label='OASIS 2DS')             
        
        plt.title(str(RunTemperature[i])+' K')
        #plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
        #plt.xlabel('Diameter, μm')      
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim([8,1000])
        plt.ylim([1E-4,1E1])
        #plt.legend()
        #plt.savefig(Path2DS+'Psd'+str(RunTemperature[i])+'K'+'.png',dpi=200)
        #plt.close(fig)
       
#        plt.plot(PSD_SizeMid,dNdD_L_CH0_colocation)
#        plt.plot(PSD_SizeMid,dNdD_L_CH1_colocation)
#        
#        HaloRun = HaloHoloData[str(RunTemperature[i])+'K']
#        HaloHoloPSD = HaloRun[1:].astype(np.float)
#        HaloHoloPSD[0:7] = np.nan
#        plt.plot(HaloholoMidsize, HaloHoloPSD)
#        
        plt.title(str(RunTemperature[i])+' K')
        #if ((i+1) % 2 != 0) :  # odd panel
        if (i+1 == 1 or i+1 == 4 or i+1 == 7) :
            plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
        #if (i+ 1 == 7 ) or (i + 1 == 6) :
        if (i+1 >= 5):
            plt.xlabel('Diameter, μm')      
#        plt.yscale('log')
#        plt.xscale('log')
#        plt.xlim([10,1000])
#        plt.ylim([1E-3,1])
    
    #plt.subplot(3, 3, 8)      
    #plt.legend((L1, L2, L3, L4, L5), ('2DS', 'CIP100','CH0','CH1','HALOHolo'))
    plt.legend(bbox_to_anchor=(1.3, 0.5), loc='center left', borderaxespad=0.)
    
    
    if (SavePlot == 1) : 
        plt.savefig(Path2DS+FileName[:-3]+'_PSDs.png',dpi=200)
        plt.close(fig)
    
    
    if (SaveData == 1) : 
    
        SavePath= Path2DS+'B895_dNdD_L_runavg.h5'
        file = h5py.File(SavePath, 'w')
        
        file.create_dataset('Run_dNdD_L_CH0_colocation_avg', data=Run_dNdD_L_CH0_colocation_avg)
        file.create_dataset('Run_dNdD_L_CH1_colocation_avg', data=Run_dNdD_L_CH1_colocation_avg)
        file.create_dataset('Run_dNdD_L_CH0_avg', data=Run_dNdD_L_CH0_avg)
        file.create_dataset('Run_dNdD_L_CH1_avg', data=Run_dNdD_L_CH1_avg)
        file.create_dataset('PSD_SizeMid_colocation', data=PSD_SizeMid_colocation)
        file.close()
    
    #return Run_dNdD_L_CH0_colocation_avg

#__________________________________________________________________________________

def ColocateB895_Field07(): 

    Path2DS = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/oasisoutput/Colocation/'
    SavePath= Path2DS+'B895_dNdD_L_runavg.h5'
    
    Data_h5 = h5py.File(SavePath, 'r')
    Run_dNdD_L_CH0_colocation_avg=np.array(Data_h5['Run_dNdD_L_CH0_colocation_avg'])
    Run_dNdD_L_CH1_colocation_avg=np.array(Data_h5['Run_dNdD_L_CH1_colocation_avg'])
    Run_dNdD_L_CH0_avg=np.array(Data_h5['Run_dNdD_L_CH0_avg'])
    Run_dNdD_L_CH1_avg=np.array(Data_h5['Run_dNdD_L_CH1_avg'])
    PSD_SizeMid=np.array(Data_h5['PSD_SizeMid_colocation'])
    Data_h5.close()

    PSD_SizeWidth = np.ones(len(PSD_SizeMid)) * 10
    
    fig=plt.figure(figsize=(8,8)) 
    plt.rcParams.update({'font.size': 12})

    #Run_dNdD_L_colocation = 
    for i in range(len(Run_dNdD_L_CH0_colocation_avg[0,:])) : 
        dNdD_L_CH0_colocation_avg = Run_dNdD_L_CH0_colocation_avg[:,i]
        PSD_y = dNdD_L_CH0_colocation_avg * 1000 # Number densitys	(#/m3 /um)
        
        PSD_y_norm, PSD_xMid_norm = Field07_PSD_Normalise(PSD_SizeMid, PSD_SizeWidth, PSD_y)
        if i == 0 : 
            plt.plot(PSD_xMid_norm, PSD_y_norm,'o',color = 'tab:gray',label='2DS colocation')
        else :
            plt.plot(PSD_xMid_norm, PSD_y_norm,'o',color = 'tab:gray')
    
    Field_Trop, Field_Mid, NewPSD = Field07_PSD(PSD_xMid_norm)
    plt.plot(PSD_xMid_norm, Field_Trop, color='r', label = 'Tropical (Field et al, 2007)')
    plt.plot(PSD_xMid_norm, Field_Mid, color='b', label = 'Mid-latitude (Field et al, 2007)')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('DM$_{2}$ / M$_{3}$')
    plt.ylabel('N(D)M$_{3}^{3}$ / M$_{2}^{4}$')
    # PSD_xMid		// x bin Midpoints		(um)
    # PSD_xWidth	// x bin Widths		(um)
    # PSD_y			// Number densitys	(#/m3 /um)
    # NormX		// flag to determine if X data should be normalised. 1=normalise. 0=no normalisation

    plt.savefig(Path2DS+'b895_Colocation_field07.png',dpi=200)
    plt.close(fig)



#__________________________________________________________________________________

# CIRCCREX flight B895 create size distributions averaged over DeltaT. Plot PSDs

def ColocateC031():
    
    ExpDate= datetime.datetime(2017, 8, 18, 0, 0, 0)
    Path2DS = 'C:/Users/Admin TEMP/Documents/Clarify/C031/'
    FileName = 'base170818115058.h5'
    # Load colocation PSDs
    Data_h5 = h5py.File(Path2DS + 'dNdD_L_Colocate_'+FileName, 'r')
    dNdD_L_CH0_colocation=np.array(Data_h5['dNdD_L_CH0'])
    dNdD_L_CH1_colocation=np.array(Data_h5['dNdD_L_CH1'])
    Counts_PSD_CH0_colocation=np.array(Data_h5['Counts_PSD_CH0'])
    Counts_PSD_CH1_colocation=np.array(Data_h5['Counts_PSD_CH1'])
    TimeMidBins_s_colocation=np.array(Data_h5['TimeMidBins_s'])
    PSD_SizeMid_colocation=np.array(Data_h5['PSD_SizeMid'])
    Data_h5.close()
    #Load all PSDs
    Data_h5 = h5py.File(Path2DS + 'dNdD_L_'+FileName, 'r')
    dNdD_L_CH0=np.array(Data_h5['dNdD_L_CH0'])
    dNdD_L_CH1=np.array(Data_h5['dNdD_L_CH1'])
    Counts_PSD_CH0=np.array(Data_h5['Counts_PSD_CH0'])
    Counts_PSD_CH1=np.array(Data_h5['Counts_PSD_CH1'])
    TimeMidBins_s=np.array(Data_h5['TimeMidBins_s'])
    PSD_SizeMid=np.array(Data_h5['PSD_SizeMid'])
    Data_h5.close()
    
    #Set up time to average psd over
    DeltaT = 300
    AvgBinsEdge = np.arange(TimeMidBins_s_colocation[0], TimeMidBins_s_colocation[-1],DeltaT)
    AvgBinsMid = AvgBinsEdge[0:-2:1] + AvgBinsEdge[1:-1:1]
    
    IDX_colocation = np.searchsorted(TimeMidBins_s_colocation,AvgBinsEdge)
    IDX = np.searchsorted(TimeMidBins_s,AvgBinsEdge)
    
    
    CoreCloudFile='core-cloud-phy_faam_20170818_v501_r0_c031_pcasp1.nc'
    CoreCloudPath='C:/Users/Admin TEMP/Documents/Clarify/C031/'
    CdpCalPath='C:/Users/Admin TEMP/Documents/Clarify/'
    CdpCalFile='IAN_CDPbins_formatted.csv'   
    
    CDPBinsEdge = [(ExpDate+datetime.timedelta(seconds=AvgBinsEdge[x])) for x in range(len(AvgBinsEdge))]       
    
    for i in range(len(AvgBinsMid)-1) : 
        StartIdx = IDX[i]
        EndIdx = IDX[i+1] 
        dNdD_L_CH0_avg = np.nanmean(dNdD_L_CH0[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_avg = np.nanmean(dNdD_L_CH1[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_error = dNdD_L_CH0_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH0[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_error = dNdD_L_CH1_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH1[StartIdx:EndIdx][:],axis=0)))
        
        StartIdx = IDX_colocation[i]
        EndIdx = IDX_colocation[i+1]
        dNdD_L_CH0_colocation_avg = np.nanmean(dNdD_L_CH0_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_colocation_avg = np.nanmean(dNdD_L_CH1_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_colocation_error = dNdD_L_CH0_colocation_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH0_colocation[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_colocation_error = dNdD_L_CH1_colocation_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH1_colocation[StartIdx:EndIdx][:],axis=0)))

        CDP_BinCentre, CDP_dNdDp_L_avg = Avg_CDP_dndDp(CDPBinsEdge[i],CDPBinsEdge[i+1],CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)

        fig=plt.figure(figsize=(7,7)) 
        
        plt.plot(CDP_BinCentre, CDP_dNdDp_L_avg,color='b', label='CDP')
        plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH0_avg,yerr=(dNdD_L_CH0_error), color = 'tab:gray', label = '2DS all')
        plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH1_avg,yerr=(dNdD_L_CH1_error), color = 'tab:gray')
        plt.errorbar(PSD_SizeMid,dNdD_L_CH0_colocation_avg,yerr=(dNdD_L_CH0_colocation_error), color = 'r', fmt='o',markersize=5, capsize=2, label = 'CH0 colocation')
        plt.errorbar(PSD_SizeMid,dNdD_L_CH1_colocation_avg,yerr=(dNdD_L_CH1_colocation_error), color = 'c', fmt='o',markersize=5, capsize=2, label = 'CH1 colocation')
        
        plt.title(str(AvgBinsEdge[i])+'s to '+str(AvgBinsEdge[i+1]))
        plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
        plt.xlabel('Diameter, μm')      
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim([1,1000])
        plt.ylim([1E-5,1E4])
        plt.legend()
        plt.savefig(Path2DS+FileName[:-3]+'_psd'+str(int(AvgBinsMid[i]))+'s.png',dpi=200)
        plt.close(fig)
        
    
#__________________________________________________________________________________

# Take OAP h5 output from oasis and determed PSDs at 1HZ

def OAP_H5toPSD():

    PathOAP = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/CIP100/20150313085414/OASIS_OUTPUT/'
    Prefix = 'Imagefile2_'
    Suffix = ".h5"
    OutputName = 'CIP100_maxD.h5'

    for filena in os.listdir(PathOAP):
        if filena.endswith(Suffix) and filena.startswith(Prefix):
            print(filena)
            TmpMeanXY_BG, TmpDiameterBGy, TmpDiameterBGx, TmpMaxDiameterBG, TmpParticleTimeSeconds = LoadFilterOAPh5(PathOAP,filena)
        
            if 'MeanXY_BG' not in dir():
                MeanXY_BG = TmpMeanXY_BG
                DiameterBGy = TmpDiameterBGy
                DiameterBGx = TmpDiameterBGx
                MaxDiameterBG = TmpMaxDiameterBG
                ParticleTimeSeconds = TmpParticleTimeSeconds
            else:    
                MeanXY_BG= np.append(MeanXY_BG,TmpMeanXY_BG,axis=0)
                DiameterBGy= np.append(DiameterBGy,TmpDiameterBGy,axis=0)
                DiameterBGx= np.append(DiameterBGx,TmpDiameterBGx,axis=0)
                MaxDiameterBG= np.append(MaxDiameterBG,TmpMaxDiameterBG,axis=0)
                ParticleTimeSeconds= np.append(ParticleTimeSeconds,TmpParticleTimeSeconds,axis=0)
                
    PSD_1hz(ParticleTimeSeconds, MaxDiameterBG,1,PathOAP,OutputName)

#__________________________________________________________________________________
 

#Go from PbP seconds and diameters to dNdD_L at 1Hz

def PSD_1hz(Seconds, Diameter,SaveFile,PathOAP,SaveName):
    
    ProbeRes= 100
    ArrayElements = 64
    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
    PSD_SizeMid=SizeBins-ProbeRes*0.5
    TAS = 100 # m/s
    c=6
    ArmSep = 70
    Lambda = 0.658

    Startime=int(Seconds[0])
    Endtime=int(np.minimum(Seconds[-1], 172800))
    
    Ns=Endtime-Startime
    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
    TimeMidBins_s=TimeBins_s-0.5 

    Time_idx = np.digitize(Seconds, TimeBins_s)
    Size_idx = np.digitize(Diameter, SizeBins)
    
    # Calculate concentration of individual particle   
    SVol_L_s = Lawson_sVol_array(Diameter,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS)    

    # Filter and bin data
    Counts_PSD = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L = np.zeros((len(TimeBins_s), len(SizeBins)))
        
    for i in range(len(Seconds)):      
        if (Seconds[i] > Startime) & (Seconds[i] < Endtime) : 
            if Size_idx[i] < ArrayElements:
                Counts_PSD[Time_idx[i],Size_idx[i]]+=1
                PSD_L[Time_idx[i],Size_idx[i]]+= 1/(SVol_L_s[i])  
    
    dNdD_L = PSD_L/ProbeRes # normalise using bin width

    if SaveFile == 1 : 
        SavePath= PathOAP+'dNdD_L_'+SaveName+'.h5'
        file = h5py.File(SavePath, 'w')
        file.create_dataset('dNdD_L', data=dNdD_L)
        file.create_dataset('Counts_PSD', data=Counts_PSD)
        file.create_dataset('TimeMidBins_s', data=TimeMidBins_s)
        file.create_dataset('PSD_SizeMid', data=PSD_SizeMid)
        file.close()






             
#__________________________________________________________________________________

# load particle stats output from oasis. filter using IAT_treshold, MaxAspectRatio, Streak


def LoadFilterOAPh5(PathOAP,filena):
    
    IAT_treshold = 1E-6
    PixelSize = 100 
    MaxAspectRatio = 10
    Streak = 5
    MaxSeconds = 172800 # max number of seconds from mignight.
    
    Data_h5 = h5py.File(PathOAP+filena, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    #TmpStats=np.array(Data_h5['PartGrayStatsWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    Data_h5.close()

    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)
    DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    MaxDiameterBG= np.sqrt(DiameterBGx**2 + DiameterBGy**2)
    StreakFlag =(np.where(np.logical_and(DiameterBGx == PixelSize, DiameterBGy >= Streak*PixelSize), 1, 0))   
    AspectRatio = DiameterBGy / DiameterBGx

    # Nan edge parAspectRatioticles
    shape = np.shape(ParticleStatsWv)
    Edge= np.zeros(shape[0])
    Edge[np.logical_or(ParticleStatsWv[:,16]!=0, ParticleStatsWv[:,17]!=0)] = 1
    IAT= GetIAT_TimeInS_vector(ParticleTimeSeconds)

    for i in range(len(ParticleTimeSeconds)): 
        if Edge[i] == 1 or IAT[i] < IAT_treshold or StreakFlag[i] == 1 or AspectRatio[i] > MaxAspectRatio  or ParticleTimeSeconds[i]>MaxSeconds:
            ParticleTimeSeconds[i] = np.nan
            #MeanXY_BG[i] = np.nan
            #DiameterBGy[i] = np.nan
            #DiameterBGx[i] = np.nan 
            #MaxDiameterBG[i] = np.nan
            #StreakFlag[i] = np.nan
            #AspectRatio[i] = np.nan
    MeanXY_BG = MeanXY_BG[~np.isnan(ParticleTimeSeconds)]
    DiameterBGy = DiameterBGy[~np.isnan(ParticleTimeSeconds)]
    DiameterBGx = DiameterBGx[~np.isnan(ParticleTimeSeconds)]
    MaxDiameterBG = MaxDiameterBG[~np.isnan(ParticleTimeSeconds)]
    ParticleTimeSeconds = ParticleTimeSeconds[~np.isnan(ParticleTimeSeconds)]  
     
    return MeanXY_BG, DiameterBGy, DiameterBGx, MaxDiameterBG, ParticleTimeSeconds


#__________________________________________________________________________________
    
# plot single image from oasis h5 image file

def PlotImage(ImagePath,FileName,ParticleBufferTime,ParticleID,ImageZD,ImageZ, SavePath ):
    #ImagePath = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/Images/'
    #FileName = 'Export_base150313084741.h5'
    PixelSize= 10
    ArrayWidth = 128
    
    Data_h5 = h5py.File(ImagePath + FileName, 'r')
    #ImageData=Data_h5['ImageData']
    ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
    ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
    ImageID =np.array(Data_h5['ImageTimes'][:,2])
       
    ImageSlices[ImageSlices<0] = np.nan
    #Find start position of image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append(0, ImagePosition)
    
    #Search for particle image
    ImageID[ImageTimes != ParticleBufferTime] = np.nan
    i = np.nanargmin((np.absolute(ImageID - ParticleID)))
    #print(ImageID[i])
    #print(ImageTimes[i])
    
    Image = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])-1]) 
    x = np.arange(-PixelSize/2,ArrayWidth*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,(ImageSlices[i]-1)*PixelSize+PixelSize/2, PixelSize)
    #Yscaling = len(y) / len(x) # TAS dimension / Array dimension
    
    #fig=plt.figure(figsize=(10*Yscaling,10))
    if 1==1 :
        fig=plt.figure(figsize=(8,16))
        plt.rcParams.update({'font.size': 12})
        plt.title(str(ImageTimes[i])+'s, ID = ' + str(ImageID[i])+', Z = '+str(ImageZ)+'mm, Zd = '+str(ImageZD))
        plt.pcolormesh(y,x,Image, cmap='Greys_r')
    #    plt2.set_clim(vmin=0, vmax=1.5)
        plt.xlabel('x, μm')
        plt.ylabel('y, μm')
    #    cbar=plt.colorbar(orientation='vertical')
        plt.ylim([0,1280])
        plt.xlim([0,1280])
        plt.savefig(SavePath+str(ImageTimes[i])+'s_ID_' + str(ImageID[i])+'.png',dpi=200)
        plt.close(fig)
    Data_h5.close()
    
    #return Image

#__________________________________________________________________________________
   
# Load single OAP image 

def LoadImage(ImagePath,FileName,ParticleBufferTime,ParticleID, PixelSize, ArrayWidth):
    #ImagePath = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/Images/'
    #FileName = 'Export_base150313084741.h5'
    #PixelSize= 10
    #ArrayWidth = 128
    
    Data_h5 = h5py.File(ImagePath + FileName, 'r')
    ImageData=Data_h5['ImageData']
    ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
    ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
    ImageID =np.array(Data_h5['ImageTimes'][:,2])
    
    #Define new ImageSearchID that is unique to each particle and monotonic
    #ImageID= (np.where(np.logical_or(0<ImageID , ImageID > 65535) , np.nan ,ImageID))
    #ImageSearchID = np.round(ImageTimes, decimals=0) + (ImageID / 65536)
    
    ImageSlices[ImageSlices<0] = np.nan
    #Find start position of image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append(0, ImagePosition)
    
    #Search for particle image
    ImageID[ImageTimes != ParticleBufferTime] = np.nan
    i = np.nanargmin((np.absolute(ImageID - ParticleID)))
    #print('Found image')
    #print(ImageID[i])
    #print(ImageTimes[i])
    
    Image = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])-1]) 
    #Image =np.array(ImageData[:,36:102]) 
    x = np.arange(-PixelSize/2,ArrayWidth*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,(ImageSlices[i]-1)*PixelSize+PixelSize/2, PixelSize)
    Yscaling = len(y) / len(x) # TAS dimension / Array dimension
       
    Data_h5.close()
    return Image, x, y
    

#__________________________________________________________________________________
# normalise a PSD by a pair of moments into a "univerisal" form. See field et al 2007 for details
# "Snow Size Distribution Parameterization for Midlatitude and Tropical Ice Clouds"
# This overwites the input data
  
def Field07_PSD_Normalise(PSD_xMid, PSD_xWidth, PSD_y):
    
    # PSD_xMid		// x bin Midpoints		(um)
    # PSD_xWidth	// x bin Widths		(um)
    # PSD_y			// Number densitys	(#/m3 /um)
    # NormX		// flag to determine if X data should be normalised. 1=normalise. 0=no normalisation
	
    #Calculate the required Moments
    M2 = Field07_PSD_CalcMoment(2, PSD_xMid, PSD_xWidth, PSD_y)
    M3 = Field07_PSD_CalcMoment(3, PSD_xMid, PSD_xWidth, PSD_y)
	
    # Use the moments to normalise the PSD
    PSD_y_norm = PSD_y * ( (M3**3) / (M2**4) )
	
    #Use the moments to normlise the bins 
    PSD_xMid_norm = PSD_xMid * (M2/M3)
        
    return PSD_y_norm, PSD_xMid_norm

#__________________________________________________________________________________

# Calculate Moment of a PSD
def Field07_PSD_CalcMoment(Mnt, PSD_xMid, PSD_xWidth, PSD_y) : 
	# Mnt			// moment of interest
	# PSD_xMid		// x bin Midpoints		(um)
	# PSD_xWidth	// x bin Widths	(um)
	# PSD_y			// Number densities	(#/Volume /um)		Volume is /L, /m3 etc and is preserved into output
	
    PSD_xMid = np.array(PSD_xMid)
    PSD_xWidth = np.array(PSD_xWidth)
    PSD_y = np.array(PSD_y)
    
    
	#var to store Result
    Moment_Calc = 0
	
	#loop through the bins to get the moment value
    
    #for i in range(len(PSD_y)): 
    #    Moment_Calc +=( PSD_y[i]*PSD_xWidth[i]*(PSD_xMid[i]**Mnt) )

    Moment_Calc = np.nansum( PSD_y*PSD_xWidth*(PSD_xMid**Mnt) )
	
	#return the result
    return Moment_Calc

#__________________________________________________________________________________


# Code to generate a lognormal PSD.
# for my example i used the following: LogNormalTest(0.7, ln(100), 100)
#def Field07_LogNormalTest(sigma, mu, Num) : 
#	variable sigma, mu, Num
#	
#	make/o/n=(128) LogNorm_x=(p+1)*10
#	make/o/n=(128) LogNorm_y= (ln(LogNorm_x[p])-mu)^2
#	LogNorm_y/=(-2*sigma*sigma)
#	LogNorm_y = exp(LogNorm_y[p])
#	LogNorm_y*=(Num/(logNorm_x[p]*sigma*sqrt(2*pi)))
#	
#  
    
#__________________________________________________________________________________

def Field07_PSD(d_norm): 
    Field_Trop1=152*np.exp(-12.4*d_norm)						# exponential component
    Field_Trop2=3.28*(d_norm**(-0.78))*np.exp(-1.94*d_norm)	# gamma component

    Field_Trop=Field_Trop1+Field_Trop2

    Field_Mid1=141*np.exp(-16.8*d_norm)						# exponential component
    Field_Mid2=102*(d_norm**(2.07))*np.exp(-4.82*d_norm)		# gamma component

    Field_Mid=Field_Mid1+Field_Mid2

    NewPSD=102*(d_norm**(2.07))*np.exp(-4.82*d_norm)		#pure gamma

    return Field_Trop, Field_Mid, NewPSD

#__________________________________________________________________________________

    
# plot image from each channel side by side single image from oasis h5 image file

def PlotImages_2Channels(ImagePath,FileName,ParticleBufferTime_Ch0,ParticleID_Ch0,ParticleBufferTime_Ch1,ParticleID_Ch1, SavePath ):
    #ImagePath = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/Images/'
    #FileName = 'Export_base150313084741.h5'
    PixelSize= 10
    ArrayWidth = 128
    Data_h5 = h5py.File(ImagePath + FileName, 'r')
    ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
    ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
    ImageID_Ch0 =np.array(Data_h5['ImageTimes'][:,2])
    ImageID_Ch1 =np.array(Data_h5['ImageTimes'][:,2])
    ImageSlices[ImageSlices<0] = np.nan
    #Find start position of image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append(0, ImagePosition)
    
    fig=plt.figure(figsize=(12,5))
    plt.rcParams.update({'font.size': 12})
    
    # Channel 0 
    #Search for particle image
    plt.subplot(1, 2, 1)
    ImageID_Ch0[ImageTimes != ParticleBufferTime_Ch0] = np.nan
    i = np.nanargmin((np.absolute(ImageID_Ch0 - ParticleID_Ch0)))
    Image = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])-1]) 
    x = np.arange(-PixelSize/2,ArrayWidth*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,(ImageSlices[i]-1)*PixelSize+PixelSize/2, PixelSize)
    plt.title(str(ImageTimes[i])+'s, ID = ' + str(ImageID_Ch0[i]))
    plt.pcolormesh(y,x,Image, cmap='Greys_r')
    plt.xlabel('y, μm')
    plt.ylabel('x, μm')
    plt.ylim([0,1280])
    plt.xlim([0,1280])
        
    # Channel 1
    #Search for particle image
    plt.subplot(1, 2, 2)
    ImageID_Ch1[ImageTimes != ParticleBufferTime_Ch1] = np.nan
    i = np.nanargmin((np.absolute(ImageID_Ch1 - ParticleID_Ch1)))
    Image = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])-1]) 
    x = np.arange(-PixelSize/2,ArrayWidth*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,(ImageSlices[i]-1)*PixelSize+PixelSize/2, PixelSize)
    plt.title(str(ImageTimes[i])+'s, ID = ' + str(ImageID_Ch1[i]))
    plt.pcolormesh(y,x,Image, cmap='Greys_r')
    plt.xlabel('y, μm')
    plt.ylabel('x, μm')
    plt.ylim([0,1280])
    plt.xlim([0,1280])
    
    
    plt.savefig(SavePath+str(ImageTimes[i])+'s_ID_' + str(ImageID_Ch0[i])+'.png',dpi=200)
    plt.close(fig)
    Data_h5.close()

#__________________________________________________________________________________

#__________________________________________________________________________________
    
# plot single greyscale image from oasis h5 image file

def PlotImage_GS(ImagePath,FileName,ParticleBufferTime,ParticleID,ImageZD,ImageZ, SavePath ):
    #ImagePath = 'C:/Users/Admin TEMP/Documents/CIRCCREX_B895/2ds/Images/'
    #FileName = 'Export_base150313084741.h5'
    PixelSize= 15
    ArrayWidth = 64
    
    Data_h5 = h5py.File(ImagePath + FileName, 'r')
    ImageData=Data_h5['ImageData']
    ImageTimes=np.array(Data_h5['ImageTimes'][:,0])
    ImageSlices  =np.array(Data_h5['ImageTimes'][:,1])
    ImageID =np.array(Data_h5['ImageTimes'][:,2])
       
    ImageSlices[ImageSlices<0] = np.nan
    #Find start position of image within ImageData
    ImagePosition = np.cumsum(ImageSlices, axis = 0)
    ImagePosition = np.append(0, ImagePosition)
    
    #Search for particle image
    ImageID[ImageTimes != ParticleBufferTime] = np.nan
    i = np.nanargmin((np.absolute(ImageID - ParticleID)))

    
    Image = np.array(Data_h5['ImageData'][:,int(ImagePosition[i]):int(ImagePosition[i+1])-1]) 
    x = np.arange(-PixelSize/2,ArrayWidth*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,(ImageSlices[i]-1)*PixelSize+PixelSize/2, PixelSize)
    #Yscaling = len(y) / len(x) # TAS dimension / Array dimension
    
    Image[Image == 255] = 0 
    
    
    #fig=plt.figure(figsize=(10*Yscaling,10))
    if 1==1 :
        fig=plt.figure(figsize=(10,10))
        plt.rcParams.update({'font.size': 12})
        plt.title(str(ImageTimes[i])+'s, ID = ' + str(ImageID[i])+', Z = '+str(ImageZ)+'mm, Zd = '+str(ImageZD))
        plt.pcolormesh(y,x,Image, cmap='gist_stern_r', vmin=0, vmax=3)
        #plt2.set_clim(vmin=0, vmax=3)
        plt.xlabel('x, μm')
        plt.ylabel('y, μm')
    #    cbar=plt.colorbar(orientation='vertical')
        plt.ylim([0,960])
        plt.xlim([0,1920])
        plt.savefig(SavePath+str(ImageTimes[i])+'s_ID_' + str(ImageID[i])+'.png',dpi=200)
        plt.close(fig)
    Data_h5.close()