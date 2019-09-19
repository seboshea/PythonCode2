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
    for filena in os.listdir(Path2DS):
        if filena.endswith(".h5") and filena.startswith('base'):
            print(filena)
            FindParticlesOnBothChannels(Path2DS,filena)
            HistMIDPosition2DS(Path2DS,filena)

#__________________________________________________________________________________
#
            
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
    
    if 1 == 1:
        fig=plt.figure(figsize=(7,7)) 
        plt.rcParams.update({'font.size': 12})
        plt.plot(HistBinsMid,MIDx_Ch0_hist, label= 'Ch = 0')
        plt.plot(HistBinsMid,MIDx_Ch1_hist, label= 'Ch = 1' )
        plt.xlabel('Array elements')
        plt.ylabel('Counts')
        plt.legend()
        plt.savefig(Path2DS+filena[:-3]+'_MIDx_edgeRej.png',dpi=200)
        plt.close(fig)

#__________________________________________________________________________________
# Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/'
# filena = 'base180207155526.h5'

def FindParticlesOnBothChannels(Path2DS,filena):
    ColocationThreshold = 5E-7
    IAT_treshold = 1E-5
    SaveFile =1
    PixelSize = 10 
    MaxAspectRatio = 10
    Streak = 5
    
    Data_h5 = h5py.File(Path2DS+filena, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    #TmpStats=np.array(Data_h5['PartGrayStatsWv'])
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
    MIDx = (ParticleStatsWv[:,10] + ParticleStatsWv[:,9])/2

    # Nan edge parAspectRatioticles
    shape = np.shape(ParticleStatsWv)
    Edge= np.zeros(shape[0])
    for i in range(shape[0]):
        if (ParticleStatsWv[i,3]==0 or ParticleStatsWv[i,4]==127):
            Edge[i] = 1

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
    MIDx_Ch0 = MIDx[Channel == 0]
    MIDx_Ch1 = MIDx[Channel == 1]

    #Get IAT times for each channel seprately
    IAT_CH0 = GetIAT_TimeInS_vector(Seconds_CH0)    
    Seconds_CH0[IAT_CH0<IAT_treshold] = np.nan
    MeanXY_BG_CH0[IAT_CH0<IAT_treshold] = np.nan
    IAT_CH1 = GetIAT_TimeInS_vector(Seconds_CH1)
    #print('No filtering ' + str(len(IAT_CH1)))
    Seconds_CH1[IAT_CH1<IAT_treshold] = np.nan
    MeanXY_BG_CH1[IAT_CH1<IAT_treshold] = np.nan
    MaxDiameterBG_Ch0[IAT_CH0<IAT_treshold] = np.nan
    MaxDiameterBG_Ch1[IAT_CH1<IAT_treshold] = np.nan
    MIDx_Ch0[IAT_CH0<IAT_treshold] = np.nan
    MIDx_Ch1[IAT_CH1<IAT_treshold] = np.nan
    
    # Remove high aspect ratio particles, likely streaking   
    Seconds_CH0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    MeanXY_BG_CH0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    Seconds_CH1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan
    MeanXY_BG_CH1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan
    IAT_CH0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    IAT_CH1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan
    MaxDiameterBG_Ch0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    MaxDiameterBG_Ch1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan     
    MIDx_Ch0[AspectRatio_Ch0 >= MaxAspectRatio] = np.nan
    MIDx_Ch1[AspectRatio_Ch1 >= MaxAspectRatio] = np.nan
    
    
    Seconds_CH0[Streak_Ch0 == 1] = np.nan
    MeanXY_BG_CH0[Streak_Ch0 == 1] = np.nan
    Seconds_CH1[Streak_Ch1 == 1] = np.nan
    MeanXY_BG_CH1[Streak_Ch1 == 1] = np.nan
    IAT_CH0[Streak_Ch0 == 1] = np.nan
    IAT_CH1[Streak_Ch1 == 1] = np.nan
    MaxDiameterBG_Ch0[Streak_Ch0 == 1] = np.nan
    MaxDiameterBG_Ch1[Streak_Ch1 == 1] = np.nan
    MIDx_Ch0[Streak_Ch0 == 1] = np.nan
    MIDx_Ch1[Streak_Ch1 == 1] = np.nan
    
    
#    print(np.nanmean(Streak_Ch0))
#    print(np.nanmean(Streak_Ch1))
#    print(np.nanmedian(AspectRatio_Ch0))
#    print(np.nanmedian(AspectRatio_Ch1))
    
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

    #remove edge particles
    Seconds_CH0[Edge_CH0 == 1] = np.nan
    MeanXY_BG_CH0[Edge_CH0 == 1] = np.nan
    Seconds_CH1[Edge_CH1 == 1] = np.nan
    MeanXY_BG_CH1[Edge_CH1 == 1] = np.nan
    IAT_CH0[Edge_CH0 == 1] = np.nan
    IAT_CH1[Edge_CH1 == 1] = np.nan
    MaxDiameterBG_Ch0[Edge_CH0 == 1] = np.nan
    MaxDiameterBG_Ch1[Edge_CH1 == 1] = np.nan
    MIDx_Ch0[Edge_CH0 == 1] = np.nan
    MIDx_Ch1[Edge_CH1 == 1] = np.nan
    
    #remove nans
    MeanXY_BG_CH0 *= MaxDiameterBG_Ch0/MaxDiameterBG_Ch0
    MeanXY_BG_CH1 *= MaxDiameterBG_Ch1/MaxDiameterBG_Ch1
    Seconds_CH1 *= MeanXY_BG_CH1/ MeanXY_BG_CH1
    Seconds_CH0 *= MeanXY_BG_CH0/ MeanXY_BG_CH0    
    IAT_CH0*= MeanXY_BG_CH0/ MeanXY_BG_CH0
    IAT_CH1*= MeanXY_BG_CH1/ MeanXY_BG_CH1
    MIDx_Ch0*= MeanXY_BG_CH0/ MeanXY_BG_CH0
    MIDx_Ch1*= MeanXY_BG_CH1/ MeanXY_BG_CH1
   
    Seconds_CH1 = Seconds_CH1[~np.isnan(Seconds_CH1)]
    MeanXY_BG_CH1 = MeanXY_BG_CH1[~np.isnan(MeanXY_BG_CH1)]
    Seconds_CH0 = Seconds_CH0[~np.isnan(Seconds_CH0)]
    MeanXY_BG_CH0 = MeanXY_BG_CH0[~np.isnan(MeanXY_BG_CH0)]
    IAT_CH0 = IAT_CH0[~np.isnan(IAT_CH0)]
    IAT_CH1 = IAT_CH1[~np.isnan(IAT_CH1)]
    MaxDiameterBG_Ch0= MaxDiameterBG_Ch0[~np.isnan(MaxDiameterBG_Ch0)]
    MaxDiameterBG_Ch1= MaxDiameterBG_Ch1[~np.isnan(MaxDiameterBG_Ch1)]
    MIDx_Ch0= MIDx_Ch0[~np.isnan(MIDx_Ch0)]
    MIDx_Ch1= MIDx_Ch1[~np.isnan(MIDx_Ch1)]
    
#    print(len(MeanXY_BG_CH0))
#    print(len(MaxDiameterBG_Ch0))
#    print(len(MeanXY_BG_CH1))
#    print(len(MaxDiameterBG_Ch1))
#    print('With filtering ' + str(len(IAT_CH1)))
    
    # Look for colocated particles. 
    ChTimeDelta_high= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX_high= np.zeros(len(Seconds_CH1))*np.nan
    ChTimeDelta_low= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX_low= np.zeros(len(Seconds_CH1))*np.nan
    ChTimeDelta= np.zeros(len(Seconds_CH1))*np.nan
    ChIDX= np.zeros(len(Seconds_CH1))*np.nan
    

    for i in range(len(Seconds_CH1)) : 
        ChIDX_high[i] = np.searchsorted(Seconds_CH0, Seconds_CH1[i], side="left")
        if (ChIDX_high[i]<len(Seconds_CH0)):
            ChTimeDelta_high[i]=np.absolute(Seconds_CH1[i] - Seconds_CH0[int(ChIDX_high[i])])
    #need to check whether i or i-1 gives lowest time delta
    ChIDX_low = ChIDX_high - 1
    ChIDX_low[ChIDX_low<0] = 0
    ChTimeDelta_low = np.array([np.absolute(Seconds_CH1[x] - Seconds_CH0[int(ChIDX_low[x])]) for x in range(len(ChIDX_low))])
    ChIDX = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChIDX_low, ChIDX_high))
    ChTimeDelta = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChTimeDelta_low, ChTimeDelta_high))


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

    ColocationIDX = (ChIDX[ChTimeDelta <= ColocationThreshold]).astype(int) # Indexes are for channel 0 (same length as channel 1)
    ColocationDelta = ChTimeDelta[ChTimeDelta < ColocationThreshold] 
    #ColocationParticleTime_CH1 = ParticleTime_CH1[ChTimeDelta < ColocationThreshold]
    ColocationMeanXY_BG_CH1 = MeanXY_BG_CH1[ChTimeDelta < ColocationThreshold]
    ColocationMaxDiameterBG_Ch1 = MaxDiameterBG_Ch1[ChTimeDelta < ColocationThreshold]
    ColocationSecondsCH1 = Seconds_CH1[ChTimeDelta < ColocationThreshold]
    ColocationMIDx_Ch1= MIDx_Ch1[ChTimeDelta < ColocationThreshold]
    
    ColocationSecondsCH0 = Seconds_CH0[ColocationIDX]
    ColocationMeanXY_BG_CH0 = MeanXY_BG_CH0[ColocationIDX]
    ColocationMaxDiameterBG_Ch0 = MaxDiameterBG_Ch0[ColocationIDX]
    ColocationMIDx_Ch0= MIDx_Ch0[ColocationIDX]

    # select inter arrival time for only colocated particles 
    IAT_CH0_colocation = IAT_CH0[ColocationIDX]
    IAT_CH1_colocation = IAT_CH1[ChTimeDelta < ColocationThreshold]
    IAT_colocation = np.minimum(IAT_CH0_colocation,IAT_CH1_colocation)
    
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
        file.create_dataset('ColocationMIDx_Ch0', data=ColocationMIDx_Ch0)
        file.create_dataset('ColocationMIDx_Ch1', data=ColocationMIDx_Ch1)
        
        file.close()
        #Create PSDs using all data
        #PSD_Colocate_1hz(Seconds_CH0, MeanXY_BG_CH0, Seconds_CH1, MeanXY_BG_CH1,SaveFile,Path2DS,filena,0)
        PSD_Colocate_1hz(Seconds_CH0, MaxDiameterBG_Ch0, Seconds_CH1, MaxDiameterBG_Ch1,SaveFile,Path2DS,filena,0)
        #Create PSDs using just colocated data
        #PSD_Colocate_1hz(ColocationSecondsCH0, ColocationMeanXY_BG_CH0, ColocationSecondsCH1, ColocationMeanXY_BG_CH1,SaveFile,Path2DS,filena,1)
        PSD_Colocate_1hz(ColocationSecondsCH0, ColocationMaxDiameterBG_Ch0, ColocationSecondsCH1, ColocationMaxDiameterBG_Ch1,SaveFile,Path2DS,filena,1)
 
        
        #ColocatePSD_1hz(Path2DS, 'Colocate_'+filena)
    
    


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

def Colocated_sVol_array(Diameter_CH0, Diameter_CH1,ArrayElements,ProbeRes,TAS ):
        # Calculate concentration of individual particle   
    ArrayWidthCH0_mm = ((ArrayElements*ProbeRes) - (Diameter_CH0))*0.001 # mm
    ArrayWidthCH1_mm = ((ArrayElements*ProbeRes) - (Diameter_CH1))*0.001 # mm
    TAS*=1000 # mm/s
    SVol_mm3_s = TAS * ArrayWidthCH0_mm * ArrayWidthCH1_mm
    SVol_L_s = SVol_mm3_s / 1E6
    
    return SVol_L_s

#__________________________________________________________________________________
 

#Create arrays with dNdD_L_CH0 and dNdD_L_CH1. 
#ColocatedFlag =0  Lawson et al sample volume 
#ColocatedFlag =1  Lawson et al sample volume 

def PSD_Colocate_1hz(SecondsCH0, DiameterCH0, SecondsCH1, DiameterCH1,SaveFile,Path2DS,filena,ColocatedFlag):
    
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
        SVolCH0_L_s = Colocated_sVol_array(DiameterCH0, DiameterCH1,ArrayElements,ProbeRes,TAS )
        SVolCH1_L_s = Colocated_sVol_array(DiameterCH0, DiameterCH1,ArrayElements,ProbeRes,TAS )
    else : 
        SVolCH0_L_s = Lawson_sVol_array(DiameterCH0,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS)    
        SVolCH1_L_s = Lawson_sVol_array(DiameterCH1,Lambda,c,ArmSep,ArrayElements,ProbeRes,TAS) 

    # Filter and bin data
    Counts_PSD_CH0 = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L_CH0 = np.zeros((len(TimeBins_s), len(SizeBins)))

    Counts_PSD_CH1 = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L_CH1 = np.zeros((len(TimeBins_s), len(SizeBins)))
        
    for i in range(len(SecondsCH0)):      
        if (SecondsCH0[i] > Startime) & (SecondsCH0[i] < Endtime) : 
            if SizeCH0_idx[i] < ArrayElements:
                Counts_PSD_CH0[TimeCH0_idx[i],SizeCH0_idx[i]]+=1
                PSD_L_CH0[TimeCH0_idx[i],SizeCH0_idx[i]]+= 1/(SVolCH0_L_s[i])
            
    for i in range(len(SecondsCH1)):      
        if (SecondsCH1[i] > Startime) & (SecondsCH1[i] < Endtime) : 
            if SizeCH1_idx[i] < ArrayElements:
                Counts_PSD_CH1[TimeCH1_idx[i],SizeCH1_idx[i]]+=1
                PSD_L_CH1[TimeCH1_idx[i],SizeCH1_idx[i]]+= 1/(SVolCH1_L_s[i])     
    
    dNdD_L_CH0 = PSD_L_CH0/ProbeRes # normalise using bin width
    #dNdD_L_CH0[:,0:2] =np.nan # nan first two bins
    #dNdD_L_CH0_error = dNdD_L_CH0 / (2*np.sqrt(Counts_PSD_CH0)) # /2 since plot plus minus

    dNdD_L_CH1 = PSD_L_CH1/ProbeRes # normalise using bin width
    
    #dNdD_L_CH1_error = dNdD_L_CH1 / (2*np.sqrt(Counts_PSD_CH1)) # /2 since plot plus minus

    if ColocatedFlag == 1:
        dNdD_L_CH0[:,0:2] =np.nan # nan first two bins
        dNdD_L_CH1[:,0:2] =np.nan # nan first two bins

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
    
##Create arrays with dNdD_L_CH0_colocation vs t
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
        
def ColocateB895():
    
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
    
#    fig=plt.figure(figsize=(9,18)) 
#    plt.rcParams.update({'font.size': 12})
    for i in range(len(IDX_Start_colocation)):
        
        
        fig=plt.figure(figsize=(7,7)) 
        StartIdx = IDX_Start_CIP100[i]
        EndIdx = IDX_End_CIP100[i]
        dNdD_L_CIP100_avg = np.nanmean(dNdD_L_CIP100[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CIP100_avg[0] =np.nan # Ditch first bin
        dNdD_L_CIP100_error = dNdD_L_CIP100_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CIP100[StartIdx:EndIdx][:],axis=0)))
        #plt.errorbar(PSD_SizeMid_CIP100, dNdD_L_CIP100_avg,yerr=(dNdD_L_CIP100_error), capsize=2, label= 'CIP100')
        
        
        StartIdx = IDX_Start[i]
        EndIdx = IDX_End[i] 
        dNdD_L_CH0_avg = np.nanmean(dNdD_L_CH0[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_avg = np.nanmean(dNdD_L_CH1[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_error = dNdD_L_CH0_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH0[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_error = dNdD_L_CH1_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH1[StartIdx:EndIdx][:],axis=0)))
        
        StartIdx = IDX_Start_colocation[i]
        EndIdx = IDX_End_colocation[i]
        dNdD_L_CH0_colocation_avg = np.nanmean(dNdD_L_CH0_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH1_colocation_avg = np.nanmean(dNdD_L_CH1_colocation[StartIdx:EndIdx][:],axis=0)
        dNdD_L_CH0_colocation_error = dNdD_L_CH0_colocation_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH0_colocation[StartIdx:EndIdx][:],axis=0)))
        dNdD_L_CH1_colocation_error = dNdD_L_CH1_colocation_avg / (2 * np.sqrt(np.nansum(Counts_PSD_CH1_colocation[StartIdx:EndIdx][:],axis=0)))

        HaloRun = HaloHoloData[str(RunTemperature[i])+'K']
        HaloHoloPSD = np.array(HaloRun[1:].astype(np.float))
        HaloHoloPSD[0:6] = np.nan
        
        plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH0_avg,yerr=(dNdD_L_CH0_error), color = 'tab:gray', label = '2DS all')
        plt.errorbar(PSD_SizeMid_colocation,dNdD_L_CH1_avg,yerr=(dNdD_L_CH1_error), color = 'tab:gray')
        plt.plot(PSD_SizeMid_CIP100, dNdD_L_CIP100_avg,color='k', label= 'CIP100')
        plt.errorbar(PSD_SizeMid,dNdD_L_CH0_colocation_avg,yerr=(dNdD_L_CH0_colocation_error), color = 'r', fmt='o',markersize=5, capsize=2, label = 'CH0 colocation')
        plt.errorbar(PSD_SizeMid,dNdD_L_CH1_colocation_avg,yerr=(dNdD_L_CH1_colocation_error), color = 'c', fmt='o',markersize=5, capsize=2, label = 'CH1 colocation')
        plt.plot(HaloholoMidsize, HaloHoloPSD,color='b', label='HALOHolo')
        
#        OASISRun = OASISData[str(RunTemperature[i])+'K']
#        OASISPSD = np.array(OASISRun[1:].astype(np.float))
#        plt.plot(OASISMidsize, OASISPSD, label='OASIS 2DS')             
        
        plt.title(str(RunTemperature[i])+' K')
        plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
        plt.xlabel('Diameter, μm')      
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim([10,1000])
        plt.ylim([1E-4,1E1])
        plt.legend()
        plt.savefig(Path2DS+'Psd'+str(RunTemperature[i])+'K'+'.png',dpi=200)
        plt.close(fig)
        
        
   
    
#        plt.subplot(4, 2, i+1)
#        plt.plot(PSD_SizeMid,dNdD_L_CH0_colocation)
#        plt.plot(PSD_SizeMid,dNdD_L_CH1_colocation)
#        
#        HaloRun = HaloHoloData[str(RunTemperature[i])+'K']
#        HaloHoloPSD = HaloRun[1:].astype(np.float)
#        HaloHoloPSD[0:7] = np.nan
#        plt.plot(HaloholoMidsize, HaloHoloPSD)
#        
#        plt.title(str(RunTemperature[i])+' K')
#        if ((i+1) % 2 != 0) :  # odd panel
#            plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
#        if (i+ 1 == 7 ) or (i + 1 == 6) : 
#            plt.xlabel('Diameter, μm')      
#        plt.yscale('log')
#        plt.xscale('log')
#        plt.xlim([10,1000])
#        plt.ylim([1E-3,1])
    

    

#__________________________________________________________________________________

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
 

#Create arrays with dNdD_L_CH0 and dNdD_L_CH1. 
#ColocatedFlag =0  Lawson et al sample volume 
#ColocatedFlag =1  Lawson et al sample volume 

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
 