# -*- coding: utf-8 -*-
"""
Created on Thu May 24 10:45:55 2018

@author: Admin TEMP
"""

import datetime
import numpy as np
from netCDF4 import Dataset
import math
import bisect
import h5py
import matplotlib.pyplot as plt 
import scipy
from scipy.optimize import fsolve
import scipy.optimize as opt
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
from MyFunctions import LoadOAP_nc, loadmat, TimeSeriesPlot, Matlab2PythonTimeArray, Matlab2PythonTime, Time2Secs,Load_HaloHolo_PSD
from MyFunctions import BinHaloHolo_v2,CalculateEffectiveDiameter, CalculateVolumeMeanDiameter,datetime_2_s_sincemidnight
#from Mass_Dimension import Minimise_Fb, AvgPSDcompositePSD,Find_a_b, CalculateReflectivityIce

from GrayscaleFunctions import GetIAT_TimeInS_vector,LoadGrayScaleStats,GetIAT_vector, ParticlesStats2PSD, GrayscaleRatios, DvsSvol

from PICASSO_functions import Avg_dndDp_CountStats

from DropGun import DimensionlessDistance, LoadModelGrayScale, ModelGS_Stats, LoadGrayScaleStats_Stage, Dof_lawson
import pandas as pd
from dateutil import parser
import matplotlib.dates as mdates

# 


Path2DS = 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/rawData/C078/2DS/Output/'
filena = 'base180207155526.h5'

ColocationThreshold = 1E-8
IAT_treshold = 1E-6
SaveFile =1
PixelSize = 10 
    
Data_h5 = h5py.File(Path2DS+filena, 'r')              
ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
Data_h5.close()

ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
Channel = ParticleTimesWv[:,4] 
MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)


    # Nan edge particles
#Constant K_OASIS_PSTAT_COL_N_EDGE0 = 13
#Constant K_OASIS_PSTAT_COL_N_EDGE1 = 14
shape = np.shape(ParticleStatsWv)
Edge= np.zeros(shape[0])
for i in range(shape[0]):
    if (ParticleStatsWv[i,3]==0 or ParticleStatsWv[i,4]==127):
        Edge[i] = 1


MeanXY_BG_CH0 = MeanXY_BG[Channel == 0]
Seconds_CH0 = ParticleTimeSeconds[Channel == 0]
Edge_CH0 = Edge[Channel == 0]

MeanXY_BG_CH1 = MeanXY_BG[Channel == 1]
Seconds_CH1 = ParticleTimeSeconds[Channel == 1]
Edge_CH1 = Edge[Channel == 1]

#remove shattering. Need to do IAT filtering for each channel
IAT_CH0 = GetIAT_TimeInS_vector(Seconds_CH0)  
Seconds_CH0[IAT_CH0<IAT_treshold] = np.nan
MeanXY_BG_CH0[IAT_CH0<IAT_treshold] = np.nan

IAT_CH1 = GetIAT_TimeInS_vector(Seconds_CH1)  
Seconds_CH1[IAT_CH1<IAT_treshold] = np.nan
MeanXY_BG_CH1[IAT_CH1<IAT_treshold] = np.nan

#remove edge particles
Seconds_CH0[Edge_CH0 == 1] = np.nan
MeanXY_BG_CH0[Edge_CH0 == 1] = np.nan
Seconds_CH1[Edge_CH1 == 1] = np.nan
MeanXY_BG_CH1[Edge_CH1 == 1] = np.nan

#remove nans
Seconds_CH1 *= MeanXY_BG_CH1/ MeanXY_BG_CH1
Seconds_CH0 *= MeanXY_BG_CH0/ MeanXY_BG_CH0


Seconds_CH1 = Seconds_CH1[~np.isnan(Seconds_CH1)]
MeanXY_BG_CH1 = MeanXY_BG_CH1[~np.isnan(MeanXY_BG_CH1)]
Seconds_CH0 = Seconds_CH0[~np.isnan(Seconds_CH0)]
MeanXY_BG_CH0 = MeanXY_BG_CH0[~np.isnan(MeanXY_BG_CH0)]

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

ChIDX_low = ChIDX_high - 1
ChIDX_low[ChIDX_low<0] = 0
ChTimeDelta_low = np.array([np.absolute(Seconds_CH1[x] - Seconds_CH0[int(ChIDX_low[x])]) for x in range(len(ChIDX_low))])


ChIDX = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChIDX_low, ChIDX_high))
ChTimeDelta = (np.where(ChTimeDelta_low < ChTimeDelta_high,ChTimeDelta_low, ChTimeDelta_high))


#    if ChIDX_low[i] > 0 : 
#        ChIDX_low[i] = ChIDX_high[i]-1
#        ChTimeDelta_low[i]=np.absolute(Seconds_CH1[i] - Seconds_CH0[int(ChIDX_low[i])])
#    else: 
#        ChIDX_low[i] = ChIDX_high[i]
#        ChTimeDelta_low[i]=ChTimeDelta_high[i]
#    
#    if ChTimeDelta_low[i] < ChTimeDelta_high[i]:
#        ChIDX[i] = ChIDX_low[i]
#        ChTimeDelta[i]=ChTimeDelta_low[i]
#    else : 
#        ChIDX[i] = ChIDX_high[i]
#        ChTimeDelta[i]=ChTimeDelta_high[i]

ColocationIDX = (ChIDX[ChTimeDelta < ColocationThreshold]).astype(int) # Indexes are for channel 0 (same length as channel 1)
ColocationDelta = ChTimeDelta[ChTimeDelta < ColocationThreshold] 
    #ColocationParticleTime_CH1 = ParticleTime_CH1[ChTimeDelta < ColocationThreshold]
ColocationMeanXY_BG_CH1 = MeanXY_BG_CH1[ChTimeDelta < ColocationThreshold]
ColocationSecondsCH1 = Seconds_CH1[ChTimeDelta < ColocationThreshold]

ColocationSecondsCH0 = Seconds_CH0[ColocationIDX]
ColocationMeanXY_BG_CH0 = MeanXY_BG_CH0[ColocationIDX]


# Counts vs size for all_accept
ExpDate= datetime.datetime(2018, 2, 7, 0, 0, 0)
PltTitle='C078_1602'
StartAvg=datetime.datetime(2018, 2, 7, 16, 2, 0)
EndAvg=datetime.datetime(2018, 2, 7, 16, 10, 0) 
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

Size_idx = np.digitize(ColocationMeanXY_BG_CH0, SizeBins)
     
# Calculate concentration of individual particle   
ArrayWidth_mm = ((ArrayElements*ProbeRes) - (ColocationMeanXY_BG_CH0))*0.001 # mm
SVol_mm3_s = TAS_mm * (ArrayWidth_mm**2)
SVol_L_s = SVol_mm3_s / 1E6
SVol_L= SVol_L_s * Ns

    # Filter and bin data
Counts_PSD_colocation = np.zeros(len(SizeBins))
PSD_L_colocation = np.zeros(len(SizeBins)) 
    
# for colocated data
for i in range(len(ColocationSecondsCH0)):      
    if (ColocationSecondsCH0[i] > Startime) & (ColocationSecondsCH0[i] < Endtime) : 
        Counts_PSD_colocation[Size_idx[i]]+=1
        PSD_L_colocation[Size_idx[i]]+= 1/(SVol_L[i])
    
#PSD_L_colocation /= Ns # convert to #/L
dNdD_L_colocation = PSD_L_colocation/ProbeRes # normalise using bind width
dNdD_L_colocation[0:2] =np.nan # nan first two bins
dNdD_L_colocation_error = dNdD_L_colocation / (2*np.sqrt(Counts_PSD_colocation)) # /2 since plot plus minus

# total counts
Size_idx = np.digitize(MeanXY_BG_CH0, SizeBins)
Counts_PSD_total = np.zeros(len(SizeBins))
PSD_L_total = np.zeros(len(SizeBins))

ArrayWidth_mm = ((ArrayElements*ProbeRes) - (MeanXY_BG_CH0))*0.001 # mm
SVol_mm3_s = ArrayWidth_mm * TAS_mm * np.minimum(Dof_lawson(MeanXY_BG_CH0,0.685,8),63)
SVol_L_s = SVol_mm3_s / 1E6
SVol_L= SVol_L_s * Ns


for i in range(len(Seconds_CH0)):
    if (Seconds_CH0[i]==Seconds_CH0[i]):     
        if (Seconds_CH0[i] > Startime) & (Seconds_CH0[i] < Endtime) : 
            if (Size_idx[i]<128):
                Counts_PSD_total[Size_idx[i]]+=1
                PSD_L_total[Size_idx[i]]+= 1/(SVol_L[i])
      
dNdD_L_total = PSD_L_total/ProbeRes # normalise using bind width

          
if 1==1 : 
    fig=plt.figure(figsize=(7,7)) 
    plt.rcParams.update({'font.size': 12})
#    plt.subplot(3, 1, 1)
#    plt.plot(ColocationMeanXY_BG_CH1, ColocationMeanXY_BG_CH0,'o',markersize=2,label='Both channels')
#    plt.plot(ColocationMeanXY_BG_CH1,ColocationMeanXY_BG_CH1,label='1 to 1')
#    plt.legend()
#    plt.xlabel('Diameter (Channel 1), μm')
#    plt.ylabel('Diameter (Channel 0), μm')
#    plt.xscale('log')
#    plt.yscale('log')
#    
#    plt.subplot(3, 1, 2)
#    plt.plot(SizeBins,Counts_PSD_total, label='All')
#    plt.plot(SizeBins,Counts_PSD_colocation, label='Both channels')
#    plt.legend()
#    plt.xlabel('Diameter, μm')
#    plt.ylabel('Counts')
#    plt.xscale('log')
#    plt.yscale('log')
    
    plt.subplot(1, 1, 1)
    plt.errorbar(PSD_SizeMid, dNdD_L_colocation, yerr =dNdD_L_colocation_error, marker='o',linestyle = 'None', label='2DS (Both channels)')
    plt.plot(PSD_SizeMid, dNdD_L_total, label='2DS (All)')
    
    HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C078/CIP15/'
    HaloName='C078_forSeb_v1.h5'
    BinMids, HaloHolo_dNdDp_L_avg, CountUncertainty=BinHaloHolo_v2(0.01, StartAvg, EndAvg, HaloName,HaloPath)
    plt.plot(BinMids,HaloHolo_dNdDp_L_avg,'o',label='HALOHolo',color= 'tab:pink')
    plt.legend()
    plt.xlabel('Diameter, μm')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    plt.yscale('log')
    plt.xscale('log')

    Path2DS='C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/'
        #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
#    plt.savefig(Path2DS+'Colocation2DS'+PltTitle,dpi=200)
#    plt.close(fig)

