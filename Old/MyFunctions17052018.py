# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:54:16 2018

@author: Admin TEMP
"""
import datetime
import numpy as np
from netCDF4 import Dataset
import math
import bisect
import h5py




#__________________________________________________________________________________________________
#Bin DataY based on DataX value


def BinLinear(DataY,DataX,StartBin,EndBin, NBins):
    Output=np.zeros((len(DataX),NBins))*np.nan
    bins = np.linspace(StartBin,EndBin, NBins)
    digitized = np.digitize(DataX, bins)
    
    for i in range(len(DataY)-1):
        if DataX[i]>=StartBin and DataX[i]<EndBin :
            Output[i,digitized[i]]=DataY[i]

    return Output, bins

#__________________________________________________________________________________________________



def FindFirstGreaterThan(Matrix,Threshold):
    return next(x[0] for x in enumerate(Matrix) if x[1] > Threshold)

#__________________________________________________________________________________________________

def remove_prefix(text, prefix):
    return text[text.startswith(prefix) and len(prefix):]
#__________________________________________________________________________________________________

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#__________________________________________________________________________________________________



def DateTime2IgorTime(DateTimeArray):
    dt_base = datetime.datetime(1904, 1, 1, 0, 0, 0)
    #IgorTime = (DateTimeArray - dt_base).total_seconds()
    
    IgorTime= [(DateTimeArray[x]-dt_base).total_seconds() for x in range(len(DateTimeArray))]   
    return IgorTime

#__________________________________________________________________________________________________

def Igor2DateTime(IgorTime):
    dt_base = datetime.datetime(1904, 1, 1, 0, 0, 0)
    Dt= [dt_base + datetime.timedelta(seconds=IgorTime[x]) for x in range(len(IgorTime))]
    return Dt

#__________________________________________________________________________________________________

def loadTAS(FullPath,FAAMCoreName):
    FAAMCore = Dataset(FullPath+FAAMCoreName)
    Time=np.array(FAAMCore['Time'][:])
    TAS_RVSM=np.array(FAAMCore['TAS_RVSM'][:])
    TAS_RVSM_FLAG=np.array(FAAMCore['TAS_RVSM_FLAG'][:])
    dt_base = datetime.datetime(1904, 1, 1, 0, 0, 0)
    FlightDate= datetime.datetime(int(FAAMCoreName[10:14]), int(FAAMCoreName[14:16]), int(FAAMCoreName[16:18]), 0, 0, 0)
    print(FlightDate)
    TotalSeconds = (FlightDate - dt_base).total_seconds()
    DateTime=TotalSeconds+Time
    TAS_RVSM[TAS_RVSM_FLAG!=0]=np.nan
    
    return TAS_RVSM, DateTime

#__________________________________________________________________________________________________
#

#FullPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
#FAAMCoreName='core_faam_20180214_v004_r0_c082_1hz.nc'

def loadFAAMCore(FullPath,FAAMCoreName):
    FAAMCore = Dataset(FullPath+FAAMCoreName)
    Time=np.array(FAAMCore['Time'][:])
    FlightDate= datetime.datetime(int(FAAMCoreName[10:14]), int(FAAMCoreName[14:16]), int(FAAMCoreName[16:18]), 0, 0, 0)
    Time_Core= [FlightDate + datetime.timedelta(seconds=int(Time[x])) for x in range(len(Time))]
    VariableNames=FAAMCore.variables.keys() 
    
    #print(VariableNames)    
    FAAMCoreDict={}
    for NameStr in VariableNames: 
        if not (NameStr.endswith('_FLAG')) and (NameStr!='Time') :
#            print(NameStr)
            CoreVariable=np.array(FAAMCore[NameStr][:])
            if NameStr+'_FLAG' in VariableNames:
                CoreVariableFlag=np.array(FAAMCore[NameStr+'_FLAG'][:])
                CoreVariable[CoreVariableFlag!=0]=np.nan
            FAAMCoreDict[NameStr]=CoreVariable

    FAAMCoreDict['Time_Core']= Time_Core 
    
    TAT_DI_R_C=np.array(FAAMCore['TAT_DI_R'][:])
    TAT_DI_R_C-=273.15
    TAT_ND_R_C=np.array(FAAMCore['TAT_ND_R'][:])
    TAT_ND_R_C-=273.15
    
    FAAMCoreDict['TAT_ND_R_C']=TAT_ND_R_C
    FAAMCoreDict['TAT_DI_R_C']=TAT_DI_R_C  
    
    
    return FAAMCoreDict   
    
    #TAT_DI_R=np.array(FAAMCore['TAT_DI_R'][:])
    #TAT_DI_R_FLAG=np.array(FAAMCore['TAT_DI_R_FLAG'][:])
    #TAT_DI_R[TAT_DI_R_FLAG!=0]=np.nan
    #TAT_DI_R



    


#__________________________________________________________________________________________________


def LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile):

    CoreCloud = Dataset(CoreCloudPath+CoreCloudFile)
    
    #load bin sizes
    BinSize = np.loadtxt(CdpCalPath+CdpCalFile,skiprows=9,usecols=range(1,31),delimiter=',')
    CDP_BinCentre=BinSize[0][:]
    CDP_BinWidth=BinSize[2][:]
    #CDP_LogCentre=BinSize[4][:]
    #CDP_LogWidth=BinSize[6][:]

    #Load CDP
    CDP_FLAG= np.array(CoreCloud['CDP_FLAG'][:])
    for i in range(1, 31):
        if(i<10):
            i_string = 'CDP_0' + str(i)
        else:
            i_string = 'CDP_'+ str(i)
            # print(i_string)         
        CDP_bin = np.array(CoreCloud[i_string][:])
        CDP_bin[CDP_FLAG!=0]=np.nan
        if (i==1):
            CDP_Matrix=np.array(CDP_bin)
        else:
            CDP_Matrix=np.column_stack([CDP_Matrix, CDP_bin])

    FlightDate= datetime.datetime(int(CoreCloudFile[20:24]), int(CoreCloudFile[24:26]), int(CoreCloudFile[26:28]), 0, 0, 0)
    CDP_TSPM= np.array(CoreCloud['CDP_TSPM'][:])
    CDP_time_mid= [FlightDate + datetime.timedelta(seconds=int(CDP_TSPM[x])) for x in range(len(CDP_TSPM))]
    CDP_cm3=np.sum(CDP_Matrix,axis=1)
    CDP_dNdDp= CDP_Matrix / CDP_BinWidth
    CDP_gcm3=np.sum(CDP_Matrix*(4/3)*math.pi*(CDP_BinCentre/2e4)**3, axis=1)

    return CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre



#__________________________________________________________________________________________________


def ChangeBaseAvg(t1,w1,tdes, max_avg_window):
    wdes=np.zeros(len(tdes))*np.nan
    #max_avg_window=1 # max difference in time that could be averaged
    counts=np.zeros(len(wdes))
    if (len(w1)==len(t1)) and (len(wdes)==len(tdes)): # check wave lengths are correct   
        for i in range(len(w1)):        
            if(w1[i]==w1[i]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                if (np.sqrt((t1[i]-tdes[match])**2) < max_avg_window)	:	
                    if (counts[match]==0):
                        wdes[match]=w1[i]
                    else:
                        wdes[match]+=w1[i]
                    counts[match]+=1

        wdes/=counts   
        return(wdes)
    else:
        print("Array lengths incorrect")

#__________________________________________________________________________________________________

# t1_dt and tdes_dt are DateTime arrays. w1 is numpy array

def ChangeTimeBaseAvg(t1_dt,w1,tdes_dt, max_avg_window):
    
    t1=DateTime2IgorTime(t1_dt)
    tdes=DateTime2IgorTime(tdes_dt)
    
    wdes=np.zeros(len(tdes))*np.nan
    #max_avg_window=1 # max difference in time that could be averaged
    counts=np.zeros(len(wdes))
    if (len(w1)==len(t1)) and (len(wdes)==len(tdes)): # check wave lengths are correct   
        for i in range(len(w1)):        
            if(w1[i]==w1[i]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                if (np.sqrt((t1[i]-tdes[match])**2) < max_avg_window)	:	
                    if (counts[match]==0):
                        wdes[match]=w1[i]
                    else:
                        wdes[match]+=w1[i]
                    counts[match]+=1

        wdes/=counts   
        return(wdes)
    else:
        print("Array lengths incorrect")


#__________________________________________________________________________________________________

#x= time dimension
#y= e.g. size   
# uses x to sort        


def ChangeTimeBase2DAvg(t1_dt,w1,tdes_dt, max_avg_window):
    
    t1=DateTime2IgorTime(t1_dt)
    tdes=DateTime2IgorTime(tdes_dt)
    
    wdes=np.zeros((len(tdes),len(w1[0,:])))*np.nan
    #max_avg_window=1 # max difference in time that could be averaged
    counts=np.zeros((len(wdes),len(w1[0,:])))
    if (len(w1)==len(t1)) and (len(wdes)==len(tdes)): # check wave lengths are correct   
        for i in range(len(w1)):        
            if(w1[i,0]==w1[i,0]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                if (np.sqrt((t1[i]-tdes[match])**2) < max_avg_window)	:	
                    if (counts[match,0]==0):
                        wdes[match,:]=w1[i,:]
                    else:
                        wdes[match,:]+=w1[i,:]
                    counts[match,:]+=1

        wdes/=counts   
        return(wdes)
    else:
        print("Array lengths incorrect")
        
        
#__________________________________________________________________________________________________

#NevPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/Nevzorov/'
#NevName='c081_nevzorov_20180213_1hz_r0.nc'

def LoadNevzorov(NevPath,NevName):
    NevData = Dataset(NevPath+NevName)
    #print(NevData.variables.keys())
    #return NevData
    Time=np.array(NevData['TIME'][:])
    TWC_g_m3=np.array(NevData['TWC'][:])
    LWC_g_m3=np.array(NevData['LWC'][:])
    FlightDate= datetime.datetime(int(NevName[14:18]), int(NevName[18:20]), int(NevName[20:22]), 0, 0, 0)
    TimeNev= [FlightDate + datetime.timedelta(seconds=int(Time[x])) for x in range(len(Time))]
    TimeNevIgor=DateTime2IgorTime(TimeNev)
    
       
    return TWC_g_m3, LWC_g_m3, TimeNev, TimeNevIgor

#__________________________________________________________________________________________________

def LoadOAP(FilePath,FileNameOAP): 
    
    DataOAP = h5py.File(FilePath+FileNameOAP, 'r')
    NC_all_x=np.array(DataOAP['NC_all_x'])
    NC_all_y=np.array(DataOAP['NC_all_y'])
    NC_all_z=np.array(DataOAP['NC_all_z'])    
    NC_S_z=np.array(DataOAP['NC_S_z']) 
    NC_LI_z=np.array(DataOAP['NC_LI_z']) 
    NC_MI_z=np.array(DataOAP['NC_MI_z']) 
    NC_HI_z=np.array(DataOAP['NC_HI_z'])
    NC_All_accept_CH0_z=np.array(DataOAP['NC_All_accept_CH0_z'])
    NC_All_accept_CH1_z=np.array(DataOAP['NC_All_accept_CH1_z'])
    DataOAP.close()
    
    NC_All_accept_total=np.sum(NC_All_accept_CH0_z,axis=1)
    NC_HI_total=np.sum(NC_HI_z,axis=1)
    NC_MI_total=np.sum(NC_MI_z,axis=1)
    NC_LI_total=np.sum(NC_LI_z,axis=1)

    #Calculate DateTime
    NC_sMidnight=(NC_all_x[:-1:1]+NC_all_x[1::1])/2
    FlightDate= datetime.datetime(int(FileNameOAP[5:9]), int(FileNameOAP[9:11]), int(FileNameOAP[11:13]), 0, 0, 0)
    NC_DateTime= [FlightDate + datetime.timedelta(seconds=int(NC_sMidnight[x])) for x in range(len(NC_sMidnight))]

    return NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_LI_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime


#__________________________________________________________________________________________________

def round_time(dt=None, dateDelta=datetime.timedelta(minutes=1), to='average'):
    """Round a datetime object to a multiple of a timedelta
    dt : datetime.datetime object, default now.
    dateDelta : timedelta object, we round to a multiple of this, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
            Stijn Nevens 2014 - Changed to use only datetime objects as variables
    """
    round_to = dateDelta.total_seconds()

    if dt == None : dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    # // is a floor division, not a comment on following line:
    if to == 'up':
        # // is a floor division, not a comment on following line (like in javascript):
        rounding = (seconds + round_to) // round_to * round_to
    elif to == 'down':
        rounding = seconds // round_to * round_to
    else:
        rounding = (seconds + round_to / 2) // round_to * round_to
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)


#__________________________________________________________________________________________________

def Average_nPts(Array,nPts):
    
    ArrayAvg=np.nanmean(np.pad(Array.astype(float), (0, nPts - Array.size%nPts), mode='constant', constant_values=np.NaN).reshape(-1, nPts), axis=1)
    return ArrayAvg

#__________________________________________________________________________________________________

def Average_nPts_datetime(DateTimeArray,nPts):
        
    Array=np.asarray(DateTime2IgorTime(DateTimeArray))
    ArrayAvg=np.nanmean(np.pad(Array.astype(float), (0, nPts - Array.size%nPts), mode='constant', constant_values=np.NaN).reshape(-1, nPts), axis=1)
    DateTimeArrayAvg=Igor2DateTime(ArrayAvg)
    
    return DateTimeArrayAvg


#__________________________________________________________________________________________________
# average every nPts points in x dimension. y dimension size remains the same


def Average_nPts_2D(Array,nPts):
    
    nx=float(Array.shape[0])
    ny=Array.shape[1]

    
    ArrayAvg=np.zeros((int(np.ceil(nx/nPts)),ny))*np.nan
    
    for i in range(int(np.ceil(nx/nPts))-1):
        Slice=Array[i*nPts:(i+1)*nPts,:]
        #SliceAvg=np.nanmean(Slice,axis=0)
        ArrayAvg[i,:]=np.nanmean(Slice,axis=0)
    
    return ArrayAvg
    
#__________________________________________________________________________________________________












