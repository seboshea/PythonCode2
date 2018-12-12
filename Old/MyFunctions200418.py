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


def ChangeBaseAvg(t1,w1,tdes):
    wdes=np.zeros(len(tdes))*np.nan
    max_avg_window=1 # max difference in time that could be averaged
    counts=np.zeros(len(wdes))
    if (len(w1)==len(t1)) and (len(wdes)==len(tdes)): # check wave lengths are correct   
        for i in range(len(w1)):        
            if(w1[i]==w1[i]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                if (np.sqrt((t1[i]-tdes[match])^2) < max_avg_window)	:	
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


def ChangeTimeBaseAvg(t1,w1,tdes):
    wdes=np.zeros(len(tdes))*np.nan
    max_avg_window=1 # max difference in time that could be averaged
    counts=np.zeros(len(wdes))
    if (len(w1)==len(t1)) and (len(wdes)==len(tdes)): # check wave lengths are correct   
        for i in range(len(w1)):        
            if(w1[i]==w1[i]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                if (np.sqrt(((t1[i]-tdes[match]).seconds)^2) < max_avg_window)	:	
                    if (counts[match]==0):
                        wdes[match]=w1[i]
                    else:
                        wdes[match]+=w1[i]
                    counts[match]+=1

        wdes/=counts   
        return(wdes)
    else:
        print("Array lengths incorrect")

