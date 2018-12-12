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
import scipy.io as sio
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.dates import DateFormatter
import os
import pandas as pd

#__________________________________________________________________________________________________
#Bin DataY based on DataX value


#def BinLinear(DataY,DataX,StartBin,EndBin, NBins):
#    Output=np.zeros((len(DataX),NBins))*np.nan
#    bins = np.linspace(StartBin,EndBin, NBins)
#    
#    digitized = np.digitize(DataX, bins)
#    
#    for i in range(len(DataY)-1):
#        if DataX[i]>=StartBin and DataX[i]<EndBin :
#            Output[i,digitized[i]]=DataY[i]
#
#    return Output, bins

#__________________________________________________________________________________________________

#Bin DataY based on DataX value. Now returns lower, upper and mid bins 


def BinLinearV2(DataY,DataX,StartBin,EndBin, NBins):
    
    BinWidth=(EndBin-StartBin)/ NBins
    
    Output=np.zeros((len(DataX),NBins))*np.nan
    
    BinsLower = np.linspace(StartBin,EndBin-BinWidth, NBins) # LowerBin
    BinsUpper = np.linspace(StartBin+BinWidth,EndBin, NBins) # UpperBin
    BinsMid= (BinsUpper + BinsLower) / 2

    digitized = np.digitize(DataX, BinsUpper)
    
    for i in range(len(DataY)-1):
        if DataX[i]>=StartBin and DataX[i]<EndBin :
            Output[i,digitized[i]]=DataY[i]

    return Output, BinsLower, BinsUpper, BinsMid






#__________________________________________________________________________________________________


def BinAverage2D(DataX,DataY,DataZ, StartBinX,EndBinX, NBinsX,StartBinY,EndBinY, NBinsY):
    
    binsX = np.linspace(StartBinX,EndBinX, NBinsX)
    digitizedX = np.digitize(DataX, binsX)

    binsY = np.linspace(StartBinY,EndBinY, NBinsY)
    digitizedY = np.digitize(DataY, binsY)

    Output=np.zeros((NBinsX,NBinsY))*np.nan
    Counts=np.zeros((NBinsX,NBinsY))
    
    for i in range(len(DataZ)-1):        
        if Counts[digitizedX[i],digitizedY[i]] == 0 :               # the -1 means that bins[i] <= x < bins[i+1]
            Output[digitizedX[i],digitizedY[i]]=DataZ[i]            
            Counts[digitizedX[i],digitizedY[i]]+=1
        else:
            Output[digitizedX[i],digitizedY[i]]+=DataZ[i]
            Counts[digitizedX[i],digitizedY[i]]+=1

    Output/=Counts
    return Output, binsX, binsY
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
        for i in range(len(w1)-1):        
            if(w1[i]==w1[i]):
                match=bisect.bisect_left(tdes, t1[i]) # assume is sorted, which it should be 
                
                if (match < len(tdes)-1):
                    if (np.absolute(t1[i]-tdes[match]) < max_avg_window)	:	
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
#

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a)) 
    # Radius of earth in kilometers is 6371
    km = 6371* c
    return km

#__________________________________________________________________________________________________


def BinMid_2_Width(BinMid):
    
    Width=np.zeros(len(BinMid))
    WidthTmp=(BinMid[1:-1:1]-BinMid[0:-2:1])/2+(BinMid[2::1]-BinMid[1:-1:1])/2
    Width[1:-1]=WidthTmp
    Width[0]=WidthTmp[0]
    Width[len(Width)-1]=WidthTmp[len(WidthTmp)-1]

    #Width=np.zeros(len(BinMid))
    return Width

#__________________________________________________________________________________________________

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    
    from: `StackOverflow <http://stackoverflow.com/questions/7008608/scipy-io-loadmat-nested-structures-i-e-dictionaries>`_
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

#__________________________________________________________________________________________________

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        
#__________________________________________________________________________________________________

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict



#print_mat_nested(matdata, nkeys=6)

#__________________________________________________________________________________________________

def Matlab2PythonTime(matlab_datenum):
    python_datetime = datetime.datetime.fromordinal(int(matlab_datenum)) + datetime.timedelta(days=matlab_datenum%1) - datetime.timedelta(days = 366)
    return python_datetime

#__________________________________________________________________________________________________

def Matlab2PythonTimeArray(matlab_datenum_array):
    F_array = np.vectorize(Matlab2PythonTime)
    PythonTimeArray=F_array(matlab_datenum_array)

    return PythonTimeArray

#__________________________________________________________________________________________________

#plt.xlim(['2018-02-13 06:35:00','2018-02-13 06:36:00'])


def TimeSeriesPlot(TimeArray, YArray, TimeStart, TimeEnd, Ymin, Ymax):
    
    
    fig=plt.figure(figsize=(10,5)) 
    plt.rcParams.update({'font.size': 10})
    formatter = DateFormatter('%H:%M')
    
    plt.subplot(1, 1, 1)
    #plt.title(FlightNumber)
    
    plt.plot(TimeArray, YArray)
    
    plt.gca().xaxis.set_major_locator(mdates.HourLocator() )
    plt.gca().xaxis.set_minor_locator(mdates.MinuteLocator(30) )
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xlabel('time') 
    #plt.xlim(['2018-02-13 06:35:00','2018-02-13 06:36:00'])
    plt.xlim([TimeStart,TimeEnd])
    plt.ylim([Ymin, Ymax])


#__________________________________________________________________________________________________

#Diameter correction from Reuter and Bakan 1997

def ReuterBakanDiameter(Diameter25,AreaFraction0): 

    C = 0.962
    D = -10.34
    E = 1.444     
    F = 7.036
    
    D_ReuterBakan=np.zeros(len(Diameter25))*np.nan
    for i in range(len(Diameter25)) :
        #D_ReuterBakan[i] = (C*Diameter25[i] + D) * AreaFraction0[i] + (E*Diameter[i] + F)
        D_ReuterBakan[i] = (Diameter25[i] - D - F*AreaFraction0[i]) / (E*AreaFraction0[i] + C)
    
    return D_ReuterBakan




#_______________________________________________________________________________________   


    
def KorolevCorrectedD(FilledArea, VoidArea,Diameter):
    
    CorrFactArray=np.empty(len(Diameter), dtype=float)
    Dspot_Dmax, Dmax_D0=GetKorolevRatios()
    for i in range(len(Diameter)):
        CorrFactArray[i]=KorolevCorrection(Dspot_Dmax, Dmax_D0,VoidArea[i],FilledArea[i])
    
    D_KorolevCorr=Diameter*CorrFactArray
    
    return D_KorolevCorr
#_______________________________________________________________________________________   
  
    
def KorolevCorrection(Dspot_Dmax, Dmax_D0,VoidAreaElement,FilledAreaElement):  
    #KorolevRatios= KorolevCorrection()
    #KorolevRatiosArray=KorolevRatios.values
    #Dspot_Dmax=KorolevRatiosArray[:,3]
    #Dmax_D0=KorolevRatiosArray[:,1]
    if VoidAreaElement>0 and FilledAreaElement>0:
        scaleF = 4/math.pi
        PixelRatio = math.sqrt(VoidAreaElement*scaleF)/math.sqrt(FilledAreaElement*scaleF)
        pos=(np.abs(Dspot_Dmax-PixelRatio)).argmin()
        CorrFac = 1/Dmax_D0[pos]
    else:
        CorrFac=1
    
    return CorrFac
    #pos = binarysearchinterp(Dspot_Dmax,pixelRatio)
#_______________________________________________________________________________________    

def GetKorolevRatios():
    KorolevRatios = pd.read_csv('C:/Users/Admin TEMP/Documents/DropletGun/Korolev07_ratios.csv')
    KorolevRatiosArray=KorolevRatios.values
    Dspot_Dmax=KorolevRatiosArray[:,3]
    Dmax_D0=KorolevRatiosArray[:,1]
    
    return Dspot_Dmax, Dmax_D0 

#_________________________________________________________________________________________________  

def GetKorolevRatios_Zd():
    KorolevRatios = pd.read_csv('C:/Users/Admin TEMP/Documents/DropletGun/Korolev07_ratios_withZd.csv')
    KorolevRatiosArray=KorolevRatios.values
    
    Z_d=KorolevRatiosArray[:,0]
    Dspot_Dmax=KorolevRatiosArray[:,1]
    Dspot_Dimg=KorolevRatiosArray[:,2]
    Dmax_D0=KorolevRatiosArray[:,3]
    Dimg_Dmax=KorolevRatiosArray[:,4]
    
    return Z_d,Dspot_Dmax,Dspot_Dimg,Dmax_D0,Dimg_Dmax


#_______________________________________________________________________________________  
 

# correction factor for diameter derived from areafraction2 
    
def Level2_D_Correction_vector(Diameter, AreaFraction2,VoidRatio) :
    
    CorrFactArray=np.empty(len(Diameter), dtype=float)
    Z_d,Dspot_Dmax,Dspot_Dimg,Dmax_D0,Dimg_Dmax= GetKorolevRatios_Zd()
    
    for i in range(len(Diameter)):
        CorrFactArray[i]=Level2_D_Correction(Z_d, Dmax_D0, AreaFraction2[i],VoidRatio[i])
    
    D_AreaFraction2Corr=Diameter*CorrFactArray
    
    return D_AreaFraction2Corr
#_______________________________________________________________________________________   
  
  
    
def Level2_D_Correction(LookUp_Zd,Dmax_D0, AreaFraction2, VoidRatio):  

    a= 1.02254
    b= -1.08683
    c= 0.840135
    
    #LookUp_areafraction2= a*LookUp_Zd**2 + b*LookUp_Zd + c 
    LookUp_areafraction2=(np.where(LookUp_Zd<0.5,a*LookUp_Zd**2 + b*LookUp_Zd + c,0))
    
    
    if AreaFraction2 > 0.55 and AreaFraction2 < 0.8 and VoidRatio == 0 :
        pos=(np.abs(LookUp_areafraction2-AreaFraction2)).argmin()
        CorrFac = 1/Dmax_D0[pos]
    else:
        CorrFac=1
    
    return CorrFac

#_______________________________________________________________________________________   




