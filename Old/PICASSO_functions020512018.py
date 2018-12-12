# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 11:00:34 2018

@author: Admin TEMP
"""
import glob, os
import datetime
import numpy as np
#from netCDF4 import Dataset
from MyFunctions import loadTAS,LoadCoreCloud,loadFAAMCore, DateTime2IgorTime, FindFirstGreaterThan,is_number  
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.dates import DateFormatter
import pandas as pd
from scipy.interpolate import interp1d
from pylab import figure, cm
from matplotlib.colors import LogNorm
import h5py
import bisect

#_________________________________________________________________________________________________

def LoadFlightSummary(FilePath,FileName):
    
    #FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/'
    #FileName='flight-sum_faam_20180214_r0_c082_fmt.csv'
    text_file = open(FilePath+FileName, "r")

    #Summary=csv.reader(text_file, delimiter=' ')
    lines = text_file.readlines()
    Summary_St=np.zeros(len(lines))* np.nan
    Summary_End=np.zeros(len(lines))* np.nan
    Summary_Name=[''] * len(lines)
    for x in range(2,len(lines)):
        Tmp=lines[x].split(',')
        if (is_number(Tmp[0])):
            Fudge=Tmp[0]   
            Summary_St[x]=int(Fudge[0:2])*3600+int(Fudge[2:4])*60+int(Fudge[4:6]) # incredibly convoluted way to do this!!!
        else:
            Summary_St[x]=np.nan       
        if (is_number(Tmp[1])):
            Fudge=Tmp[1]
            Summary_End[x]=int(Fudge[0:2])*3600+int(Fudge[2:4])*60+int(Fudge[4:6]) # incredibly convoluted way to do this!!!
        else:
            Summary_End[x]=np.nan
        Summary_Name[x]=Tmp[2]
    text_file.close()
    return Summary_St, Summary_End, Summary_Name

#_________________________________________________________________________________________________

def PSD_runs():
    FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/'
    FileName='flight-sum_faam_20180214_r0_c082_fmt.csv'
    
    #Location of 2DS, CIP100, HVPS, CIP15 data files
    DataFilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C082/'
    
    Summary_St, Summary_End, Summary_Name=LoadFlightSummary(FilePath,FileName)
    #for x in range(14,17):
    for x in range(len(Summary_St)):
        if ((Summary_St[x]==Summary_St[x]) and (Summary_End[x]==Summary_End[x]) and (len(Summary_Name[x])>0)):
            StartTime=datetime.datetime(int(FileName[16:20]), int(FileName[20:22]), int(FileName[22:24]), 0, 0, 0)+datetime.timedelta(seconds=int(Summary_St[x]))
            EndTime=datetime.datetime(int(FileName[16:20]), int(FileName[20:22]), int(FileName[22:24]), 0, 0, 0)+datetime.timedelta(seconds=int(Summary_End[x]))
            #PSD_comparison(StartTime,EndTime,DataFilePath,Summary_Name[x]+'_LI','NC_LI')
            PSD_comparison(StartTime,EndTime,DataFilePath,Summary_Name[x],'NC_All_accept_CH0')

#    StartTime=datetime.datetime(2018, 2, 14, 17, 20, 0)
#    EndTime=datetime.datetime(2018, 2, 14, 17, 30, 0)

#    Figurename='test'


#_________________________________________________________________________________________________
def PSD_comparison(StartTime,EndTime,FilePath,Figurename,Prefix):            
    #Prefix='NC_All_accept'
    
    #2DS
    FileName2DS = 'C082_20180214_2DS_AllAcceptCentre.h5'
    #FileName2DS = 'C082_20180214_2DS_AllIn.h5'
    PSD_name=Prefix+'_z'
    NC_SizeMid_2DS, NC_z_dNdDp_2DS= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileName2DS,PSD_name)
    
    #HVPS
    FileNameHVPS = 'C082_20180214_HVPS_AllIn.h5'
    PSD_name='NC_All_accept_CH1_z'
    #PSD_name=Prefix+'_z'
    NC_SizeMid_HVPS, NC_z_dNdDp_HVPS= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileNameHVPS,PSD_name)
    
    #3V-2DS
    #FileName3V = 'C082_20180214_3V2DS_AllIn.h5'
    FileName3V = 'C082_20180214_3V2DS_CentreIn.h5'
    #PSD_name=Prefix+'_z'
    PSD_name='NC_All_accept_CH0_z'
    NC_SizeMid_3VCHO, NC_z_dNdDp_3VCH0= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileName3V,PSD_name)
    PSD_name='NC_All_accept_CH1_z'
    NC_SizeMid_3VCH1, NC_z_dNdDp_3VCH1= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileName3V,PSD_name)
    
    #CIP100
    FileNameCIP100 = 'C082_20180214_CIP100_GS_AllIn.h5'
    PSD_name=Prefix+'_z'
    NC_SizeMid_CIP100, NC_z_dNdDp_CIP100= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileNameCIP100,PSD_name)
    
    #CIP15
    FileNameCIP15 = 'C082_20180214_CIP15_GS_AllAcceptCentre.h5'
    #FileNameCIP15 = 'C082_20180214_CIP15_GS_AllIn.h5'
    PSD_name=Prefix+'_z'
    NC_SizeMid_CIP15, NC_z_dNdDp_CIP15= Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileNameCIP15,PSD_name)
    
    #CDP
    CoreCloudPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    CoreCloudFile='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    CdpCalPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
    CdpCalFile='C082_water_diameters.csv'
    CDP_BinCentre, CDP_dNdDp_L_avg=Avg_CDP_dndDp(StartTime,EndTime,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    
    #Temperature
    FullPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    FAAMCoreName='core_faam_20180214_v004_r0_c082_1hz.nc'    
    FAAMCoreDict=loadFAAMCore(FullPath,FAAMCoreName)
    Time_Core=FAAMCoreDict['Time_Core']
    TAT_DI_R=FAAMCoreDict['TAT_DI_R']
    TAT_ND_R=FAAMCoreDict['TAT_ND_R']
    StartIdx=bisect.bisect_left(Time_Core, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Core, EndTime) # assume time is sorted, which it should be
    TAT_DI_R_avg=np.around(np.nanmean(TAT_DI_R[StartIdx:EndIdx])-273.15, decimals=1)
    TAT_ND_R_avg=np.around(np.nanmean(TAT_ND_R[StartIdx:EndIdx])-273.15, decimals=1)
    
    #Plot    
    fig=plt.figure(figsize=(10,10)) 
    fig, ax = plt.subplots()
    plt.title('Temperature DI='+str(TAT_DI_R_avg)+'°C, ND='+ str(TAT_ND_R_avg)+'°C')
    plt.plot(NC_SizeMid_2DS,NC_z_dNdDp_2DS, label='2DS')
    plt.plot(NC_SizeMid_HVPS,NC_z_dNdDp_HVPS, label='HVPS')
    plt.plot(NC_SizeMid_CIP100,NC_z_dNdDp_CIP100, label='CIP100')
    plt.plot(NC_SizeMid_CIP15,NC_z_dNdDp_CIP15, label='CIP15')
    plt.plot(NC_SizeMid_3VCHO,NC_z_dNdDp_3VCH0, label='3VCH0')
    plt.plot(NC_SizeMid_3VCH1,NC_z_dNdDp_3VCH1, label='3VCH1')
    plt.plot(CDP_BinCentre,CDP_dNdDp_L_avg, label='CDP')
    plt.ylabel('dNdDp, L$^{-1}$ μm$^{-1}$')
    plt.xlabel('Diameter, μm')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
        
    #plt.show()
    plt.savefig(FilePath+Figurename+'.png',dpi=200)
    plt.close(fig)
    

#_________________________________________________________________________________________________
   
def Avg_CDP_dndDp(StartTime,EndTime,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile):
   

    CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    
    #Search for start and stop times
    StartIdx=bisect.bisect_left(CDP_time_mid, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(CDP_time_mid, EndTime) # assume time is sorted, which it should be

    print(StartIdx)
    print(EndIdx)

    CDP_dNdDp_L_avg=np.nanmean(CDP_dNdDp[StartIdx:EndIdx][:],axis=0)
    
    CDP_dNdDp_L_avg*=1000 # Convert to L
    
    return CDP_BinCentre, CDP_dNdDp_L_avg
    
#_________________________________________________________________________________________________
def Avg_OAP_dndDp(StartTime,EndTime,FilePath,FileName,PSD_name):            

    Data_h5 = h5py.File(FilePath+FileName, 'r')              
    NC_all_x=np.array(Data_h5['NC_all_x'])
    NC_all_y=np.array(Data_h5['NC_all_y'])
    NC_z=np.array(Data_h5[PSD_name])
    
    
    #NC_all_z=np.array(Data2DS['NC_all_z'])    
    #NC_S_z=np.array(Data2DS['NC_S_z']) 
    #NC_LI_z=np.array(Data2DS['NC_LI_z']) 
    #NC_MI_z=np.array(Data2DS['NC_MI_z']) 
    #NC_HI_z=np.array(Data2DS['NC_HI_z'])
    #NC_All_accept_CH0_z=np.array(Data2DS['NC_All_accept_CH0_z'])

    Data_h5.close()
    
    #NC_All_accept_total=np.sum(NC_All_accept_CH0_z,axis=1)
    #NC_HI_total=np.sum(NC_HI_z,axis=1)
    #NC_MI_total=np.sum(NC_MI_z,axis=1)
    #NC_LI_total=np.sum(NC_LI_z,axis=1)
  
    NC_sMidnight=(NC_all_x[:-1:1]+NC_all_x[1::1])/2
    FlightDate= datetime.datetime(int(FileName[5:9]), int(FileName[9:11]), int(FileName[11:13]), 0, 0, 0)
    NC_DateTime= [FlightDate + datetime.timedelta(seconds=int(NC_sMidnight[x])) for x in range(len(NC_sMidnight))]

    #Search for start and stop times
    StartIdx=bisect.bisect_left(NC_DateTime, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(NC_DateTime, EndTime) # assume time is sorted, which it should be 

    NC_z_PSD=np.nanmean(NC_z[StartIdx:EndIdx][:],axis=0)

    NC_SizeMid=(NC_all_y[:-1:1]+NC_all_y[1::1])/2
    NC_SizeWidth=(NC_all_y[1::1]-NC_all_y[:-1:1])
    NC_z_dNdDp=NC_z_PSD/NC_SizeWidth
    
    return NC_SizeMid, NC_z_dNdDp

#_________________________________________________________________________________________________

def LoadOAPFiles():
    OAPPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/'

    #FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/'
    #FileNameOAP = 'C082_20180214_2DS.h5'
    
    for filena in os.listdir(OAPPath):
        if filena.endswith(".h5"):
            print(filena)
            PlotOAP(OAPPath,filena)    

#_________________________________________________________________________________________________

def PlotOAP(FilePath,FileNameOAP): 
    
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

    fig=plt.figure(figsize=(10,10))
    formatter = DateFormatter('%H:%M')
    ax1=fig.add_subplot(4, 1, 1)
    plt.title(FileNameOAP)   

    plt.ylabel('All accept')  # we already handled the x-label with ax1
    plt.plot(NC_DateTime, NC_All_accept_total, '-')
             
    plt.subplot(4, 1, 2, sharex=ax1)
    plt.plot(NC_DateTime, NC_HI_total, '-')
    plt.ylabel('HI, L$^{-1}$')
    
    fig.add_subplot(4, 1, 3, sharex=ax1)
    plt.plot(NC_DateTime, NC_MI_total, '-')
    plt.ylabel('MI, L$^{-1}$')
    
    fig.add_subplot(4, 1, 4, sharex=ax1)
    plt.plot(NC_DateTime, NC_LI_total, '-')
    plt.ylabel('LI, L$^{-1}$')
      
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xlabel('time')
    plt.gca().set_xlim([NC_DateTime[0],NC_DateTime[-1]]) 

    plt.tight_layout()       
    Figurename=FilePath+FileNameOAP.replace('.h5','.png')
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)


#_________________________________________________________________________________________________

def BatchLoadTAS():
    BatchTAS_RVSM=[]
    BatchDateTime=[]


   # CorePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/'
   # for filena in os.listdir(CorePath):
   #     FullPath=CorePath+filena+'/core_processed/'
    FullPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'    
    for filena in os.listdir(FullPath):
        if filena.endswith("1hz.nc"):     
            print(FullPath)
            print(filena)
            TAS_RVSM, DateTime = loadTAS(FullPath,filena)
            BatchTAS_RVSM= np.append(BatchTAS_RVSM,TAS_RVSM)
            BatchDateTime= np.append(BatchDateTime,DateTime) 
    plt.plot(BatchDateTime,BatchTAS_RVSM)
    np.savetxt(FullPath+'AllFlightsTAS.txt', np.column_stack([BatchDateTime,BatchTAS_RVSM]), delimiter=',')
    
#__________________________________________________________________________________________________

def BatchPlotCDPCore():
        
    DataPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/'
    #Flights2Plot=['c075', 'c076', 'c077','c078','c079','c080','c081','c082']
    Flights2Plot=['c075', 'c076', 'c077','c078','c079','c080','c081','c082']
    #Flights2Plot=['c082']
    FullPath=DataPath+'AllCore/'
    CoreCloudPath=DataPath+'AllCoreCloud/'
    CdpCalPath=DataPath+'AllCDPcals/'
    
    
    for FlightNumber in Flights2Plot:
        CorePathAndName=glob.glob(FullPath+'*'+FlightNumber+'_1hz.nc')
        FAAMCoreName=os.path.basename(CorePathAndName[0])
        #print(FAAMCoreName)
        CorePathAndName=glob.glob(CoreCloudPath+'core-cloud-phy_faam_*'+FlightNumber+'.nc')
        CoreCloudFile=os.path.basename(CorePathAndName[0])
        #print(CoreCloudFile)
        CorePathAndName=glob.glob(CdpCalPath+FlightNumber+'_water_diameters.csv')
        CdpCalFile=os.path.basename(CorePathAndName[0])
        #print(CdpCalFile)
        LoadCoreCDP(FullPath,FAAMCoreName,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)

#FullPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
#FAAMCoreName='core_faam_20180214_v004_r0_c082_1hz.nc'    

 #CoreCloudPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
 #CoreCloudFile='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    
 #CdpCalPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
 #CdpCalFile='C082_water_diameters.csv'

#__________________________________________________________________________________________________

    
def LoadCoreCDP(FullPath,FAAMCoreName,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile):  
    # Load core data
    #FullPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    #FAAMCoreName='core_faam_20180214_v004_r0_c082_1hz.nc'    
    #global FAAMCoreDict
    FAAMCoreDict=loadFAAMCore(FullPath,FAAMCoreName)
    
    Time_Core=FAAMCoreDict['Time_Core']
    TAT_DI_R=FAAMCoreDict['TAT_DI_R']
    TAT_DI_R_C=TAT_DI_R-273.15
    TAT_ND_R=FAAMCoreDict['TAT_ND_R']
    TAT_ND_R_C=TAT_ND_R-273.15
    ALT_GIN=FAAMCoreDict['ALT_GIN']
    PS_RVSM=FAAMCoreDict['PS_RVSM']
    Time_Core_Igor=DateTime2IgorTime(Time_Core)
    
    #Load CDP data
    #CoreCloudPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    #CoreCloudFile='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    #CdpCalPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
    #CdpCalFile='C082_water_diameters.csv'
    CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    CDP_time_mid_Igor=DateTime2IgorTime(CDP_time_mid)    
    
    # Change CDP to core timebase 
    CDP_cm3_core=np.interp(Time_Core_Igor,CDP_time_mid_Igor,CDP_cm3)
    CDP_gcm3_core=np.interp(Time_Core_Igor,CDP_time_mid_Igor,CDP_gcm3)
    CDP_gm3_core=CDP_gcm3_core*1e6
    
    RadialFlag=CreateRadialFlag(FAAMCoreDict)
    
    PlotCoreCDP(FAAMCoreName,CoreCloudPath,Time_Core, ALT_GIN,TAT_DI_R_C,CDP_gm3_core, CDP_cm3_core,TAT_ND_R_C,RadialFlag)
     
    
    
#__________________________________________________________________________________________________

def PlotCoreCDP(FAAMCoreName,CoreCloudPath,Time_Core, ALT_GIN,TAT_DI_R_C,CDP_gm3_core, CDP_cm3_core,TAT_ND_R_C,RadialFlag):
    
    fig=plt.figure(figsize=(10,10))
    
    formatter = DateFormatter('%H:%M')
    ax1=fig.add_subplot(4, 1, 1)
    plt.title(FAAMCoreName[0:31])



    #plt.gca().set_ylim([0, 8000]) 

    plt.ylabel('245° radial', color='tab:red')  # we already handled the x-label with ax1
    plt.fill(Time_Core, RadialFlag, color='tab:red')
    plt.tick_params(axis='245° radial', labelcolor='tab:red')
    
    ax2 = ax1.twinx()
    ax2.plot(Time_Core, ALT_GIN, '-',color='tab:blue')
    ax2.set_ylabel('Altitude, m', color='tab:blue')
    
          
    plt.subplot(4, 1, 2, sharex=ax1)
    plt.plot(Time_Core, TAT_ND_R_C, '-')
    plt.plot(Time_Core, TAT_DI_R_C, '-')
    plt.gca().set_ylim([-45, 10])
    plt.ylabel('Temperature, °C')
    
    fig.add_subplot(4, 1, 3, sharex=ax1)
    plt.plot(Time_Core, CDP_cm3_core, '-')
    plt.gca().set_ylim([0, 120])
    plt.ylabel('CDP, cm$^{-3}$')
    
    fig.add_subplot(4, 1, 4, sharex=ax1)
    plt.plot(Time_Core, CDP_gm3_core, '-')
    plt.gca().set_ylim([0, 0.7])
    plt.ylabel('CDP, g m$^{-3}$')
      
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xlabel('time')
    StartPt=FindFirstGreaterThan(ALT_GIN,150)
    plt.gca().set_xlim([Time_Core[StartPt],Time_Core[-1]]) #C082
    #plt.gca().set_xlim([Time_Core[StartPt], datetime.datetime(2018, 2, 14,20,30,00)]) #C082

    plt.tight_layout()       
    Figurename=CoreCloudPath+'CDPplot_'+FAAMCoreName[10:31]+'.png'
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)
#__________________________________________________________________________________________________
    
    
def PlotCDP():
    
    CoreCloudPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    CoreCloudFile='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    CdpCalPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
    CdpCalFile='C082_water_diameters.csv'
    
    
            
    CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    #print(CDP_BinCentre)
    
    #formatter = DateFormatter('%H:%M')
    #plt.subplot(2, 1, 1)
    #plt.plot(CDP_time_mid, CDP_cm3, '-')
    #plt.title('C082')
    #plt.ylabel('CDP, cm3')
    #plt.gcf().axes[0].xaxis.set_major_formatter(formatter)


    #plt.subplot(2, 1, 2)
    #plt.plot(CDP_time_mid, CDP_gcm3, '-')
    #plt.xlabel('time')
    #plt.ylabel('CDP_gcm3')
    #plt.gcf().axes[1].xaxis.set_major_formatter(formatter)
    
    #plt.show()
        
    #plot dNdDp
    x_lims = mdates.date2num(CDP_time_mid)
    Min_dNdDp=np.nanmin(CDP_dNdDp[CDP_dNdDp>0])
    print(Min_dNdDp)
    Max_dNdDp=np.nanmax(CDP_dNdDp[CDP_dNdDp>0])
    print(Max_dNdDp)
    fig, ax = plt.subplots()
    ax.imshow(CDP_dNdDp,extent = [x_lims[0], x_lims[len(x_lims)-1], CDP_BinCentre[0] , CDP_BinCentre[-1]], aspect='auto',norm=colors.LogNorm(Min_dNdDp, Max_dNdDp), cmap='PuBu_r')
    #ax.imshow(x_lims, CDP_BinCentre,CDP_dNdDp, aspect='auto',extent = [x_lims[0], x_lims[len(x_lims)-1],  4, 55] )
    
    #ax.imshow(CDP_time_mid, CDP_BinCentre, CDP_dNdDp, extent = [x_lims[0], x_lims[-1],  4, 55], aspect='auto')
    ax.xaxis_date()
    date_format = mdates.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(date_format)
    fig.autofmt_xdate()
    plt.yscale('log')
    #plt.colorbar(CDP_dNdDp, orientation='horizontal')
    plt.show()
    #return CDP_dNdDp
    
    
#_____________________________________________________________________________________________________

def CreateRadialFlag(FAAMCoreDict):
    HDG_GIN=np.array(FAAMCoreDict['HDG_GIN'])
    LAT_GIN=np.array(FAAMCoreDict['LAT_GIN'])
    LON_GIN=np.array(FAAMCoreDict['LON_GIN'])
       
    RadialFlag=np.where(((HDG_GIN>=230) & (HDG_GIN<=251)) | ((HDG_GIN>=60) & (HDG_GIN<=80)),1,0)
    RadialFlag=np.where(((LAT_GIN>=50.81) & (LAT_GIN<=51.15)) & ((LON_GIN>=-2.6) & (LON_GIN<=-1.4)),RadialFlag,0)
    
    #plt.plot(HDG_GIN, color='tab:red')
    
    
    #plt.fill(RadialFlag)
    #plt.plot(RadialFlag)
    
    
    return RadialFlag






    