# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 11:00:34 2018

@author: Admin TEMP
"""
import glob, os
import datetime
import numpy as np
#from netCDF4 import Dataset
from MyFunctions import ChangeTimeBase2DAvg,ChangeTimeBaseAvg,loadTAS,LoadCoreCloud,loadFAAMCore, DateTime2IgorTime, FindFirstGreaterThan,is_number
from MyFunctions import LoadNevzorov,LoadOAP,round_time, Average_nPts_datetime,Average_nPts,Average_nPts_2D
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

# Create dictionary with all 1hz data on core time base. Option to include psds but there on a different time base


def DataPaths():
    PathDict={}
   
    PathDict['C098','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C098','FAAMCoreName']='core_faam_20180424_v004_r0_c098_1hz.nc'
    PathDict['C098','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C098','CoreCloudFile']='core-cloud-phy_faam_20180424_v501_r0_c098.nc.nc'
    PathDict['C098','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C098','CdpCalFile']='c098_water_diameters.csv'
    PathDict['C098','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C098','NevName']=''
    PathDict['C098','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/'
    PathDict['C098','FileName2DS']='C098_20180424_2DS_Allin.h5'
    PathDict['C098','FileNameHVPS'] = 'C098_20180424_HVPS_AllIn.h5' 
    
    PathDict['C097','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C097','FAAMCoreName']='core_faam_20180423_v004_r0_c097_1hz.nc'
    PathDict['C097','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C097','CoreCloudFile']='core-cloud-phy_faam_20180423_v501_r0_c097.nc'
    PathDict['C097','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C097','CdpCalFile']='c097_water_diameters.csv'
    PathDict['C097','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C097','NevName']=''
    PathDict['C097','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C097/'
    PathDict['C097','FileName2DS']='C097_20180423_2DS_AllIn.h5'
    PathDict['C097','FileNameHVPS'] = 'C097_20180423_HVPS_AllIn.h5'    
      
    PathDict['C082','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    PathDict['C082','FAAMCoreName']='core_faam_20180214_v004_r0_c082_1hz.nc'
    PathDict['C082','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_processed/'
    PathDict['C082','CoreCloudFile']='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    PathDict['C082','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
    PathDict['C082','CdpCalFile']='C082_water_diameters.csv'
    PathDict['C082','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C082','NevName']='c082_nevzorov_20180214_1hz_r0.nc'
    PathDict['C082','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C082/'
    PathDict['C082','FileName2DS']='C082_20180214_2DS_AllIn.h5'
    PathDict['C082','FileNameHVPS'] = 'C082_20180214_HVPS.h5'
    
    PathDict['C081','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C081','FAAMCoreName']='core_faam_20180213_v004_r0_c081_1hz.nc'
    PathDict['C081','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C081','CoreCloudFile']='core-cloud-phy_faam_20180213_v501_r0_c081.nc'
    PathDict['C081','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C081','CdpCalFile']='c081_water_diameters.csv'
    PathDict['C081','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C081','NevName']='c081_nevzorov_20180213_1hz_r0.nc'
    PathDict['C081','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'
    PathDict['C081','FileName2DS']='C081_20180213_2DS_Allin.h5'
    PathDict['C081','FileNameHVPS'] = 'C081_20180213_HVPS_AllIn.h5'
   
    PathDict['C080','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C080','FAAMCoreName']='core_faam_20180208_v004_r0_c080_1hz.nc'
    PathDict['C080','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C080','CoreCloudFile']='core-cloud-phy_faam_20180208_v501_r0_c080.nc'
    PathDict['C080','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C080','CdpCalFile']='c080_water_diameters.csv'
    PathDict['C080','NevPath']=''
    PathDict['C080','NevName']=''
    PathDict['C080','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C080/'
    PathDict['C080','FileName2DS']='C080_20180208_2DS_AllIn.h5'
    PathDict['C080','FileNameHVPS'] = 'C080_20170208_HVPS_AllIn.h5'   
    
    PathDict['C079','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C079','FAAMCoreName']='core_faam_20180208_v004_r0_c079_1hz.nc'
    PathDict['C079','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C079','CoreCloudFile']='core-cloud-phy_faam_20180208_v501_r0_c079.nc'
    PathDict['C079','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C079','CdpCalFile']='c079_water_diameters.csv'
    PathDict['C079','NevPath']=''
    PathDict['C079','NevName']=''
    PathDict['C079','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C079/'
    PathDict['C079','FileName2DS']='C079_20180208_2DS_AllIn.h5'
    PathDict['C079','FileNameHVPS'] = 'C079_20180208_HVPS_AllIn'     
    
    PathDict['C078','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C078','FAAMCoreName']='core_faam_20180207_v004_r0_c078_1hz.nc'
    PathDict['C078','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C078','CoreCloudFile']='core-cloud-phy_faam_20180207_v501_r0_c078.nc'
    PathDict['C078','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C078','CdpCalFile']='c078_water_diameters.csv'
    PathDict['C078','NevPath']=''
    PathDict['C078','NevName']=''
    PathDict['C078','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C078/'
    PathDict['C078','FileName2DS']='C078_20180207_2DS_Allin.h5'
    PathDict['C078','FileNameHVPS'] ='' 
    
    PathDict['C077','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C077','FAAMCoreName']='core_faam_20180129_v004_r0_c077_1hz.nc'
    PathDict['C077','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C077','CoreCloudFile']='core-cloud-phy_faam_20180129_v501_r0_c077.nc'
    PathDict['C077','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C077','CdpCalFile']='c077_water_diameters.csv'
    PathDict['C077','NevPath']=''
    PathDict['C077','NevName']=''
    PathDict['C077','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C077/'
    PathDict['C077','FileName2DS']='C077_20180129_2DS_AllIn.h5'
    PathDict['C077','FileNameHVPS'] = '' 
    
    PathDict['C076','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C076','FAAMCoreName']='core_faam_20180124_v004_r0_c076_1hz.nc'
    PathDict['C076','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C076','CoreCloudFile']='core-cloud-phy_faam_20180124_v501_r0_c076.nc'
    PathDict['C076','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C076','CdpCalFile']='c076_water_diameters.csv'
    PathDict['C076','NevPath']=''
    PathDict['C076','NevName']=''
    PathDict['C076','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C076/'
    PathDict['C076','FileName2DS']='C076_20180124_2DS_AllIn.h5'
    PathDict['C076','FileNameHVPS'] = ''     
    
    
      
    return PathDict


#_________________________________________________________________________________________________

def PlotFlight(PathDict,Flight):
    #PlotOAP(PathDict[Flight,'FilePathOAP'],PathDict[Flight,'FileName2DS'])
    PlotOAP(PathDict[Flight,'FilePathOAP'],PathDict[Flight,'FileNameHVPS'])

#_________________________________________________________________________________________________

# Create dictionary with all 1hz data on core time base. Option to include psds but there on a different time base


def LoadFlightData2Dict(IncPSDs,FlightNumber,PathDict):
    
    #Paths 
    FullPath=PathDict[FlightNumber,'FullPath']
    FAAMCoreName=PathDict[FlightNumber,'FAAMCoreName']     
    CoreCloudPath=PathDict[FlightNumber,'CoreCloudPath']
    CoreCloudFile=PathDict[FlightNumber,'CoreCloudFile']
    CdpCalPath=PathDict[FlightNumber,'CdpCalPath']
    CdpCalFile=PathDict[FlightNumber,'CdpCalFile']
    NevPath=PathDict[FlightNumber,'NevPath']
    NevName=PathDict[FlightNumber,'NevName']
    FilePathOAP=PathDict[FlightNumber,'FilePathOAP']
    FileName2DS=PathDict[FlightNumber,'FileName2DS']
    FileNameHVPS = PathDict[FlightNumber,'FileNameHVPS']
    
    # Load core data
    FlightDict=loadFAAMCore(FullPath,FAAMCoreName)
    Time_Core=FlightDict['Time_Core']
    #Time_Core_Igor=DateTime2IgorTime(Time_Core)
    #TAT_DI_R=FAAMCoreDict['TAT_DI_R']
    #TAT_DI_R_C=TAT_DI_R-273.15
    #TAT_ND_R=FAAMCoreDict['TAT_ND_R']
    #TAT_ND_R_C=TAT_ND_R-273.15
    #ALT_GIN=FAAMCoreDict['ALT_GIN']
    #PS_RVSM=FAAMCoreDict['PS_RVSM']
    
    #Load CDP data
    CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    #CDP_time_mid_Igor=DateTime2IgorTime(CDP_time_mid)    
    # Change CDP to core timebase 
    #CDP_cm3_core=np.interp(Time_Core_Igor,CDP_time_mid_Igor,CDP_cm3)
    CDP_cm3_core=ChangeTimeBaseAvg(CDP_time_mid,CDP_cm3,Time_Core,1)
    #CDP_gcm3_core=np.interp(Time_Core_Igor,CDP_time_mid_Igor,CDP_gcm3)
    CDP_gcm3_core=ChangeTimeBaseAvg(CDP_time_mid,CDP_gcm3,Time_Core,1)
         
    CDP_gm3_core=CDP_gcm3_core*1e6
    FlightDict['CDP_cm3']= CDP_cm3_core
    FlightDict['CDP_gm3']= CDP_gm3_core
    
    if IncPSDs==1 :
        CDP_dNdDp_core=ChangeTimeBase2DAvg(CDP_time_mid,CDP_dNdDp,Time_Core, 1)        
        FlightDict['CDP_dNdDp']=CDP_dNdDp_core
        #FlightDict['CDP_time_mid']=CDP_time_mid
        FlightDict['CDP_size']=CDP_BinCentre
    
    #Create chilbolton warden hill radial flag
    RadialFlag=CreateRadialFlag(FlightDict)
    FlightDict['RadialFlag']= RadialFlag
    
    #Nevzorov
    TWC_g_m3, LWC_g_m3, TimeNev, TimeNevIgor= LoadNevzorov(NevPath,NevName)
    IWC_g_m3=TWC_g_m3-LWC_g_m3
    
    #TWC_g_m3_core=np.interp(Time_Core_Igor,TimeNevIgor,TWC_g_m3)
    TWC_g_m3_core=ChangeTimeBaseAvg(TimeNev,TWC_g_m3,Time_Core,1)
    #LWC_g_m3_core=np.interp(Time_Core_Igor,TimeNevIgor,LWC_g_m3)
    LWC_g_m3_core=ChangeTimeBaseAvg(TimeNev,LWC_g_m3,Time_Core,1)
    IWC_g_m3_core=ChangeTimeBaseAvg(TimeNev,IWC_g_m3,Time_Core,1)
    
    FlightDict['TWC_g_m3']= TWC_g_m3_core
    FlightDict['LWC_g_m3']= LWC_g_m3_core
    FlightDict['IWC_g_m3']= IWC_g_m3_core
    
    #2DS
    NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_LI_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime= LoadOAP(FilePathOAP,FileName2DS)    
    
    NC_HI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_HI_total,Time_Core,1)
    NC_MI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_MI_total,Time_Core,1)
    NC_LI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_LI_total,Time_Core,1)
    NC_All_accept_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_All_accept_total,Time_Core,1)
    
    FlightDict['NC_HI_2DS']= NC_HI_2DS
    FlightDict['NC_MI_2DS']= NC_MI_2DS   
    FlightDict['NC_LI_2DS']= NC_LI_2DS 
    FlightDict['NC_All_accept_2DS']= NC_All_accept_2DS   
    
    if IncPSDs==1 :
        NC_All_accept_CH0_core=ChangeTimeBase2DAvg(NC_DateTime,NC_All_accept_CH0_z,Time_Core, 1)       
        NC_HI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_HI_z,Time_Core, 1) 
        NC_MI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_MI_z,Time_Core, 1) 
        NC_LI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_LI_z,Time_Core, 1) 
        
        FlightDict['NC_All_accept_2DS_psd']=NC_All_accept_CH0_core
        FlightDict['NC_HI_z_2DS_psd']=NC_HI_z_core
        FlightDict['NC_MI_z_2DS_psd']=NC_MI_z_core
        FlightDict['NC_LI_z_2DS_psd']=NC_LI_z_core
        FlightDict['NC_2DS_size']=(NC_all_y[:-1:1]+NC_all_y[1::1])/2
    
    
    #HVPS
    NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_LI_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime= LoadOAP(FilePathOAP,FileNameHVPS)    
    NC_All_accept_total=np.sum(NC_All_accept_CH1_z,axis=1) # need to use channel 1 for HVPS
    #NC_DateTimeIgor=DateTime2IgorTime(NC_DateTime)
    #NC_HI_HVPS=np.interp(Time_Core_Igor,NC_DateTimeIgor,NC_HI_total)
    #NC_MI_HVPS=np.interp(Time_Core_Igor,NC_DateTimeIgor,NC_MI_total)
    #NC_LI_HVPS=np.interp(Time_Core_Igor,NC_DateTimeIgor,NC_LI_total)
    #NC_All_accept_HVPS=np.interp(Time_Core_Igor,NC_DateTimeIgor,NC_All_accept_total)
   
    NC_HI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_HI_total,Time_Core,1)
    NC_MI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_MI_total,Time_Core,1)
    NC_LI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_LI_total,Time_Core,1)
    NC_All_accept_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_All_accept_total,Time_Core,1)
    
    FlightDict['NC_HI_HVPS']= NC_HI_HVPS
    FlightDict['NC_MI_HVPS']= NC_MI_HVPS   
    FlightDict['NC_LI_HVPS']= NC_LI_HVPS 
    FlightDict['NC_All_accept_HVPS']= NC_All_accept_HVPS

    if IncPSDs==1 :
        NC_All_accept_CH0_core=ChangeTimeBase2DAvg(NC_DateTime,NC_All_accept_CH0_z,Time_Core, 1)       
        NC_HI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_HI_z,Time_Core, 1) 
        NC_MI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_MI_z,Time_Core, 1) 
        NC_LI_z_core=ChangeTimeBase2DAvg(NC_DateTime,NC_LI_z,Time_Core, 1) 
        
        FlightDict['NC_All_accept_HVPS_psd']=NC_All_accept_CH0_core
        FlightDict['NC_HI_z_HVPS_psd']=NC_HI_z_core
        FlightDict['NC_MI_z_HVPS_psd']=NC_MI_z_core
        FlightDict['NC_LI_z_HVPS_psd']=NC_LI_z_core
        FlightDict['NC_HVPS_size']=(NC_all_y[:-1:1]+NC_all_y[1::1])/2    
    
    return FlightDict
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
    
    NC_All_accept_total=np.sum(NC_All_accept_CH1_z,axis=1)
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

    plt.ylabel('All accept, L$^{-1}$')  # we already handled the x-label with ax1
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
      
    plt.gca().xaxis.set_major_locator(mdates.HourLocator() )
    plt.gca().xaxis.set_minor_locator(mdates.MinuteLocator(30) )

    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xlabel('time') 
    #plt.gca().set_xlim([round_time(dt=NC_DateTime[0], dateDelta=datetime.timedelta(minutes=15),to='down'),NC_DateTime[-1]])
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
    #Flights2Plot=['c075', 'c076', 'c077','c078','c079','c080','c081','c082']
    Flights2Plot=['c097','c098']
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
    plt.gca().xaxis.set_major_locator(mdates.HourLocator() )
    plt.gca().xaxis.set_minor_locator(mdates.MinuteLocator(30) )
    
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

#_____________________________________________________________________________________________________
    

#Plot flight summary

def PlotFlightSummary(FlightDict,FlightNumber,SavePath):
   
    Time_Core=FlightDict['Time_Core']
    ALT_GIN=FlightDict['ALT_GIN']
    TAT_DI_R_C=FlightDict['TAT_DI_R_C']
    TAT_ND_R_C=FlightDict['TAT_ND_R_C']
    CDP_cm3=FlightDict['CDP_cm3']
    CDP_gm3=FlightDict['CDP_gm3']
    IWC_g_m3=FlightDict['IWC_g_m3']
    LWC_g_m3=FlightDict['LWC_g_m3']
    NC_HI_2DS=FlightDict['NC_HI_2DS']
    NC_MI_2DS=FlightDict['NC_MI_2DS']
    NC_LI_2DS=FlightDict['NC_LI_2DS']
    NC_HI_HVPS=FlightDict['NC_HI_HVPS']
    NC_MI_HVPS=FlightDict['NC_MI_HVPS']
    NC_LI_HVPS=FlightDict['NC_LI_HVPS']      
 
    fig=plt.figure(figsize=(10,10)) 
    formatter = DateFormatter('%H:%M')
    ax1=fig.add_subplot(5, 1, 1)
    plt.title(FlightNumber)

    #plt.gca().set_ylim([0, 8000]) 
    plt.ylabel('Altitude, m',color='tab:gray')  # we already handled the x-label with ax1
    plt.plot(Time_Core, ALT_GIN, color='tab:gray', label='Altitude')
    #plt.tick_params(axis='245° radial', labelcolor='tab:red')
    ax2 = ax1.twinx()
    ax2.plot(Time_Core, TAT_ND_R_C, '-',color='tab:blue', label='T non deiced')
    ax2.plot(Time_Core, TAT_DI_R_C, '-',color='tab:red', label='T deiced')
    ax2.set_ylabel('Temperature, °C')#, color='tab:blue')
    
    
    #plt.ylabel(r'\textcolor{red}{Today} '+ r'\textcolor{green}{is} '+ r'\textcolor{blue}{cloudy.}')        
    #ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    plt.gca().set_ylim([-45, 10])
    
    ax3=plt.subplot(5, 1, 2, sharex=ax1)
    plt.plot(Time_Core, CDP_gm3, '-',color='tab:gray')
    plt.ylabel('CDP, g m$^{-3}$',color='tab:gray')
    plt.gca().set_ylim([0, 0.6])
    ax4=ax3.twinx()
    ax4.plot(Time_Core, CDP_cm3, '-',color='tab:red')
    plt.ylabel('CDP, cm$^{-3}$',color='tab:red')
    plt.gca().set_ylim([0, 120])
    
    ax5=plt.subplot(5, 1, 3, sharex=ax1)
    plt.plot(Time_Core, LWC_g_m3, '-',color='tab:gray')
    plt.ylabel('LWC Nevzorov, g m$^{-3}$',color='tab:gray')
    ax6=ax5.twinx()
    ax6.plot(Time_Core, IWC_g_m3, '-',color='tab:red')
    plt.ylabel('IWC Nevzorov, g m$^{-3}$',color='tab:red')
    
    ax7=plt.subplot(5, 1, 4, sharex=ax1)  
    plt.plot(Time_Core, NC_LI_2DS, '-',color='tab:gray', label='LI')
    plt.ylabel('2DS LI, L$^{-1}$',color='tab:gray')
    ax7.set_yscale('log')
    ax8=ax7.twinx()
    ax8.plot(Time_Core, NC_MI_2DS, '-',color='tab:red', label='MI')
    ax8.plot(Time_Core, NC_HI_2DS, '-',color='tab:blue', label='HI')
    ax8.legend(loc='upper right')
    ax8.set_yscale('log')
    plt.ylabel('2DS, L$^{-1}$')
    
    
    ax9=plt.subplot(5, 1, 5, sharex=ax1)  
    plt.plot(Time_Core, NC_LI_HVPS, '-',color='tab:gray', label='LI')
    plt.ylabel('HVPS LI, L$^{-1}$',color='tab:gray')
    ax9.set_yscale('log')
    ax10=ax9.twinx()
    ax10.plot(Time_Core, NC_MI_HVPS, '-',color='tab:red', label='MI')
    ax10.plot(Time_Core, NC_HI_HVPS, '-',color='tab:blue', label='HI')
    ax10.legend(loc='upper right')
    ax10.set_yscale('log')
    plt.ylabel('HVPS, L$^{-1}$')
    
    #plt.gca().set_ylim([0, 0.7])
    #plt.ylabel('CDP, g m$^{-3}$')
 
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    StartPt=FindFirstGreaterThan(ALT_GIN,150)
    #roundTime(dt=Time_Core[StartPt], roundTo=30*60)
    #print(Time_Core[StartPt])
    plt.gca().set_xlim([round_time(dt=Time_Core[StartPt], dateDelta=datetime.timedelta(minutes=15),to='down'),Time_Core[-1]]) #C082
    plt.xlabel('Time')
    plt.tight_layout()       
    #plt.show()
    Figurename='Summary_'+FlightNumber+'.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)        
    


#_____________________________________________________________________________________________________
# doesn't average PSDs



def AverageFlightDict(FlightDict,nPts):
    
    VariableNames=list(FlightDict.keys()) 
    
    #print(VariableNames)    
    FlightDictAvg={}
    for NameStr in VariableNames:
        Array=FlightDict[NameStr]
        if (NameStr=='Time_Core'):
            FlightDictAvg[NameStr]=Average_nPts_datetime(Array,nPts)
        else:    
            if NameStr.endswith('size'):
                FlightDictAvg[NameStr]=Array
            else:              
                if len(Array.shape)==1 : 
                    FlightDictAvg[NameStr]=Average_nPts(Array,nPts)    
                if len(Array.shape)==2 :    
                    FlightDictAvg[NameStr]=Average_nPts_2D(Array,nPts) 
                    
                    
            #if len(Array.shape)==1 and not NameStr.endswith('size'): #don't average psd and bin sizes
            #    FlightDictAvg[NameStr]=Average_nPts(Array,nPts)
            
            
                       
    return FlightDictAvg
    

#_____________________________________________________________________________________________________
#mass diameter relationship from brown and francis 1995

def BrownFrancis(PSD,MidSize):
    	
    MidMass=np.zeros(len(MidSize))
    for i in range(len(MidSize)):
        if(MidSize[i]<100):
            prefix = 4.82
            exponent = -13
            power = 3.0
        else:
            prefix = 7.38
            exponent = -11
            power = 1.9
        MidMass[i] = prefix * (MidSize[i]**power) * (10**exponent)
	
    PSD_mass= PSD * np.transpose(MidMass)      
    
    TotalMass = np.sum(PSD_mass,axis=1)
    
    return PSD_mass, TotalMass  

#_____________________________________________________________________________________________________


def CompareIWC(FlightDict,FlightNumber,SavePath):
    
    Time_Core=FlightDict['Time_Core']
    IWC_g_m3_Nev=FlightDict['IWC_g_m3'] # Nev data 
    NC_2DS_size=FlightDict['NC_2DS_size']
    NC_HI_2DS_psd=FlightDict['NC_HI_z_2DS_psd']
    NC_MI_2DS_psd=FlightDict['NC_MI_z_2DS_psd']
    NC_HVPS_size=FlightDict['NC_HVPS_size']
    NC_HI_HVPS_psd=FlightDict['NC_HI_z_HVPS_psd']
    NC_MI_HVPS_psd=FlightDict['NC_MI_z_HVPS_psd']    
        
    Mass_HI_2DS_psd, Mass_HI_2DS=BrownFrancis(NC_HI_2DS_psd,NC_2DS_size)
    Mass_HI_2DS_g_m3=Mass_HI_2DS *1000
    Mass_MI_2DS_psd, Mass_MI_2DS=BrownFrancis(NC_MI_2DS_psd,NC_2DS_size)
    Mass_HI_MI_2DS_g_m3=(Mass_HI_2DS+Mass_MI_2DS) *1000   
    
    Mass_HI_HVPS_psd, Mass_HI_HVPS=BrownFrancis(NC_HI_HVPS_psd,NC_HVPS_size)
    Mass_HI_HVPS_g_m3=Mass_HI_HVPS *1000
    Mass_MI_HVPS_psd, Mass_MI_HVPS=BrownFrancis(NC_MI_HVPS_psd,NC_2DS_size)
    Mass_HI_MI_HVPS_2DS_g_m3=(Mass_HI_2DS+Mass_MI_2DS+Mass_HI_HVPS+Mass_MI_HVPS) *1000
    
    fig=plt.figure(figsize=(10,5)) 
    formatter = DateFormatter('%H:%M')
    
    plt.subplot(1, 2, 1)
    plt.title(FlightNumber)
    
    plt.plot(Time_Core,IWC_g_m3_Nev,label='Nevzorov')
    plt.plot(Time_Core,Mass_HI_MI_2DS_g_m3,label='2DS HI + MI (Brown & Francis)')
    #plt.plot(Time_Core,Mass_HI_2DS_g_m3,label='2DS HI (Brown & Francis)')
    #plt.plot(Time_Core,Mass_HI_HVPS_g_m3,label='HVPS HI (Brown & Francis)')
    plt.plot(Time_Core,Mass_HI_MI_HVPS_2DS_g_m3,label='HVPS + 2DS + HI + MI (Brown & Francis)')
    plt.xlabel('Time')
    plt.ylabel('IWC, g m$^{-3}$')
    plt.ylim([0, 0.3])

    plt.legend()
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    
    plt.subplot(1, 2, 2)
    plt.plot(IWC_g_m3_Nev, Mass_HI_MI_HVPS_2DS_g_m3,'o')
    plt.ylabel('Nevzorov, g m$^{-3}$')
    plt.xlabel('2DS + HVPS, g m$^{-3}$')
    
    plt.ylim([0, 0.3])
    plt.xlim([0, 0.3])
    
    Figurename='IWC_'+FlightNumber+'.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)   


    