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
from MyFunctions import Matlab2PythonTimeArray,loadmat,Matlab2PythonTime,BinLinearV2,LoadNevzorov,LoadOAP,round_time, Average_nPts_datetime,Average_nPts,Average_nPts_2D, haversine,BinMid_2_Width,BinAverage2D
from MyFunctions import LoadOAP_nc
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
from mpl_toolkits.basemap import Basemap
import math
from netCDF4 import Dataset
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import scipy.io as sio
from scipy import optimize

#_________________________________________________________________________________________________

# Create dictionary with all 1hz data on core time base. Option to include psds but there on a different time base


def DataPaths():
    PathDict={}
   
    PathDict['C098','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C098','FAAMCoreName']='core_faam_20180424_v004_r0_c098_1hz.nc'
    PathDict['C098','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C098','CoreCloudFile']='core-cloud-phy_faam_20180424_v501_r0_c098.nc'
    PathDict['C098','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C098','CdpCalFile']='c098_water_diameters.csv'
    PathDict['C098','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C098','NevName']=''
    PathDict['C098','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/'
    #PathDict['C098','FileName2DS']='C098_20180424_2DS_Allin.h5'
    PathDict['C098','FileName2DS']='UMAN_2DS_20180424_r0_C098.nc'
    #PathDict['C098','FileNameHVPS'] = 'C098_20180424_HVPS_AllIn.h5' 
    PathDict['C098','FileNameHVPS'] = 'UMAN_HVPS_20180424_r0_C098.nc'
    
    
    
    
    PathDict['C097','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C097','FAAMCoreName']='core_faam_20180423_v004_r0_c097_1hz.nc'
    PathDict['C097','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C097','CoreCloudFile']='core-cloud-phy_faam_20180423_v501_r0_c097.nc'
    PathDict['C097','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C097','CdpCalFile']='c097_water_diameters.csv'
    PathDict['C097','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C097','NevName']=''
    PathDict['C097','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C097/'
    #PathDict['C097','FileName2DS']='C097_20180423_2DS_AllIn.h5'
    PathDict['C097','FileName2DS']='UMAN_2DS_20180423_r0_C097.nc'
    #PathDict['C097','FileNameHVPS'] = 'C097_20180423_HVPS_AllIn.h5' 
    PathDict['C097','FileNameHVPS'] = 'UMAN_HVPS_20180423_r0_C097.nc'
      
    PathDict['C082','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C082','FAAMCoreName']='core_faam_20180214_v004_r0_c082_1hz.nc'
    PathDict['C082','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C082','CoreCloudFile']='core-cloud-phy_faam_20180214_v501_r0_c082.nc'
    PathDict['C082','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c082-feb-14/core_raw/core-cloud-phy_faam_20180214_r0_c082_rawpads/CDP_cal_processed/'
    PathDict['C082','CdpCalFile']='C082_water_diameters.csv'
    PathDict['C082','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C082','NevName']='c082_nevzorov_20180214_1hz_r0.nc'
    PathDict['C082','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C082/'
    #PathDict['C082','FileName2DS']='C082_20180214_2DS_AllIn.h5'
    PathDict['C082','FileName2DS']='UMAN_2DS_20180214_r0_C082.nc'
    #PathDict['C082','FileNameHVPS'] = 'C082_20180214_HVPS.h5'
    PathDict['C082','FileNameHVPS'] = 'UMAN_HVPS_20180214_r0_C082.nc'
      
    PathDict['C081','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C081','FAAMCoreName']='core_faam_20180213_v004_r0_c081_1hz.nc'
    PathDict['C081','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C081','CoreCloudFile']='core-cloud-phy_faam_20180213_v501_r0_c081.nc'
    PathDict['C081','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C081','CdpCalFile']='c081_water_diameters.csv'
    PathDict['C081','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
    PathDict['C081','NevName']='c081_nevzorov_20180213_1hz_r0.nc'
    PathDict['C081','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'
    #PathDict['C081','FileName2DS']='C081_20180213_2DS_Allin.h5'
    PathDict['C081','FileName2DS']='UMAN_2DS_20180213_r0_C081.nc'
    #PathDict['C081','FileNameHVPS'] = 'C081_20180213_HVPS_AllIn.h5'
    PathDict['C081','FileNameHVPS'] = 'UMAN_HVPS_20180213_r0_C081.nc'
    PathDict['C081','SEA_path']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/SEA/SEAmerge_c081/'
    PathDict['C081','SEA_filename']='seaprobe_core_processed_20180213_c081_MO-Nev.csv'
    
    
    
    
    PathDict['C080','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C080','FAAMCoreName']='core_faam_20180208_v004_r0_c080_1hz.nc'
    PathDict['C080','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C080','CoreCloudFile']='core-cloud-phy_faam_20180208_v501_r0_c080.nc'
    PathDict['C080','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C080','CdpCalFile']='c080_water_diameters.csv'
    PathDict['C080','NevPath']=''
    PathDict['C080','NevName']=''
    PathDict['C080','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C080/'
    #PathDict['C080','FileName2DS']='C080_20180208_2DS_AllIn.h5'
    PathDict['C080','FileName2DS']='UMAN_2DS_20180208_r0_C080.nc'
    #PathDict['C080','FileNameHVPS'] = 'C080_20180208_HVPS_AllIn.h5'   
    PathDict['C080','FileNameHVPS'] = 'UMAN_HVPS_20180208_r0_C080.nc'
    
    PathDict['C079','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C079','FAAMCoreName']='core_faam_20180208_v004_r0_c079_1hz.nc'
    PathDict['C079','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C079','CoreCloudFile']='core-cloud-phy_faam_20180208_v501_r0_c079.nc'
    PathDict['C079','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C079','CdpCalFile']='c079_water_diameters.csv'
    PathDict['C079','NevPath']=''
    PathDict['C079','NevName']=''
    PathDict['C079','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C079/'
    PathDict['C079','FileName2DS']='UMAN_2DS_20180208_r0_C079.nc'
    #PathDict['C079','FileNameHVPS'] = 'C079_20180208_HVPS_Allin.h5'     
    PathDict['C079','FileNameHVPS'] = 'UMAN_HVPS_20180208_r0_C079.nc'   
    
    PathDict['C078','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C078','FAAMCoreName']='core_faam_20180207_v004_r0_c078_1hz.nc'
    PathDict['C078','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C078','CoreCloudFile']='core-cloud-phy_faam_20180207_v501_r0_c078.nc'
    PathDict['C078','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C078','CdpCalFile']='c078_water_diameters.csv'
    PathDict['C078','NevPath']=''
    PathDict['C078','NevName']=''
    PathDict['C078','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C078/'
    PathDict['C078','FileName2DS']='UMAN_2DS_20180207_r0_C078.nc'
    PathDict['C078','FileNameHVPS'] ='' 
    
    PathDict['C077','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C077','FAAMCoreName']='core_faam_20180129_v004_r0_c077_1hz.nc'
    PathDict['C077','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C077','CoreCloudFile']='core-cloud-phy_faam_20180129_v501_r0_c077.nc'
    PathDict['C077','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C077','CdpCalFile']='c077_water_diameters.csv'
    PathDict['C077','NevPath']=''
    PathDict['C077','NevName']=''
    PathDict['C077','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C077/'
    PathDict['C077','FileName2DS']='UMAN_2DS_20180129_r0_C077.nc'
    PathDict['C077','FileNameHVPS'] = '' 
    
    PathDict['C076','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C076','FAAMCoreName']='core_faam_20180124_v004_r0_c076_1hz.nc'
    PathDict['C076','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C076','CoreCloudFile']='core-cloud-phy_faam_20180124_v501_r0_c076.nc'
    PathDict['C076','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C076','CdpCalFile']='c076_water_diameters.csv'
    PathDict['C076','NevPath']=''
    PathDict['C076','NevName']=''
    PathDict['C076','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C076/'
    PathDict['C076','FileName2DS']='UMAN_2DS_20180124_r0_C076.nc'
    PathDict['C076','FileNameHVPS'] = ''

    PathDict['C075','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C075','FAAMCoreName']='core_faam_20180123_v004_r0_c075_1hz.nc'
    PathDict['C075','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C075','CoreCloudFile']='core-cloud-phy_faam_20180123_v501_r0_c075.nc'
    PathDict['C075','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C075','CdpCalFile']='c075_water_diameters.csv'
    PathDict['C075','NevPath']=''
    PathDict['C075','NevName']=''
    PathDict['C075','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C075/'
    PathDict['C075','FileName2DS']='UMAN_2DS_20180123_r0_C075.nc'
    PathDict['C075','FileNameHVPS'] = ''
     
    PathDict['C074','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C074','FAAMCoreName']='core_faam_20180118_v004_r0_c074_1hz.nc'
    PathDict['C074','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C074','CoreCloudFile']='core-cloud-phy_faam_20180118_v501_r0_c074.nc'
    PathDict['C074','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C074','CdpCalFile']='c075_water_diameters.csv'
    PathDict['C074','NevPath']=''
    PathDict['C074','NevName']=''
    PathDict['C074','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C074/'
    PathDict['C074','FileName2DS']='UMAN_2DS_20180118_r0_C074.nc'
    PathDict['C074','FileNameHVPS'] = ''    
    
    PathDict['C073','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C073','FAAMCoreName']='core_faam_20180109_v004_r0_c073_1hz.nc'
    PathDict['C073','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C073','CoreCloudFile']='core-cloud-phy_faam_20180109_v501_r0_c073.nc'
    PathDict['C073','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C073','CdpCalFile']='c075_water_diameters.csv'
    PathDict['C073','NevPath']=''
    PathDict['C073','NevName']=''
    PathDict['C073','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C073/'
    PathDict['C073','FileName2DS']='C073_20180109_2DS_Allin.h5'
    PathDict['C073','FileNameHVPS'] = ''  
    
    PathDict['C072','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
    PathDict['C072','FAAMCoreName']='core_faam_20171214_v004_r0_c072_1hz.nc'
    PathDict['C072','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
    PathDict['C072','CoreCloudFile']='core-cloud-phy_faam_20171214_v501_r0_c072.nc'
    PathDict['C072','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
    PathDict['C072','CdpCalFile']='c075_water_diameters.csv'
    PathDict['C072','NevPath']=''
    PathDict['C072','NevName']=''
    PathDict['C072','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C072/'
    PathDict['C072','FileName2DS']='C072_20171214_2DS_ALLin.h5'
    PathDict['C072','FileNameHVPS'] = '' 
    
    
      
    return PathDict


#_________________________________________________________________________________________________

# loops through 



def PlotFlight(PathDict,FlightList):

    for i in range(len(FlightList)):    
        Flight=FlightList[i]
        #PlotOAP(PathDict[Flight,'FilePathOAP'],PathDict[Flight,'FileName2DS'])
        #LoadOAP_nc('2DS',PathDict[Flight,'FilePathOAP'],PathDict[Flight,'FileName2DS'],1)
        LoadOAP_nc('HVPS',PathDict[Flight,'FilePathOAP'],PathDict[Flight,'FileNameHVPS'],1)
        #FlightDict=LoadFlightData2Dict(0,0,0,Flight,PathDict)
        #TrackPlot(FlightDict,PathDict[Flight,'FullPath'],Flight)
        #LoadCoreCDP(PathDict[Flight,'FullPath'],PathDict[Flight,'FAAMCoreName'],PathDict[Flight,'CoreCloudPath'],PathDict[Flight,'CoreCloudFile'],PathDict[Flight,'CdpCalPath'],PathDict[Flight,'CdpCalFile'])
        #PlotCrossSections2CHB(FlightDict,PathDict[Flight,'FullPath'],Flight)

#_________________________________________________________________________________________________

# Create dictionary with all 1hz data on core time base. Option to include psds but there on a different time base


def LoadFlightData2Dict(IncSEA,IncHVPS,IncNev,IncPSDs,FlightNumber,PathDict):
    
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
    #RadialFlag=CreateRadialFlag(FlightDict)
    #FlightDict['RadialFlag']= RadialFlag
    
    #Nevzorov
    if IncNev==1 :
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
    #NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_LI_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime= LoadOAP(FilePathOAP,FileName2DS)    
    NC_DateTime, Time_mid, Time_edge, Size_mid, Size_edge, PSD_Num_S, PSD_Num_LI, PSD_Num_MI, PSD_Num_HI, PSD_Num_All =LoadOAP_nc('2DS',FilePathOAP,FileName2DS,0)
    
    NC_All_total=np.nansum(PSD_Num_All,axis=1)
    NC_HI_total=np.nansum(PSD_Num_HI,axis=1)
    NC_MI_total=np.nansum(PSD_Num_MI,axis=1)
    NC_LI_total=np.nansum(PSD_Num_LI,axis=1)
    
    NC_HI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_HI_total,Time_Core,1)
    NC_MI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_MI_total,Time_Core,1)
    NC_LI_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_LI_total,Time_Core,1)
    NC_All_accept_2DS=ChangeTimeBaseAvg(NC_DateTime,NC_All_total,Time_Core,1)
    
    FlightDict['NC_HI_2DS']= NC_HI_2DS
    FlightDict['NC_MI_2DS']= NC_MI_2DS   
    FlightDict['NC_LI_2DS']= NC_LI_2DS 
    FlightDict['NC_All_2DS']= NC_All_accept_2DS   
    
    if IncPSDs==1 :
        NC_All_accept_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_All,Time_Core, 1)       
        NC_HI_z_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_HI,Time_Core, 1) 
        NC_MI_z_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_MI,Time_Core, 1) 
        NC_LI_z_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_LI,Time_Core, 1) 
        
        FlightDict['PSD_Num_All_2DS']=NC_All_accept_core
        FlightDict['PSD_Num_HI_2DS']=NC_HI_z_core
        FlightDict['PSD_Num_MI_2DS']=NC_MI_z_core
        FlightDict['PSD_Num_LI_2DS']=NC_LI_z_core
        #FlightDict['NC_2DS_size']=(NC_all_y[:-1:1]+NC_all_y[1::1])/2 # Bin centre
        FlightDict['Size_mid_2DS']=Size_mid
        FlightDict['Size_edge_2DS']=Size_edge
    #HVPS
    
    if IncHVPS == 1:
        #NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_LI_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime= LoadOAP(FilePathOAP,FileNameHVPS)    
        #NC_All_accept_total=np.sum(NC_All_accept_CH1_z,axis=1) # need to use channel 1 for HVPS
        
        #NC_DateTime, Time_mid, Time_edge, Size_mid, Size_edge, PSD_Num_S, PSD_Num_LI, PSD_Num_MI, PSD_Num_HI, PSD_Num_All=LoadOAP_nc(FilePathOAP,FileNameHVPS,0)
        
        NC_DateTime, Time_mid, Time_edge, Size_mid, Size_edge, PSD_Num_S, PSD_Num_LI, PSD_Num_MI, PSD_Num_HI, PSD_Num_All =LoadOAP_nc('HVPS',FilePathOAP,FileNameHVPS,0)
        
        NC_All=np.sum(PSD_Num_All,axis=1)
        NC_HI=np.sum(PSD_Num_LI,axis=1)
        NC_MI=np.sum(PSD_Num_MI,axis=1)
        NC_LI=np.sum(PSD_Num_HI,axis=1)
        
        NC_HI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_HI,Time_Core,1)
        NC_MI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_MI,Time_Core,1)
        NC_LI_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_LI,Time_Core,1)
        NC_All_HVPS=ChangeTimeBaseAvg(NC_DateTime,NC_All,Time_Core,1)
        FlightDict['NC_HI_HVPS']= NC_HI_HVPS
        FlightDict['NC_MI_HVPS']= NC_MI_HVPS   
        FlightDict['NC_LI_HVPS']= NC_LI_HVPS 
        FlightDict['NC_All_HVPS']= NC_All_HVPS

        if IncPSDs==1 :
            PSD_Num_All_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_All,Time_Core, 1)       
            PSD_Num_HI_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_HI,Time_Core, 1) 
            PSD_Num_MI_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_MI,Time_Core, 1) 
            PSD_Num_LI_core=ChangeTimeBase2DAvg(NC_DateTime,PSD_Num_LI,Time_Core, 1) 
        
            FlightDict['PSD_Num_All_HVPS']=PSD_Num_All_core
            FlightDict['PSD_Num_HI_HVPS']=PSD_Num_HI_core
            FlightDict['PSD_Num_MI_HVPS']=PSD_Num_MI_core
            FlightDict['PSD_Num_LI_HVPS']=PSD_Num_LI_core
            FlightDict['Size_mid_HVPS']=Size_mid
            FlightDict['Size_edge_HVPS']=Size_edge
    
    if IncSEA == 1 : 
        SEA_path = PathDict[FlightNumber,'SEA_path']
        SEA_filename = PathDict[FlightNumber,'SEA_filename']
        df = pd.read_csv(SEA_path+SEA_filename)
        TimeStartStr=df['TIMESTAMP']
        SEA_Time=[]
        #SEA_Time=[parser.parse(TimeStartStr[x]) for x in range(len(TimeStartStr))]
        SEA_Time= [datetime.datetime.strptime(TimeStartStr[x], "%Y-%m-%d %H:%M:%S.%f") for x in range(len(TimeStartStr))]
        LWC_021_core=ChangeTimeBaseAvg(SEA_Time,df['021_LWC'],Time_Core,1) 
        IWC_021_core=ChangeTimeBaseAvg(SEA_Time,df['021_IWC'],Time_Core,1)
        TWC_021_core=ChangeTimeBaseAvg(SEA_Time,df['021_TWC'],Time_Core,1)
        LWC_083_core=ChangeTimeBaseAvg(SEA_Time,df['083_LWC'],Time_Core,1)
        IWC_083_core=ChangeTimeBaseAvg(SEA_Time,df['083_IWC'],Time_Core,1)
        TWC_083_core=ChangeTimeBaseAvg(SEA_Time,df['083_TWC'],Time_Core,1)
        #021_LWC=ChangeTimeBaseAvg(SEA_Time,021_LWC,Time_Core,1)            
        FlightDict['LWC_021']=LWC_021_core
        FlightDict['IWC_021']=IWC_021_core
        FlightDict['TWC_021']=TWC_021_core
        FlightDict['LWC_083']=LWC_083_core
        FlightDict['IWC_083']=IWC_083_core
        FlightDict['TWC_083']=TWC_083_core
    
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

#    PathDict['C098','FullPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCore/'
#    PathDict['C098','FAAMCoreName']='core_faam_20180424_v004_r0_c098_1hz.nc'
#    PathDict['C098','CoreCloudPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCoreCloud/'
#    PathDict['C098','CoreCloudFile']='core-cloud-phy_faam_20180424_v501_r0_c098.nc'
#    PathDict['C098','CdpCalPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/AllCDPcals/'
#    PathDict['C098','CdpCalFile']='c098_water_diameters.csv'
#    PathDict['C098','NevPath']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/Nevzorov/'
#    PathDict['C098','NevName']=''
#    PathDict['C098','FilePathOAP']='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/'
#    PathDict['C098','FileName2DS']='C098_20180424_2DS_Allin.h5'
#    PathDict['C098','FileNameHVPS'] = 'C098_20180424_HVPS_AllIn.h5'

#FlightStr='C081'

def PSD_runs_dict(PathDict, FlightStr,CSVPath,CSVName):
    
    FilePath= PathDict[FlightStr,'FilePathOAP']
    FlightDict=LoadFlightData2Dict(1,0,1,FlightStr,PathDict)
    #CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
    #CSVName='flight-sum_faam_20180213_r0_c081.csv'
    StartTime, EndTime, RunNumber,ProfileNumber=LoadFlightSummaryCSV(CSVPath,CSVName)

    for i in range(len(StartTime)):
        if (RunNumber[i] > 0) : 
            PSD_comparison_Dict(FlightDict, StartTime[i],EndTime[i],FilePath,'PSD_Num_All_run_'+str(int(RunNumber[i])),'PSD_Num_All')            


#_________________________________________________________________________________________________


def LoadFlightSummaryCSV(FilePath,FileName):
    df = pd.read_csv(FilePath+FileName)
    DateStr=df['Date']
    TimeStartStr=df['TimeStart']
    TimeStopStr=df['TimeStop']
    RunNumber=df['RunNumber']
    ProfileNumber=df['ProfileNumber']
    
    StartTime=[]#*len(DateStr)
    EndTime=[]
    
    for i in range(len(DateStr)):
        StartTime.append(datetime.datetime.strptime(DateStr[i]+' '+TimeStartStr[i], "%Y/%m/%d %H:%M:%S"))
        EndTime.append(datetime.datetime.strptime(DateStr[i]+' '+TimeStopStr[i], "%Y/%m/%d %H:%M:%S"))
    return StartTime, EndTime, RunNumber,ProfileNumber    

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

# Plot a psd for the probes using data from FlightDict


def PSD_comparison_Dict(FlightDict, StartTime,EndTime,FilePath,Figurename,Prefix):            

    
    
    #2DS
    NC_SizeMid_2DS, NC_z_dNdDp_2DS= Avg_OAP_dndDp_Dict(FlightDict, 'Size_mid_2DS', 'Time_Core', Prefix+'_2DS', StartTime,EndTime,1)   
    
    #HVPS
    NC_SizeMid_HVPS, NC_z_dNdDp_HVPS= Avg_OAP_dndDp_Dict(FlightDict, 'Size_mid_HVPS', 'Time_Core', Prefix+'_HVPS', StartTime,EndTime,1)
    
    #print(np.nanmax(NC_z_dNdDp_HVPS))
       
    #CDP
    CDP_BinCentre=FlightDict['CDP_size']
    SizeMid, CDP_dNdDp_L_avg=Avg_CDP_dndDp_Dict(FlightDict,StartTime,EndTime)
    
    #Temperature
    Time_Core=FlightDict['Time_Core']
    TAT_DI_R=FlightDict['TAT_DI_R']
    TAT_ND_R=FlightDict['TAT_ND_R']
    StartIdx=bisect.bisect_left(Time_Core, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Core, EndTime) # assume time is sorted, which it should be
    TAT_DI_R_avg=np.around(np.nanmean(TAT_DI_R[StartIdx:EndIdx])-273.15, decimals=1)
    TAT_ND_R_avg=np.around(np.nanmean(TAT_ND_R[StartIdx:EndIdx])-273.15, decimals=1)
    
    #Nevzorov
    #IWC_gm3=FlightDict['IWC_g_m3']    
    #IWC_avg=np.around(np.nanmean(IWC_gm3[StartIdx:EndIdx]), decimals=3)
    
    #Plot    
    fig=plt.figure(figsize=(10,10)) 
    fig, ax = plt.subplots()
    #plt.title('Temperature DI='+str(TAT_DI_R_avg)+'°C, ND='+ str(TAT_ND_R_avg)+'°C, IWC='+str(IWC_avg)+' g m$^{-3}$ ')
    plt.title('Temperature DI='+str(TAT_DI_R_avg)+'°C, ND='+ str(TAT_ND_R_avg)+'°C')
    plt.plot(NC_SizeMid_2DS,NC_z_dNdDp_2DS, label='2DS')
    plt.plot(NC_SizeMid_HVPS,NC_z_dNdDp_HVPS, label='HVPS')
    
    #if ProbeFlag==1 :
    #    plt.plot(NC_SizeMid_CIP100,NC_z_dNdDp_CIP100, label='CIP100')
    #    plt.plot(NC_SizeMid_CIP15,NC_z_dNdDp_CIP15, label='CIP15')
    #    plt.plot(NC_SizeMid_3VCHO,NC_z_dNdDp_3VCH0, label='3VCH0')
    #    plt.plot(NC_SizeMid_3VCH1,NC_z_dNdDp_3VCH1, label='3VCH1')
        
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
    
   
    ProbeFlag==0
    
    if ProbeFlag==1 :
    
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
   
def Avg_CDP_dndDp_Dict(FlightDict,StartTime,EndTime):
   

    Time_Core=np.array(FlightDict['Time_Core'])
    SizeMid=np.array(FlightDict['CDP_size']) # Bin Midpoint
    PSD=np.array(FlightDict['CDP_dNdDp']) 
    
    
    #Search for start and stop times
    StartIdx=bisect.bisect_left(Time_Core, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Core, EndTime) # assume time is sorted, which it should be

    CDP_dNdDp_L_avg=np.nanmean(PSD[StartIdx:EndIdx][:],axis=0)
    
    CDP_dNdDp_L_avg*=1000 # Convert to L
    
    return SizeMid, CDP_dNdDp_L_avg
    




#_________________________________________________________________________________________________

def Avg_dndDp(StartTime,EndTime,PSD, Time_Mid, Size_Mid):
    
    
    #Search for start and stop times
    StartIdx=bisect.bisect_left(Time_Mid, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Mid, EndTime) # assume time is sorted, which it should be 

    PSD_avg=np.nanmean(PSD[StartIdx:EndIdx][:],axis=0)

    #PSD_avg=np.nansum(PSD[StartIdx:EndIdx][:],axis=0)
    
    SizeWidth=BinMid_2_Width(Size_Mid)
    
    PSD_avg_dNdDp=PSD_avg/SizeWidth
    
    return SizeWidth, PSD_avg_dNdDp, PSD_avg

#_________________________________________________________________________________________________
#Average PSD between StartTime and EndTime

# PSD_avg_CountStats, PSD_avg_dNdDp_CountStats = mean / sqrt(N)

def Avg_dndDp_CountStats(StartTime,EndTime,PSD,Counts_PSD, Time_Mid, Size_Mid):
    
    
    #Search for start and stop times
    StartIdx=bisect.bisect_left(Time_Mid, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Mid, EndTime) # assume time is sorted, which it should be 

    PSD_avg = np.nanmean(PSD[StartIdx:EndIdx][:],axis=0)
    Counts = np.nansum(Counts_PSD[StartIdx:EndIdx][:],axis=0)
    
    #PSD_avg=np.nansum(PSD[StartIdx:EndIdx][:],axis=0)
    
    SizeWidth=BinMid_2_Width(Size_Mid)
    
    PSD_avg_dNdDp=PSD_avg/SizeWidth
    
    #PSD_avg_dNdDp_CountStats = PSD_avg_dNdDp / np.sqrt(Counts_avg)
    #PSD_avg_CountStats = PSD_avg / np.sqrt(Counts_avg)
    
    return SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts


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

# Flag should be same length as TimeOAP. Flag 1 = include 



def Avg_OAP_dndDp_Dict(FlightDict, SizeName, TimeName, dN_Name, StartTime,EndTime,Flag):            

    TimeOAP=np.array(FlightDict[TimeName])
    SizeMid=np.array(FlightDict[SizeName]) # Bin Midpoint
    PSD=np.array(FlightDict[dN_Name]) # dN
        
    if not isinstance(Flag,int) :
        PSD[Flag!=1,:]=np.nan
        
    #Search for start and stop times
    StartIdx=bisect.bisect_left(TimeOAP, StartTime) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(TimeOAP, EndTime) # assume time is sorted, which it should be 

    PSD_mean=np.nanmean(PSD[StartIdx:EndIdx][:],axis=0)

    BinWidth=BinMid_2_Width(SizeMid)

    dNdDp=PSD_mean/BinWidth 
    
    return SizeMid, dNdDp




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

#Plot OAP timeseries


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
    
    #MVD=CalculateMVD(PSD,NC_all_x)
    
    
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
    Flights2Plot=['c075', 'c076', 'c077','c078','c079','c080','c081','c082','c097','c098']
    #Flights2Plot=['c097','c098']
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

# plot some core data and cdp


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
    Figurename=CoreCloudPath+FAAMCoreName[27:31]+'_CDP_'+FAAMCoreName[10:26]+'.png'
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

# Creates a flag which is approx when on the chilbolton warden hill line

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
    plt.gca().set_ylim([0,0.5])
    ax6=ax5.twinx()
    ax6.plot(Time_Core, IWC_g_m3, '-',color='tab:red')
    plt.gca().set_ylim([0, 0.5])
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
# average every nPts of flight dictionary 



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

# Compare IWC from nev and 2ds+HVPS


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
    plt.rcParams.update({'font.size': 10})
    formatter = DateFormatter('%H:%M')
    
    plt.subplot(1, 1, 1)
    plt.title(FlightNumber)
    
    plt.plot(Time_Core,IWC_g_m3_Nev,label='Nevzorov')
    plt.plot(Time_Core,Mass_HI_MI_2DS_g_m3,label='2DS HI + MI (Brown & Francis)')
    #plt.plot(Time_Core,Mass_HI_2DS_g_m3,label='2DS HI (Brown & Francis)')
    #plt.plot(Time_Core,Mass_HI_HVPS_g_m3,label='HVPS HI (Brown & Francis)')
    plt.plot(Time_Core,Mass_HI_MI_HVPS_2DS_g_m3,label='HVPS + 2DS + HI + MI (Brown & Francis)')
    plt.xlabel('Time')
    plt.ylabel('IWC, g m$^{-3}$')
    plt.ylim([0, 0.5])

    plt.legend()
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    
    #plt.subplot(1, 2, 2)
    #plt.plot(IWC_g_m3_Nev, Mass_HI_MI_HVPS_2DS_g_m3,'o')
    #plt.ylabel('Nevzorov, g m$^{-3}$')
    #plt.xlabel('2DS + HVPS, g m$^{-3}$')
    
    #plt.ylim([0, 0.3])
    #plt.xlim([0, 0.3])
    
    Figurename='IWC_'+FlightNumber+'.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)   


#_____________________________________________________________________________________________________

#Plot track plot over map


def TrackPlot(FlightDict,SavePath,FlightNumber):
    LON=FlightDict['LON_GIN']
    LAT=FlightDict['LAT_GIN']
    
    
    fig=plt.figure(figsize=(10,5))
    plt.rcParams.update({'font.size': 16})
    #ax=fig.add_axes([0.1,0.1,0.8,0.8])
    # setup mercator map projection.
    #map = Basemap(llcrnrlon=-8.,llcrnrlat=49.,urcrnrlon=2.,urcrnrlat=53.,resolution='l',projection='merc',lat_0=40.,lon_0=-20.,lat_ts=20.)
    map = Basemap(llcrnrlon=-7.,llcrnrlat=49.5,urcrnrlon=2.,urcrnrlat=53.,resolution='l',projection='merc',lat_0=40.,lon_0=-20.,lat_ts=20., epsg=5520)

    

    #map.drawcoastlines()
    #map.fillcontinents()
    x,y=map(LON,LAT)
    map.plot(x, y, 'o',markersize=0.5,linewidth=0, color='tab:red')
    
    Landmark_Lat=[51.145037, 50.819167, 52.069791]
    Landmark_Lon=[ -1.438416, -2.556333, -0.619135]
    Labels = ['Chilbolton', 'Wardon Hill', 'Cranfield']



    x2,y2=map(Landmark_Lon,Landmark_Lat)
    map.plot(x2, y2, '+',markersize=10,linewidth=0, color='tab:green')
    
    for i in range(len(Labels)):
        plt.annotate(
            Labels[i],
            xy=(x2[i], y2[i]), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    
    
    map.drawparallels(np.arange(10,90,1),labels=[1,1,0,0])
    # draw meridians
    map.arcgisimage(service='ESRI_Imagery_World_2D', xpixels = 1500, verbose= True)
    map.drawmeridians(np.arange(-180,180,1),labels=[0,1,0,1])
    #ax.set_title('Great Circle from New York to London')
    Figurename=FlightNumber+'FlightPath.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)
    #plt.show()
   
#_____________________________________________________________________________________________________

# Plot lat lon cross sections


def PlotCrossSections(FlightDict,SavePath,FlightNumber):
    
    NC_HI_2DS=FlightDict['NC_HI_2DS'] 
    CDP_cm3=FlightDict['CDP_cm3']
    LON=FlightDict['LON_GIN']
    LAT=FlightDict['LAT_GIN']    
    ALT=FlightDict['ALT_GIN']
    
    ALT[ (LAT > 51.5) | (LON > -1.0) ] =np.nan
    NC_HI_2DS[ (LAT > 51.5) | (LON > -1.0) ] =np.nan
    CDP_cm3[ (LAT > 51.5) | (LON > -1.0) ] =np.nan
    
    fig=plt.figure(figsize=(10,10)) 
    
    plt.title(FlightNumber)
    ax1=plt.subplot(2, 1, 1)
    plt.title(FlightNumber)
    cm = plt.cm.get_cmap('jet')
    plt.plot(LON, ALT,color='gray',markersize=0,zorder=0)
    CDP_cm3[(CDP_cm3<=0)]=np.nan
    sc = plt.scatter(LON, ALT, c=CDP_cm3, cmap=cm,s=10,zorder=1) 
    plt.ylabel('Altitude, m')
    plt.xlabel('Longitude')
    #plt.ylim([0, 0.3])
    #plt.xlim([0, 0.3])
    cbar=plt.colorbar(sc)
    cbar.set_label('CDP, cm$^{-3}$')
    
    plt.subplot(2, 1, 2,sharex=ax1)
    cm = plt.cm.get_cmap('jet')
    plt.plot(LON, ALT,color='gray',markersize=0,zorder=0)
    NC_HI_2DS[(NC_HI_2DS<=0)]=np.nan
    sc = plt.scatter(LON, ALT, c=NC_HI_2DS, norm=colors.LogNorm(vmin=0.01, vmax=100),cmap=cm,s=10,zorder=1) 
    plt.ylabel('Altitude, m')
    plt.xlabel('Longitude')
    #plt.ylim([0, 0.3])
    #plt.xlim([0, 0.3])
    cbar=plt.colorbar(sc)
    cbar.set_label('2DS HI, L$^{-1}$')
    
    #plt.show()
    #plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    Figurename=FlightNumber+'_CrossSection.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)


#_____________________________________________________________________________________________________

#Plot alt vs distance from chilbolton


def PlotCrossSections2CHB(FlightDict,SavePath,FlightNumber):
    
    NC_HI_2DS=FlightDict['NC_HI_2DS'] 
    CDP_cm3=FlightDict['CDP_cm3']
    LON=FlightDict['LON_GIN']
    LAT=FlightDict['LAT_GIN']    
    ALT=FlightDict['ALT_GIN']
       
    ALT[ (LAT > 51.5) | (LON > -1.4) ] =np.nan
    NC_HI_2DS[ (LAT > 51.5) | (LON > -1.4) ] =np.nan
    CDP_cm3[ (LAT > 51.5) | (LON > -1.4) ] =np.nan
    
    Distance2CHB=np.zeros(len(LAT))*np.nan
    
    #Landmark_Lat=[51.145037, 50.819167, 52.069791]
    #Landmark_Lon=[ -1.438416, -2.556333, -0.619135]
    lonCHB = -1.438416
    latCHB = 51.145037
    
    for x in range(len(ALT)-1):    
        Distance2CHB[x] = haversine(LON[x], LAT[x], lonCHB, latCHB)
    
    fig=plt.figure(figsize=(10,10)) 
    plt.rcParams.update({'font.size': 16})
    plt.title(FlightNumber)
    ax1=plt.subplot(2, 1, 1)
    plt.title(FlightNumber)
    cm = plt.cm.get_cmap('jet')
    plt.plot(Distance2CHB, ALT,color='gray',markersize=0,zorder=0)
    CDP_cm3[(CDP_cm3<=0)]=np.nan
    sc = plt.scatter(Distance2CHB, ALT, c=CDP_cm3, cmap=cm,s=10,zorder=1, norm=colors.Normalize(vmin=0, vmax=100)) 
    plt.ylabel('Altitude, m')
    plt.xlabel('Distance from Chilbolton, km')
    #plt.ylim([0, 0.3])
    #plt.xlim([0, 0.3])
    cbar=plt.colorbar(sc)
    cbar.set_label('CDP (3 to 50 μm), cm$^{-3}$')
    
    plt.subplot(2, 1, 2,sharex=ax1)
    cm = plt.cm.get_cmap('jet')
    plt.plot(Distance2CHB, ALT,color='gray',markersize=0,zorder=0)
    NC_HI_2DS[(NC_HI_2DS<=0)]=np.nan
    sc = plt.scatter(Distance2CHB, ALT, c=NC_HI_2DS, norm=colors.LogNorm(vmin=0.01, vmax=100),cmap=cm,s=10,zorder=1) 
    plt.ylabel('Altitude, m')
    plt.xlabel('Distance from Chilbolton, km')
    #plt.ylim([0, 0.3])
    #plt.xlim([0, 0.3])
    cbar=plt.colorbar(sc)
    cbar.set_label('2DS Highly Irregular, L$^{-1}$')
    
    #plt.show()
    #plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    Figurename=FlightNumber+'_CrossSection2CHB.png'
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)




#_____________________________________________________________________________________________________

# Ouput combine 2DS and HVPS size distribution for a given run



#    CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
#    CSVName='flight-sum_faam_20180213_r0_c081.csv'


def SingleRunPSD(RunNumber,FlightDict,FlightTrack_ZED_H,AltitudeBins,ZED_H_grid,CSVPath,CSVName):

    #Load flight log and choose run number
    StartTime, EndTime, RunNumberArray, ProfileNumber = LoadFlightSummaryCSV(CSVPath,CSVName)
    
    #RunNumber=2
    
    tmp=RunNumberArray.tolist()
    i=tmp.index(RunNumber)
    
    print ('Run number '+str(RunNumberArray[i]))
    RunStr='RunNumber'+str(int(RunNumberArray[i]))

    NC_HI_2DS=FlightDict['NC_HI_2DS']
    
    IWC_gm3=FlightDict['IWC_g_m3']  
    Reflectivity_1D= FlightDict['Reflectivity']
    

    #Flag to only average the insitu data when colocated with the radar RHIs
    #RadarFlag=np.zeros(len(FlightTrack_ZED_H))
    #RadarFlag[FlightTrack_ZED_H==FlightTrack_ZED_H]=1 # 1=radar colocated 

    #Flag to only average the insitu data when colocated with the radar tracking
    RadarFlag=np.zeros(len(IWC_gm3))
    RadarFlag=(np.where((Reflectivity_1D==Reflectivity_1D) & (IWC_gm3==IWC_gm3) & (NC_HI_2DS>0),1,0)) 


    Time_Core=FlightDict['Time_Core']
    StartIdx=bisect.bisect_left(Time_Core, StartTime[i]) # assume time is sorted, which it should be 
    EndIdx=bisect.bisect_left(Time_Core, EndTime[i]) # assume time is sorted, which it should be  

    Zed_colocate_avg=np.nanmean(FlightTrack_ZED_H[StartIdx:EndIdx])   
   
    #Prefix='NC_All_accept'    
    #2DS
    NC_SizeMid_2DS, NC_z_dNdDp_2DS= Avg_OAP_dndDp_Dict(FlightDict, 'NC_2DS_size', 'Time_Core', 'NC_HI_z_2DS_psd', StartTime[i],EndTime[i],RadarFlag)    
    #HVPS
    NC_SizeMid_HVPS, NC_z_dNdDp_HVPS= Avg_OAP_dndDp_Dict(FlightDict, 'NC_HVPS_size', 'Time_Core', 'NC_All_accept_HVPS_psd', StartTime[i],EndTime[i],RadarFlag)
    #CombinePSDs
    CompositeSize, Composite_dN, Composite_dNdDp=CompositePSD(NC_SizeMid_2DS, NC_z_dNdDp_2DS, NC_SizeMid_HVPS, NC_z_dNdDp_HVPS)
    
    #plt.plot(CompositeSize, Composite_dN)
    #plt.xscale('log')
    #plt.yscale('log')
    
    #Radar_ref, IWC=CalculateReflectivityIce(CompositeSize, Composite_dN, 7.38E-11, 1.9) # brown and francis
    #Radar_ref, IWC=CalculateReflectivityIce(CompositeSize, Composite_dN, 1.6547E-12, 2.47368)
    
    #Nevzorov
    IWC_avg=np.nanmean(IWC_gm3[StartIdx:EndIdx])
    IWC_gm3[RadarFlag!=1]=np.nan
    IWC_avg_colocate=np.nanmean(IWC_gm3[StartIdx:EndIdx])
      
    #Reflectivity
    #Reflectivity_1D[RadarFlag!=1]=np.nan
    Zed_1D_avg=np.nanmean( Reflectivity_1D[StartIdx:EndIdx])
    
    # Average altitude
    ALT_GIN=FlightDict['ALT_GIN'] 
    ALT_avg=np.nanmean(ALT_GIN[StartIdx:EndIdx])
    Zed_Alt_avg=CalcAltitudeAvgReflectivity(ALT_avg,AltitudeBins,ZED_H_grid)
    
    print('Reflectivity from altitude average = '+str(Zed_Alt_avg)+' dBz')
    print('Reflectivity from colocate = '+str(Zed_colocate_avg)+' dBz')
    print('Reflectivity from tracking = '+str(Zed_1D_avg)+' dBz')
    print('IWC = '+str(IWC_avg)+' g/m3')
    print('IWC (colocate) = '+str(IWC_avg_colocate)+' g/m3')
    
    
    return IWC_avg_colocate, IWC_avg, Zed_Alt_avg, Zed_1D_avg, Zed_colocate_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr


#_____________________________________________________________________________________________________

# Merge 2DS and HVPS PSD. Weighted average between 700-1000 um.    


def CompositePSD(Size2DS, PSD_2DS, SizeHVPS, PSD_HVPS):
    
    BinWidth=10
    
    # <700um 2DS
    # 700 to 1000 um 2DS + HVPS
    # >1000 um HVPS
    
    CompositeSize=np.linspace(10,19200,1920)
    
    CompositeHVPS=np.interp(CompositeSize,SizeHVPS, PSD_HVPS,left=np.nan,right=np.nan)
    Composite2DS=np.interp(CompositeSize,Size2DS, PSD_2DS,left=np.nan,right=np.nan)
    
    Composite_dNdDp=np.zeros(len(CompositeSize))
    
    for i in range(len(CompositeSize)):
        
        if CompositeSize[i]<700 :
            Composite_dNdDp[i]=Composite2DS[i]
            
        if CompositeSize[i]>1000 :
            Composite_dNdDp[i]=CompositeHVPS[i]
    
        if (CompositeSize[i]>=700) & (CompositeSize[i]<=1000):
            ScaleFactor2DS=(1000-CompositeSize[i])/(1000-700)
            ScaleFactorHVPS=1-ScaleFactor2DS
            Composite_dNdDp[i]=ScaleFactor2DS*Composite2DS[i]+ScaleFactorHVPS*CompositeHVPS[i]
    
    Composite_dN=Composite_dNdDp*BinWidth
    
    if 1==2 :
        plt.plot(CompositeSize, CompositeHVPS)
        plt.plot(CompositeSize, Composite2DS)
        plt.plot(CompositeSize, Composite_dNdDp)
        plt.yscale('log')
        plt.xscale('log')
    
    return CompositeSize, Composite_dN, Composite_dNdDp



#_____________________________________________________________________________________________________    
# Calculate MVD. THIS NEEDS CHECKING!!!!!!!!!!!!!!!!!!!!



def CalculateMVD(PSD,MidSize):
    
    MidMass=np.zeros(len(MidSize))
    for i in range(len(MidSize)):
    #    if(MidSize[i]<100):
    #        prefix = 4.82
    #        exponent = -13
    #        power = 3.0
    #    else:
    #        prefix = 7.38
    #        exponent = -11
    #        power = 1.9
    #    MidMass[i] = prefix * (MidSize[i]**power) * (10**exponent)	
        MidMass[i] = 1/6 * np.pi * MidSize[i] ** 3 # Convert D to radius
    
    PSD_mass= PSD * np.transpose(MidMass)      
    TotalMass = np.sum(PSD_mass,axis=1)
    
    #CumPSD=np.zeros(np.size(PSD_mass, axis=0),np.size(PSD_mass, axis=1))
    
    Pro= PSD_mass / TotalMass[:,None]
    CumPSD= np.cumsum(PSD_mass, axis=1,)/ TotalMass[:,None]
    
    #PSD_mass / TotalMass[:,None]
    Idx=np.zeros(np.size(PSD_mass, axis=0))*np.nan
    Idx=np.argmax(CumPSD>0.5,1) #- (np.amin(CumPSD<=0.5,1))
    
    #Idx[Idx==0]=np.nan
    #MVD = MidSize[Idx]+ ((0.5 - CumPSD[Idx-1])/Pro[Idx])*(MidSize[Idx+1] - MidSize[Idx])
    MVD=np.zeros(len(Idx))*np.nan
    
    for i in range(len(Idx)):    
        if Idx[i]==Idx[i]:
            MVD[i]= MidSize[Idx[i]-1]+((0.5-CumPSD[i,Idx[i]-1])/(CumPSD[i,Idx[i]]-CumPSD[i,Idx[i]-1]))*(MidSize[Idx[i]]-MidSize[Idx[i]-1])
    
    return MVD #CumPSD, Idx 

#_____________________________________________________________________________________________________

#load radar RHI as 1D array. Distance_1D is distance from chilbolton.


def LoadRadar(Filename):
    Path='C:/Users/Admin TEMP/Documents/PICASSO/RadarData/20180213_v1/'
    #Filename='20180213_CAMRA_run7_up.nc'
    

    RadarData = Dataset(Path+Filename)
    print(RadarData.variables.keys())

    Range=np.array(RadarData['range'][:])
    Range/=1000 # Convert to Km
    
    Elevation=np.array(RadarData['elevation'][:])
    ZED_H=np.array(RadarData['ZED_H'][:])
    Height=RadarData['height']
    ChilboltonHeight=Height[0]/1000
   
    ZED_H_1D=np.zeros(len(Range)*len(Elevation))*np.nan
    Distance_1D=np.zeros(len(Range)*len(Elevation))*np.nan 
    Altitude_1D=np.zeros(len(Range)*len(Elevation))*np.nan 
    
    for i in range(len(Range)):
        for j in range(len(Elevation)):
            ZED_H_1D[i+(j*len(Range))]=ZED_H[i,j]
            Distance_1D[i+(j*len(Range))]=getDist(Range[i],Elevation[j])
            Altitude_1D[i+(j*len(Range))]=getAlt(Range[i],Elevation[j])+ChilboltonHeight
    

    return RadarData, ZED_H_1D,Distance_1D,Altitude_1D, Range, Elevation





#_____________________________________________________________________________________________________
    
#def GridRadar(ZED_H_1D,Distance_1D,Altitude_1D,Distance2CHB, ALT):      
#    ZED_H_grid, DistanceBins, AltitudeBins = BinAverage2D(Distance_1D,Altitude_1D,ZED_H_1D, 0,180,360,0,10, 500)
#    plt.pcolor(DistanceBins,AltitudeBins,np.transpose(ZED_H_grid), cmap='jet', vmin=-20, vmax=30)
#    plt.plot(Distance2CHB,ALTRun,'+',color='black')

#    return ZED_H_grid, DistanceBins, AltitudeBins

#_____________________________________________________________________________________________________
    
# Average the reflectivity at altitude Alt using the bins AltitudeBins


def CalcAltitudeAvgReflectivity(Alt,AltitudeBins,ZED_H_grid):
    Alt/=1000 #convert to km
    
    Alt_idx = np.digitize(Alt, AltitudeBins)
    ZED_H_AltAvg=np.nanmean(ZED_H_grid[:,Alt_idx-1])

    #print(Alt)
    #print(AltitudeBins[Alt_idx-1])
    #print(Alt_idx)
    #print(ZED_H_AltAvg)

    return ZED_H_AltAvg
#_____________________________________________________________________________________________________  

# Co locate aircraft with radar RHI. The radar and aircraft data is put on the same 2D grid. Radar then extracted along the flight path.
# Means there are data gaps in the extracted radar time series.

def Colocate(FlightDict,ZED_H_1D,Distance_1D,Altitude_1D,PltTitle):
        
    #CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
    #CSVName='flight-sum_faam_20180213_r0_c081.csv'    
    #StartTime, EndTime, RunNumber, ProfileNumber = LoadFlightSummaryCSV(CSVPath,CSVName)
    #i=2
    
    Time_Core=FlightDict['Time_Core']
    LON=FlightDict['LON_GIN']
    LAT=FlightDict['LAT_GIN']    
    ALT=FlightDict['ALT_GIN']
    
    ALT_km=ALT/1000 #convert to km
    
    
    #StartIdx=bisect.bisect_left(Time_Core, StartTime[i]) # assume time is sorted, which it should be 
    #EndIdx=bisect.bisect_left(Time_Core, EndTime[i]) # assume time is sorted, which it should be 
    
    #TimeRun=Time_Core[StartIdx:EndIdx]
    #LONRun=LON[StartIdx:EndIdx]
    #LATRun=LAT[StartIdx:EndIdx]
    #ALTRun=ALT[StartIdx:EndIdx]    
    
    LON[ (LAT > 51.15) | (LON > -1.43) ] =np.nan
    LAT[ (LAT > 51.15) | (LON > -1.43) ] =np.nan
    ALT_km[ (LAT > 51.15) | (LON > -1.43) ] =np.nan
    
    Distance2CHB=np.zeros(len(LON))*np.nan
    FlightTrack_ZED_H=np.zeros(len(LON))*np.nan
    
    lonCHB = -1.438416
    latCHB = 51.144501
    
    for x in range(len(LAT)-1):    
        Distance2CHB[x] = haversine(LON[x], LAT[x], lonCHB, latCHB)
            
    ZED_H_grid, DistanceBins, AltitudeBins = BinAverage2D(Distance_1D,Altitude_1D,ZED_H_1D, 0,180,180,0,10,500)
        
    digitizedX = np.digitize(Distance2CHB, DistanceBins)
    digitizedY = np.digitize(ALT_km, AltitudeBins)

    FlightTrack_ZED_H= ZED_H_grid[digitizedX,digitizedY]
   
    

    
    if 1==2 :    
        #font = {'family' : 'normal',
        #'weight' : 'bold',
        #'size'   : 16}

        #plt.rc('font', **font)
        plt.rcParams.update({'font.size': 16})
        #plt.style.use('presentation')
        plt.style.use('classic')
        fig=plt.figure(figsize=(10,6))
        plt.pcolor(DistanceBins,AltitudeBins,np.transpose(ZED_H_grid), cmap='jet', vmin=-20, vmax=30)
        plt.plot(Distance2CHB,ALT_km,markersize=0,linewidth=2,color='black')
        plt.ylabel('Altitude, km')
        plt.xlabel('Distance from Chilbolton, km')
        cbar=plt.colorbar(orientation='horizontal')
        cbar.set_label('Reflectivity, dBZ')
        
        
        FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'
        #plt.title(PltTitle)
        #plt.savefig(FilePath+PltTitle+'.png',dpi=200)
        #plt.close(fig)
        plt.show()
        
        
    return Distance2CHB, FlightTrack_ZED_H, LAT, LON, ALT_km,Time_Core, ZED_H_grid, DistanceBins, AltitudeBins       
    
#_____________________________________________________________________________________________________  
# Plot radar    
    
def PlotRadar(ZED_H_1D,Distance_1D,Altitude_1D):   
    cm = plt.cm.get_cmap('RdYlBu')
    plt.scatter(Distance_1D, Altitude_1D, c=ZED_H_1D,vmin=-20, vmax=20,cmap=cm)
    plt.ylim([0,12])
    #plt.pcolor(Distance,Altitude,ZED_H, cmap='jet')

#_____________________________________________________________________________________________________
# Gets horizontal distance from source for a target at given range and elevation
def getDist(Range,Elevation):
	
    r_e=6371	#earth radius in km
    #Distance=np.zeros(len(Elevation))
    Distance= r_e * np.arctan(	( Range*np.cos(Elevation*math.pi/180) ) / ( r_e + (Range*np.sin(Elevation*math.pi/180) ) )	)
    return Distance

#_____________________________________________________________________________________________________
# Gets altitude above source level for a target at given range and elevation

def getAlt(Range,Elevation):
	
    r_e=6371	#earth radius in km
    Altitude = np.sqrt(	(r_e**2)	+	(Range**2)	+	(2*r_e*Range*np.sin(Elevation*math.pi/180))	) - r_e
    return Altitude

#_____________________________________________________________________________________________________

# Plot Temperature vs HI box plots

def Plot_TvsHI(FlightDict,FilePath,PltTitle):
    
    #FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'
    #PltTitle='C081_HI_IWC'
    
    NC_HI_2DS=FlightDict['NC_HI_2DS'] 
    TAT_ND_R_C=FlightDict['TAT_ND_R_C']
    IWC_g_m3=FlightDict['IWC_g_m3']
    
    NC_HI_2DS[NC_HI_2DS==0]=np.nan
    IWC_g_m3[NC_HI_2DS==0 ]=np.nan
    IWC_g_m3[IWC_g_m3<=0]=np.nan
    IWC_g_m3[NC_HI_2DS != NC_HI_2DS]=np.nan
    
    NC_HI_2DS_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(NC_HI_2DS,TAT_ND_R_C,-36,0,18)
    mask = ~np.isnan(NC_HI_2DS_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(NC_HI_2DS_bin.T, mask.T)]
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 16})
    ax1=fig.add_subplot(1, 2, 1)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid,vert=0)
    #means = [np.nanmean(x) for x in filtered_data]
    #plt.T.boxplot(vert=False)
    plt.ylabel('Temperature, °C')
    plt.xlabel('2DS Highly Irregular, L$^{-1}$')
    #ax1.tick_params(axis='top', which='both', length=0)
    
    if 1==1 :
        CooperCurve=np.zeros(len(BinsMid))
        Cooper_INP_array = np.vectorize(Cooper_INP)
        CooperCurve=Cooper_INP_array(BinsMid)
        plt.plot(CooperCurve, BinsMid,linewidth=1, color='tab:red', label ='INP Cooper (1986)')
        plt.gca().set_xlim([0,25])
        ax1.legend()
    
    #ax2 = ax1.twiny()
    plt.subplot(1, 2, 2, sharey=ax1)

    IWC_g_m3_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(IWC_g_m3,TAT_ND_R_C,-36,0,18)
    mask = ~np.isnan(IWC_g_m3_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(IWC_g_m3_bin.T, mask.T)]  
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid,vert=0)
    plt.xlabel('IWC, g m$^{-3}$')
    
    plt.gca().set_ylim([0,-40])
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().yaxis.set_major_locator(majorLocator)
    plt.gca().yaxis.set_major_formatter(majorFormatter)
    plt.gca().yaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter


    #plt.title(PltTitle)
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    
    
#_____________________________________________________________________________________________________

# Calculate number of INP from Cooper 1986


def Cooper_INP(T_C):    
    T=T_C+273.15    
    if T < 233 :
        INcooper = 0.005*math.exp(0.304*(273.15-233))
    else :
        INcooper = 0.005*math.exp(0.304*(273.15-T))
        
    return INcooper
    

#____________________________________________________________________________________________________

# Add tracking radar .mat files to flight dictionary and put on the time base Time_Core 

def AddRadar2FlightDict(FlightDict):

    Time_Core=FlightDict['Time_Core']
    BatchTimeRadar,BatchReflectivity=BatchRadarMat()
        
    ReflectivityCore=ChangeTimeBaseAvg(BatchTimeRadar,BatchReflectivity,Time_Core,1)
    FlightDict['Reflectivity']= ReflectivityCore

    #plt.plot(BatchTimeRadar,BatchReflectivity)
    #plt.plot(Time_Core,ReflectivityCore)
#____________________________________________________________________________________________________

# Load tracking radar data

def BatchRadarMat():

    BatchTimeRadar = []
    BatchReflectivity = [] 
    #Path='C:/Users/Admin TEMP/Documents/PICASSO/RadarData/Time series/'
    Path='C:/Users/Admin TEMP/Documents/PICASSO/RadarData/C081_TimeSeries/'
    for filena in os.listdir(Path):
        if filena.endswith(".mat"):     
            print(Path)
            print(filena)
            TimeRadar,Reflectivity = LoadRadarMatFiles(Path,filena)
            BatchTimeRadar = np.append(BatchTimeRadar,TimeRadar)
            BatchReflectivity = np.append(BatchReflectivity,Reflectivity) 
    
    return BatchTimeRadar,BatchReflectivity


#____________________________________________________________________________________________________

# Load individual tracking radar .mat file

def LoadRadarMatFiles(Path,FileName):
    #Path='C:/Users/Admin TEMP/Documents/PICASSO/RadarData/Time series/'
    #FileName='Timeseries_Down_20180213065439.mat'
    
    matdata = loadmat(Path+FileName)
    Reflectivity=matdata['data']['Reflectivity']
    matlab_datenum=np.array(matdata['data']['TimeDatenum'])
    TimeRadar=Matlab2PythonTimeArray(matlab_datenum)
    
    return TimeRadar,Reflectivity 

    #plt.plot(TimeRadar,Reflectivity)
  
    #plt.show()
    #mat_contents = sio.loadmat(Path+FileName)    
    #data_struct = mat_contents['data']
    #print(mat_contents.keys())
    #return Test
    #Reflectivity = data_struct['Reflectivity'][0,0]
    #return Reflectivity 


