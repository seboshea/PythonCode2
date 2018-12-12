# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:39:15 2018

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
import gc
from PICASSO_functions import Avg_dndDp, Avg_CDP_dndDp,Avg_dndDp_CountStats
from MyFunctions import LoadCoreCloud, KorolevCorrectedD, GetKorolevRatios_Zd, ReuterBakanDiameter,DimensionlessDistance
from MyFunctions import BinLinearV2, Load_HaloHolo_PSD, BinHaloHolo_v2, datetime_2_s_sincemidnight

#_______________________________________________________________________________________

def LoadGrayScaleStats(GSpath,ExpDate):

    # Load grayscale stats
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/60um_25_50_75/20171018/20171018091312/OASIS_Out/'
    for filena in os.listdir(GSpath):
        if filena.endswith(".h5"):
            print(filena)
            Data_h5 = h5py.File(GSpath+filena, 'r')              
            TmpTimesWv=np.array(Data_h5['ParticleTimesWv'])
            #TmpStats=np.array(Data_h5['PartGrayStatsWv'])
            TmpStats=np.array(Data_h5['ParticleStatsWv'])
            Data_h5.close()
            gc.collect()       
            if 'ParticleTimesWv' not in dir():
                ParticleTimesWv=TmpTimesWv
                #PartGrayStatsWv=TmpStats
                ParticleStatsWv=TmpStats
            else:    
                ParticleTimesWv= np.append(ParticleTimesWv,TmpTimesWv,axis=0)
                #PartGrayStatsWv= np.append(PartGrayStatsWv,TmpStats,axis=0)
                ParticleStatsWv= np.append(ParticleStatsWv,TmpStats,axis=0)
            
            
    #ExpDate= datetime.datetime(2017, 10, 18, 0, 0, 0)
    ParticleTime=[ExpDate + datetime.timedelta(seconds=int(ParticleTimesWv[x][0])) for x in range(len(ParticleTimesWv))]
    

    #plt.plt(StageTime,StageY )
    return ParticleStatsWv, ParticleTime,ParticleTimesWv


#_______________________________________________________________________________________

# Calculate greyscale ratios

# Need to add stage positon flag, so can accept stage position array. If stage position present need to nan Stats when missing stage data.  

# If don't have stage positions call GrayscaleRatios(0, 0, ParticleStatsWv,0). When stage positions are included ParticleStatsWv is nanned where missing stage data

# Constant to define the single particle stats output
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

 

def GrayscaleRatios(Particle_X, Particle_Y, ParticleStatsWv,MeanXY, ThresholdFlag):    

    #Remove edge particles and where there are gaps in the stage data
    #PartGrayStats=np.array(ParticleStatsWv)
    
    
    if isinstance(Particle_X,int):
        for i in range(len(ParticleStatsWv)):
            if (ParticleStatsWv[i,23]==0 or ParticleStatsWv[i,27]==0 or ParticleStatsWv[i,31]==0 or ParticleStatsWv[i,24]==63 or ParticleStatsWv[i,28]==63 or ParticleStatsWv[i,32]==63):
                ParticleStatsWv[i,:]=np.nan
    else :
        for i in range(len(ParticleStatsWv)):
            if (ParticleStatsWv[i,23]==0 or ParticleStatsWv[i,27]==0 or ParticleStatsWv[i,31]==0 or ParticleStatsWv[i,24]==63 or ParticleStatsWv[i,28]==63 or ParticleStatsWv[i,32]==63 or Particle_X[i]!=Particle_X[i]):            
                ParticleStatsWv[i,:]=np.nan

    #for i in range(len(ParticleStatsWv)):    
    #    if (ParticleStatsWv[i,23]==0 or ParticleStatsWv[i,27]==0 or ParticleStatsWv[i,31]==0 or ParticleStatsWv[i,24]==63 or ParticleStatsWv[i,28]==63 or ParticleStatsWv[i,32]==63 or Particle_X[i]!=Particle_X[i]):
    #        ParticleStatsWv[i,:]=np.nan



    FilledArea=ParticleStatsWv[:,1]
    VoidArea=ParticleStatsWv[:,15]
    
    TotalArea=ParticleStatsWv[:,22]+ParticleStatsWv[:,21]+ParticleStatsWv[:,20]
    AreaFraction0=ParticleStatsWv[:,20]/TotalArea[:]
    AreaFraction1=ParticleStatsWv[:,21]/TotalArea[:]
    AreaFraction2=ParticleStatsWv[:,22]/TotalArea[:]
    
    #AreaFraction0=ParticleStatsWv[:,20]/TotalArea[:]
    #AreaFraction1=(ParticleStatsWv[:,21]+ParticleStatsWv[:,22])/TotalArea[:]
    #AreaFraction2=ParticleStatsWv[:,22]/TotalArea[:]
    
    AreaRatio_2_1=ParticleStatsWv[:,22]/ParticleStatsWv[:,21]
    AreaRatio_1_0=ParticleStatsWv[:,21]/ParticleStatsWv[:,20] 
        
    #Use mean x y for diameter
    #DiameterLevel0=15*(PartGrayStats[:,4]-PartGrayStats[:,3]+PartGrayStats[:,6]-PartGrayStats[:,5])/2
    #DiameterLevel1=15*(PartGrayStats[:,8]-PartGrayStats[:,7]+PartGrayStats[:,10]-PartGrayStats[:,9])/2
    #DiameterLevel2=15* (PartGrayStats[:,12]-PartGrayStats[:,11]+PartGrayStats[:,14]-PartGrayStats[:,13])/2
    
    # Just use x for diameter
    #DiameterLevel0=15*(ParticleStatsWv[:,24]-ParticleStatsWv[:,23])
    #DiameterLevel1=15*(ParticleStatsWv[:,28]-ParticleStatsWv[:,27])
    #DiameterLevel2=15* (ParticleStatsWv[:,32]-ParticleStatsWv[:,31])

    # 15+ needed if 1 pixel then 0=(Max-Min)
    if MeanXY==1:
        DiameterLevel0x=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32]),ParticleStatsWv[:,24])-np.minimum(np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]),ParticleStatsWv[:,23]))
        DiameterLevel1x=15+15*(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32])-np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]))
        DiameterLevel2x=15+15*(ParticleStatsWv[:,32]-ParticleStatsWv[:,31])    
        DiameterBGx=15+15*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
        
        DiameterLevel0y=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,30],ParticleStatsWv[:,34]),ParticleStatsWv[:,26])-np.minimum(np.minimum(ParticleStatsWv[:,29],ParticleStatsWv[:,33]),ParticleStatsWv[:,25]))
        DiameterLevel1y=15+15*(np.maximum(ParticleStatsWv[:,30],ParticleStatsWv[:,34])-np.minimum(ParticleStatsWv[:,29],ParticleStatsWv[:,33]))
        DiameterLevel2y=15+15*(ParticleStatsWv[:,34]-ParticleStatsWv[:,33])    
        DiameterBGy=15+15*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11]) #x diameter of largerst paricle in bounding box 

        DiameterLevel0= (DiameterLevel0x + DiameterLevel0y)/2
        DiameterLevel1= (DiameterLevel1x + DiameterLevel1y)/2
        DiameterLevel2= (DiameterLevel2x + DiameterLevel2y)/2  
        DiameterBG= (DiameterBGx + DiameterBGy)/2        
              
    else:
        DiameterLevel0=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32]),ParticleStatsWv[:,24])-np.minimum(np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]),ParticleStatsWv[:,23]))
        DiameterLevel1=15+15*(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32])-np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]))
        DiameterLevel2=15+15*(ParticleStatsWv[:,32]-ParticleStatsWv[:,31])    
        DiameterBG=15+15*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    
    
    
    #AspectRatio0=(PartGrayStatsWv[:,4]-PartGrayStatsWv[:,3])/(PartGrayStatsWv[:,6]-PartGrayStatsWv[:,5])

    DiameterLevel0[DiameterLevel0<0]= np.nan
    DiameterLevel1[DiameterLevel1<0]= np.nan
    DiameterLevel2[DiameterLevel2<0]= np.nan
        
    Zd_fromRatios=GetZd(AreaRatio_1_0,AreaFraction2,AreaFraction1,AreaFraction0, ThresholdFlag)
    
    return AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG, Zd_fromRatios




#_______________________________________________________________________________________

def GetIAT(ParticleTimesWv,row,NumRow):
    
	
    t1= -np.inf
    t2= ParticleTimesWv[row,0]+(ParticleTimesWv[row,1]/1000000000)
    t3= np.inf
	
    if(row>0):
        t1=ParticleTimesWv[row-1,0]+(ParticleTimesWv[row-1,1]/1000000000)
	
    if(row<(NumRow-1)):
        t3=ParticleTimesWv[row+1,0]+(ParticleTimesWv[row+1,1]/1000000000)
        
    IAT1=t2-t1
    IAT2=t3-t2
	
    return min(IAT1,IAT2)

#_______________________________________________________________________________________  

    
def GetIAT_vector(ParticleTime,ParticleTimesWv):
    IAT=np.zeros(len(ParticleTime))
    for i in range (len(IAT)):
        IAT[i]=GetIAT(ParticleTimesWv,i,len(IAT)-1)
    return IAT



#_______________________________________________________________________________________
    
def GetZd(AreaRatio_1_0,AreaFraction2,AreaFraction1,AreaFraction0,ThresholdFlag):
    
    Zd=np.zeros(len(AreaRatio_1_0))*np.nan 
    
    #Zd=np.where(AreaRatio_1_0, np.where(AreaFraction2==0,3.5,0), np.nan
    
    if ThresholdFlag == 255075 : 
    
        for i in range(len(AreaRatio_1_0)):
            if (AreaRatio_1_0[i] == AreaRatio_1_0[i]) & (AreaFraction2[i]==AreaFraction2[i]) :
                if AreaFraction2[i] > 0 :
                    Zd[i]=3.5
                else :
                    if AreaFraction1[i] == 0 :
                        Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                    else : 
                        Zd[i]= (AreaRatio_1_0[i] - 3.2048) / (-0.3772)
                
                
                
                #if AreaRatio_1_0[i] == 0 :  
                #    Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                #else:
                #    if AreaFraction2[i] == 0 :
                #        Zd[i]= (AreaRatio_1_0[i] - 3.2048) / (-0.3772)
                #    else : 
                #        Zd[i]=3.5
    
    if ThresholdFlag == 405070 : 
    
        for i in range(len(AreaRatio_1_0)): 
            if (AreaRatio_1_0[i] == AreaRatio_1_0[i]) & (AreaFraction2[i]==AreaFraction2[i]) :
                if AreaFraction2[i] > 0 :    
                    Zd[i]=3.5
                else:
                    if AreaFraction1[i] == 0 :
                        Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                    else:
                        Zd[i]= (AreaRatio_1_0[i] - 7.8282) / (-0.9507)
                        
                
            
            #if (AreaRatio_1_0[i] == AreaRatio_1_0[i]) & (AreaFraction2[i]==AreaFraction2[i]) :
                #if AreaRatio_1_0[i] == 0 :  
                #    Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                #else:
                #    if AreaFraction2[i] == 0 :
                #        Zd[i]= (AreaRatio_1_0[i] - 7.8282) / (-0.9507)
                #    else : 
                #        Zd[i]=3.5    
    
    if ThresholdFlag == 255067 : 
    
        for i in range(len(AreaRatio_1_0)):
            if (AreaRatio_1_0[i] == AreaRatio_1_0[i]) & (AreaFraction2[i]==AreaFraction2[i]) :
                if AreaFraction2[i] > 0 :
                    Zd[i]=3.5
                else:
                    if AreaFraction1[i] == 0 :
                        Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                    else :
                        Zd[i]= (AreaRatio_1_0[i] - 3.3619) / (-0.4009)
                        
                
                #if AreaRatio_1_0[i] == 0 :  
                #    Zd[i]= 8.5 # if AreaRatio_1_0 == 0 Zd will be greater than approx 8.5
                #else:
                #    if AreaFraction2[i] == 0 :
                #        Zd[i]= (AreaRatio_1_0[i] - 3.3619) / (-0.4009)
                #    else : 
                #        Zd[i]=3.5   
    
    
    
    return Zd



#_______________________________________________________________________________________
# This won't work if the data spans more than 172800 seconds
 
# if MaxParticlesPerImage==0 don't use this to filter particles

#NFrac = StatWv[row][K_OASIS_PSTAT_COL_AREA]/(StatWv[row][K_OASIS_PSTAT_COL_AREA]+StatWv[row][K_OASIS_PSTAT_COL_AREA2])

def ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,Diameter,Zd_fromRatios,Zd_max,DOF_c,MaxParticlesPerImage, ArmSep, DLow_D0_Zd_Table):
    
    IATmin=1E-5
    #IATmin=0
    #Zd_max=6 # maximum Zd to include in PSD
    #DOF_c= 8 # c value for DOF calculation, normally should be same as Zd_max
    TAS=100
    Lambda = 0.658	# in um 
    Lambda_mm=  Lambda / 1000
    ProbeRes= 15
    #ArmSep = 70
    ArrayElements = 64
    
    NParticles=np.array(ParticleStatsWv[:,0])

    # To calculate effective array width need uncorrected image x diameter
    #DiameterLevel0x=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32]),ParticleStatsWv[:,24])-np.minimum(np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]),ParticleStatsWv[:,23]))
    
    StartDate=datetime.datetime.combine(ParticleTime[0].date(),datetime.time(second=0))
    IAT= GetIAT_vector(ParticleTime,ParticleTimesWv)    
    ParticleTime_sinceMidnight= ParticleTimesWv[:,0]+(ParticleTimesWv[:,1]/1000000000)
    Startime=int(min(ParticleTimesWv[:,0]))
    Endtime=int(max(ParticleTimesWv[:,0]))
    
    if Endtime > 172800 : # 
        Endtime = 172800
    Ns=Endtime-Startime
    
    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
    
    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
    PSD_SizeMid=SizeBins-ProbeRes*0.5
    
    Time_idx = np.digitize(ParticleTime_sinceMidnight, TimeBins_s)
    Size_idx = np.digitize(Diameter, SizeBins)
    
    # Calculate concentration of individual particle
    DLow_D0_model, Zd_model=Get_DLow_D0_model(DLow_D0_Zd_Table)
    ConcentrationL= np.zeros(len(ParticleTime_sinceMidnight))
    for i in range(len(ParticleTime_sinceMidnight)):
        #ConcentrationL[i]=1/Allin_SVol_lawson(Diameter[i],DiameterLevel0x[i],Lambda,DOF_c,TAS,ProbeRes,ArmSep, ArrayElements)
        ConcentrationL[i]=1/sVol_EAW_INT(TAS,Diameter[i],DOF_c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model)
    
    # Filter and bin data
    Counts_PSD = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L = np.zeros((len(TimeBins_s), len(SizeBins))) 
    
    if MaxParticlesPerImage==0 : 
        for i in range(len(Diameter)-1) : 
            if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (Diameter[i] > (0.5*ProbeRes)) & (Diameter[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) & (Zd_fromRatios[i] <Zd_max): 
                #if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (DiameterLevel1[i] > (0.5*ProbeRes)) & (DiameterLevel1[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) :
                Counts_PSD[Time_idx[i],Size_idx[i]]+=1
                PSD_L[Time_idx[i],Size_idx[i]]+= ConcentrationL[i]
    else : 
        for i in range(len(Diameter)-1) : 
            if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (Diameter[i] > (0.5*ProbeRes)) & (Diameter[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) & (Zd_fromRatios[i] <Zd_max) & (NParticles[i] <= MaxParticlesPerImage): 
                #if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (DiameterLevel1[i] > (0.5*ProbeRes)) & (DiameterLevel1[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) :
                Counts_PSD[Time_idx[i],Size_idx[i]]+=1    
                PSD_L[Time_idx[i],Size_idx[i]]+= ConcentrationL[i]
    
    #SVol=np.zeros(len(SizeBins))
    #DoF=np.zeros(len(SizeBins))
    #for i in range(len(SizeBins)-1) : 
    #    SVol[i],DoF[i]=Allin_SVol_lawson(PSD_SizeMid[i],Lambda,DOF_c,TAS,ProbeRes,ArmSep, ArrayElements)
        
    #PSD_L=Counts_PSD/SVol[None,:]
   
    
    #PSD_MidTime=[datetime.datetime.combine(StartDate,datetime.time(second=int(TimeBins_s[x]-0.5))) for x in range(len(TimeBins_s))]
    PSD_MidTime=[StartDate + datetime.timedelta(seconds=int(TimeBins_s[x]-0.5)) for x in range(len(TimeBins_s))]
   
    return PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime#, SVol,DoF
    
    #Zd = GetZd(AreaRatio_1_0,AreaFraction2)
    #return Zd_fromRatios, IAT
     
#_______________________________________________________________________________________
## From lawson et al 2006  equation 1
# Lambda in um, CIP15 Lambda = 0.658	# in um
#Diameter in um    
#TAS in m/s
#ArmSep in mm 

def Allin_SVol_lawson(ParticleDiameter,DiameterLevel0x, Lambda,c,TAS,ProbeRes,ArmSep, ArrayElements):
    radius=ParticleDiameter/2 # um
    DoF=(2*c*radius**2)/Lambda ## 2* because there is a +/- in the equation
    DoF/=1000 # return in mm
    
    if DoF > ArmSep : 
        DoF = ArmSep
    
    ArrayWidth = ((ArrayElements*ProbeRes) - (DiameterLevel0x))*0.001 # mm
    
    TAS*=1000 # mm/s
    
    sVol = DoF * ArrayWidth * TAS # mm/s
   
    sVol /= 1E6 #L/s
    
    return sVol#, DoF


#_______________________________________________________________________________________

def DvsSvol(DLow_D0_Zd_Table):

    DLow_D0_model, Zd_model=Get_DLow_D0_model(DLow_D0_Zd_Table)
    
    D_array=np.linspace(15,960,64)
    SVol_array= np.zeros(64)

    for i in range(len(SVol_array)):
        #ConcentrationL[i]=1/Allin_SVol_lawson(Diameter[i],DiameterLevel0x[i],Lambda,DOF_c,TAS,ProbeRes,ArmSep, ArrayElements)
        SVol_array[i]=sVol_EAW_INT(100,D_array[i],6,64,15,0.000658, 70, DLow_D0_model, Zd_model)
    
    #SVol_array=sVol_EAW_INT(100,D_array[any()],6,64,15,0.000658, 70, DLow_D0_model, Zd_model)
    
    plt.plot(D_array, SVol_array)
    plt.xscale('log')
    #plt.yscale('log')
    
    
#_______________________________________________________________________________________
#TAS in m/s

def sVol_EAW_INT(TAS,D0,c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model):
    
    SA_mm2=SA_EAW_INT(D0,c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model)
    TAS*=1000 # mm/s
    
    sVol = SA_mm2 * TAS # mm/s
   
    sVol /= 1E6 #L/s
    
    return sVol#, DoF
  
#_______________________________________________________________________________________
# Needed for Svol calculation

def Get_DLow_D0_model(DLow_D0_Zd_Table):
    #DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5'
    Data_h5 = h5py.File(DLow_D0_Zd_Table, 'r')              
    DLow_D0_model=np.array(Data_h5['DLow_D0'])
    Zd_model=np.array(Data_h5['Zd'])
    Data_h5.close()
    return DLow_D0_model, Zd_model

#_______________________________________________________________________________________
#AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre

def SA_EAW_INT(D0,c,ArrayElements,ProbeRes,Lambda_mm, ArmSep,DLow_D0_model,Zd_model):
    
    #D25_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5'
    #Data_h5 = h5py.File(D25_D0_Zd_Table, 'r')              
    #D25_D0_model=np.array(Data_h5['D25_D0'])
    #Zd_model=np.array(Data_h5['Zd'])
    #Data_h5.close()
    
    #D25_D0_model
    #Zd_model
    
    
    DLow_model = DLow_D0_model * D0     
    #Lambda_mm = Lambda/1000	# in mm
    D0_mm=D0/1E3 # convert to mm    
    Z_model = (Zd_model * D0_mm**2) / (4 * Lambda_mm)
    
    EAW = ((ArrayElements*ProbeRes) - (DLow_model))*0.001 # mm
    EAW = np.where(EAW<0, 0, EAW) # remove zeros
    EAW = np.where( np.absolute(Zd_model) > c, 0, EAW) # remove positions outside DoF
    EAW = np.where( np.absolute(Z_model) > (ArmSep/2), 0, EAW) # remove positions outside Arms
    
    #plt.plot(Z_model, EAW, 'o')       
    SA = np.trapz(EAW,Z_model) ## Integrate EAW over Z 
    #print ('Integration')
    #print(SA)
    
    #DoF=(c*D0_mm**2)/(4*Lambda_mm) ## not included +/- in the equation
    #if DoF > (ArmSep/2) : 
    #    DoF = (ArmSep/2)
    #print('DoF')
    #print(DoF)

    #SA_NoImageSizeCorrection= (2*DoF) * ((ArrayElements*ProbeRes) - (50))*0.001 # mm
    #print('SA no size correction')
    #print(SA_NoImageSizeCorrection)
    
    return SA
#_______________________________________________________________________________________

#AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre


def DLow_D0_Zd_Table(Model_Stats, SaveFile, FileName):
    
    Lambda = 0.658
    CentreDof= 50    
    #D0_array=np.linspace(10,)
    #D0 = 90 
    
    for D0 in range(50,155,5): 
        DLow = np.array(Model_Stats[str(D0)][5][:]) #um
        Z = np.array(Model_Stats[str(D0)][6][:])-CentreDof #mm
        Zd_tmp= DimensionlessDistance(Z,D0,0)
        DLow_D0_tmp= DLow / D0
        if D0 == 50 : 
            DLow_D0=DLow_D0_tmp
            Zd= Zd_tmp            
        else :
            DLow_D0=np.append(DLow_D0, DLow_D0_tmp)
            Zd= np.append(Zd, Zd_tmp)
    
    D_2stage_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DLow_D0,Zd,-15.125,15.125,121)
    mask = ~np.isnan(D_2stage_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_2stage_bin.T, mask.T)]
    DLow_D0_means = [np.nanmean(x) for x in filtered_data]
    
    
    BinsMid_int= np.linspace(-15,15,3001)
    DLow_D0_means_int=np.interp(BinsMid_int,BinsMid,DLow_D0_means )
    
    plt.plot(BinsMid,DLow_D0_means,'o', color='b')
    plt.plot(BinsMid_int,DLow_D0_means_int,'o', color='r')
    #plt.xlim([-1,1])
    
    if SaveFile == 1 : 
        #SavePath= 'C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5'
        SavePath= 'C:/Users/Admin TEMP/Documents/DropletGun/'+FileName
        file = h5py.File(SavePath, 'w')
        file.create_dataset('DLow_D0', data=DLow_D0_means_int)
        file.create_dataset('Zd', data=BinsMid_int)
        file.close()
    

    
    
    return DLow_D0_means_int,  BinsMid_int

#_______________________________________________________________________________________
    
## Clarify flight C052  ********************************
    
def Clarify_C052_grayscalePSD():
    DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5'
    ExpDate= datetime.datetime(2017, 9, 5, 0, 0, 0)
    
    GSpath='C:/Users/Admin TEMP/Documents/Clarify/C052/OasisOut/Imagefile_1CIP Grayscale_20170905163107/'
    PltTitle='C052_1642'
    StartAvg=datetime.datetime(2017, 9, 5, 16, 42, 10)
    EndAvg=datetime.datetime(2017, 9, 5, 16, 43, 15)
    #EndAvg=datetime.datetime(2017, 9, 5, 16, 42, 33)
    
    
    #GSpath='C:/Users/Admin TEMP/Documents/Clarify/C052/OasisOut/Imagefile_1CIP Grayscale_20170905140314/'
    #PltTitle='C052_1608'
    #StartAvg=datetime.datetime(2017, 9, 5, 16, 8, 40)
    #StartAvg=datetime.datetime(2017, 9, 5, 16, 10, 11)
    #EndAvg=datetime.datetime(2017, 9, 5, 16, 11, 10)
    
    #PltTitle='C052_1556'
    #StartAvg=datetime.datetime(2017, 9, 5, 15, 55, 56)
    #StartAvg=datetime.datetime(2017, 9, 5, 15, 56, 22)
    #EndAvg=datetime.datetime(2017, 9, 5, 15, 56, 27)
    
    
    #ParticleStatsWv, ParticleTime,ParticleTimesWv= LoadGrayScaleStats(GSpath,ExpDate)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(0,0,ParticleStatsWv,1,405070)
    #D_KorolevCorr= KorolevCorrectedD(FilledArea, VoidArea,DiameterBG)
    #D_ReuterBakan= ReuterBakanDiameter(DiameterLevel0,AreaFraction0)

    #SaveGrayscaleOutput(GSpath,ExpDate,405070)
    
    Data_h5 = h5py.File(GSpath+'GrayscaleOutput.h5', 'r')              
    
    
    #NParticles=np.array(Data_h5['NParticles'])
    #h5f.create_dataset('ParticleTime', data=ParticleTime)
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleTime=[ExpDate + datetime.timedelta(seconds=int(ParticleTimesWv[x][0])) for x in range(len(ParticleTimesWv))]
    D_KorolevCorr=np.array(Data_h5['D_KorolevCorr'])
    #DiameterBG=np.array(Data_h5['DiameterBG'])
    Zd_fromRatios=np.array(Data_h5['Zd_fromRatios'])
    Data_h5.close()
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 16})

    # Zd filtering
    for i in range(4,9,1):
        PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,i,i,0,70,DLow_D0_Zd_Table)
        SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
        CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
        #plt.plot(PSD_SizeMid,Counts,label='Zd= '+str(i))
        plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Z$_{d}$ < '+str(i))        

    # No Zd filtering
    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,1E6,8.4,0,70,DLow_D0_Zd_Table)
    SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
    #plt.plot(PSD_SizeMid,Counts,label='c=8, no max Zd')
    plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Z$_{d}$ < 8.4') 

    # Particles per image filtering
    #PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,1E6,7,1,70)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    #plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='c=7, no max Zd, 1 particle per image')


    #NParticles=np.array(PartGrayStatsWv[:,0])

    # Reuter and Bakan diameter
    #PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_ReuterBakan,Zd_fromRatios,1E6,8,0,70,DLow_D0_Zd_Table)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    #plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Reuter & Bakan (1997)') 


    CoreCloudFile='core-cloud-phy_faam_20170905_v501_r0_c052_pcasp1.nc'
    CoreCloudPath='C:/Users/Admin TEMP/Documents/Clarify/C052/'
    CdpCalPath='C:/Users/Admin TEMP/Documents/Clarify/'
    CdpCalFile='IAN_CDPbins_formatted.csv'   
    CDP_BinCentre, CDP_dNdDp_L_avg= Avg_CDP_dndDp(StartAvg,EndAvg,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    plt.plot(CDP_BinCentre, CDP_dNdDp_L_avg,'o',label='CDP')

    #plt.plot()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Diameter, μm')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    #plt.ylabel('Counts')
    plt.legend()
    #plt.title(PltTitle)
    plt.savefig(GSpath+PltTitle,dpi=200)
    
#    return CDP_BinCentre, CDP_dNdDp_L_avg

#_______________________________________________________________________________________
    
## PICASSO flight C098  ********************************
    
def PICASSO_C098_grayscalePSD():
    DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D40_D0_Zd.hdf5'
    GSpath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/CIP15/'

    ExpDate= datetime.datetime(2018, 4, 24, 0, 0, 0)
    
    
    #PltTitle='C098_0939'   
    #StartAvg=datetime.datetime(2018, 4, 24, 9, 39, 35)    
    #EndAvg=datetime.datetime(2018, 4, 24, 9, 44, 11)
    
    #PltTitle='C098_run3'
    #StartAvg=datetime.datetime(2018, 4, 24, 9, 11, 53)    
    #EndAvg=datetime.datetime(2018, 4, 24, 9, 18, 31)
    #EndAvg=datetime.datetime(2018, 4, 24, 9, 25, 00)
    #HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/Run3/'
    #HaloName='C098_forSeb_091000_v2.h5'
    
    PltTitle='C098_run1'
    StartAvg=datetime.datetime(2018, 4, 24, 8, 48, 00)
    EndAvg=datetime.datetime(2018, 4, 24, 8, 55, 58)  
    HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/Run1/'
    HaloName='C098_084800_vNOunclassifed_particle_ID_class-ID_forSeb_v1_byPredict.h5'
    
    #PltTitle='C098_run6'
    #StartAvg=datetime.datetime(2018, 4, 24, 9, 39, 35)
    #EndAvg=datetime.datetime(2018, 4, 24, 9, 44, 11)  
    #HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/Run6/'
    ##HaloName='C098_forSeb_093935_OnlyTree084800.h5'   
    #HaloName='C098_forSeb_093935_Tree084800optimized.h5'
    
    
    
    
    #SaveGrayscaleOutput(GSpath,ExpDate,255075)
    
    Data_h5 = h5py.File(GSpath+'GrayscaleOutput.h5', 'r')              
     
    #NParticles=np.array(Data_h5['NParticles'])
    #ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    #h5f.create_dataset('ParticleTime', data=ParticleTime)
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleTime=[ExpDate + datetime.timedelta(seconds=int(ParticleTimesWv[x][0])) for x in range(len(ParticleTimesWv))]
    
    #AreaFraction0=np.array(Data_h5['AreaFraction0'])
    #AreaFraction1=np.array(Data_h5['AreaFraction1'])
    #AreaFraction2=np.array(Data_h5['AreaFraction2'])
    #AreaRatio_2_1=np.array(Data_h5['AreaRatio_2_1'])
    #AreaRatio_1_0=np.array(Data_h5['AreaRatio_1_0'])
    #DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
    #DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
    #DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
    #FilledArea=np.array(Data_h5['FilledArea'])
    #VoidArea=np.array(Data_h5['VoidArea'])
    DiameterBG=np.array(Data_h5['DiameterBG'])
    Zd_fromRatios=np.array(Data_h5['Zd_fromRatios'])
    #D_ReuterBakan=np.array(Data_h5['D_ReuterBakan'])
    #D_KorolevCorr=np.array(Data_h5['D_KorolevCorr'])
    
    #ParticleStatsWv, ParticleTime,ParticleTimesWv= LoadGrayScaleStats(GSpath,ExpDate)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(0,0,ParticleStatsWv,1,255075)
    #D_ReuterBakan= ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
    
    Data_h5.close()
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 16})
    #CoreCloudFile='core-cloud-phy_faam_20170905_v501_r0_c052_pcasp1.nc'
    #CoreCloudPath='C:/Users/Admin TEMP/Documents/Clarify/C052/'
    #CdpCalPath='C:/Users/Admin TEMP/Documents/Clarify/'
    #CdpCalFile='IAN_CDPbins_formatted.csv'

    #CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    
    #CDP_BinCentre, CDP_dNdDp_L_avg= Avg_CDP_dndDp(StartAvg,EndAvg,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)

    #plt.plot(CDP_BinCentre, CDP_dNdDp_L_avg)


    # Grayscale Plots
#    for i in range(4,9,1):
#        PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime,SVol,DoF=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,i,i,0)
#        SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
#        plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='Zd= '+str(i))

#    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime,SVol,DoF=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,1E6,8,0)
#    SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
#    plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='c=8, no max Zd')

    # Grayscale Plots
    for i in range(4,9,1):
        PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,i,i,0,40,DLow_D0_Zd_Table)
        #SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
        SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
        CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
        #plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='Zd= '+str(i))
        plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Z$_{d}$ < '+str(i))        

    #No Zd filtering    
    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,1E6,8,0,40,DLow_D0_Zd_Table)
    SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='All Z$_{d}$') 

    # Reuter and Bakan diameter
    #PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_ReuterBakan,Zd_fromRatios,1E6,8,0)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    #plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='c=8, no max Zd')
    #plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Reuter & Bakan (1997)') 

    BinMids, HaloHolo_dNdDp_L_avg, CountUncertainty=BinHaloHolo_v2(0.01, StartAvg, EndAvg, HaloName,HaloPath)
    plt.plot(BinMids,HaloHolo_dNdDp_L_avg,'o',label='HALOHolo')
    #plt.errorbar(BinMids, HaloHolo_dNdDp_L_avg,yerr=(CountUncertainty/2),'o',label='HALOHolo')  
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Diameter, μm')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    plt.legend()
    #plt.title(PltTitle)
    plt.savefig(GSpath+PltTitle,dpi=200)
    
    
    #return CDP_BinCentre, CDP_dNdDp_L_avg
  
#_______________________________________________________________________________________
    
## Clarify flight C052  ********************************
    
def ICED_B931_grayscalePSD():

    ExpDate= datetime.datetime(2015, 8, 19, 0, 0, 0)
    
    #GSpath='C:/Users/Admin TEMP/Documents/ICED/B931/20150819/OasisOut/'
    #PltTitle='B931_1356'
    #StartAvg=datetime.datetime(2015, 8, 19, 14, 11, 9)
    #EndAvg=datetime.datetime(2015, 8, 19, 14, 12, 30)
    
    
    GSpath='C:/Users/Admin TEMP/Documents/Clarify/C052/OasisOut/Imagefile_1CIP Grayscale_20170905140314/'
    PltTitle='C052_1608'
    StartAvg=datetime.datetime(2017, 9, 5, 16, 8, 40)
    EndAvg=datetime.datetime(2017, 9, 5, 16, 11, 10)
    
    
    ParticleStatsWv, ParticleTime,ParticleTimesWv= LoadGrayScaleStats(GSpath,ExpDate)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(0,0,ParticleStatsWv,1,255075)
    D_KorolevCorr= KorolevCorrectedD(FilledArea, VoidArea,DiameterBG)

    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 16})


    # Zd filtering
    for i in range(4,9,1):
        PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,i,i,0,70)
        #SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
        SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
        CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
        #plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='Zd= '+str(i))
        plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Z$_{d}$ < '+str(i))        

    # No Zd filtering
    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,1E6,8,0,70)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
    SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    #plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='c=8, no max Zd')
    plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='c=8, no max Z$_{d}$') 

    # Particles per image filtering
    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD(ParticleStatsWv, ParticleTime,ParticleTimesWv,D_KorolevCorr,Zd_fromRatios,1E6,7,1,70)
    #SizeWidth, PSD_avg_dNdDp, PSD_avg= Avg_dndDp(StartAvg,EndAvg,PSD_L, PSD_MidTime, PSD_SizeMid)
    SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    #plt.plot(PSD_SizeMid,PSD_avg_dNdDp,label='c=8, no max Zd')
    plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='c=7, no max Zd, 1 particle per image')


    #NParticles=np.array(PartGrayStatsWv[:,0])



    CoreCloudFile='core-cloud-phy_faam_20150819_v501_r1_b931.nc'
    CoreCloudPath='C:/Users/Admin TEMP/Documents/ICED/B931/'
    CdpCalPath='C:/Users/Admin TEMP/Documents/CDP/'
    CdpCalFile='master_water_diameters.csv'
#    CDP_time_mid, CDP_cm3, CDP_dNdDp, CDP_gcm3, CDP_BinCentre= LoadCoreCloud(CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)    
    CDP_BinCentre, CDP_dNdDp_L_avg= Avg_CDP_dndDp(StartAvg,EndAvg,CoreCloudPath,CoreCloudFile,CdpCalPath,CdpCalFile)
    plt.plot(CDP_BinCentre, CDP_dNdDp_L_avg,'o',label='CDP')

    #plt.plot()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Diameter, μm')
    plt.ylabel('dN/dDp')
    plt.legend()
    #plt.title(PltTitle)
    plt.savefig(GSpath+PltTitle,dpi=200)
    
    return CDP_BinCentre, CDP_dNdDp_L_avg


#_______________________________________________________________________________________
    


def SaveGrayscaleOutput(GSpath,ExpDate,ThresholdFlag):
    FileName= 'GrayscaleOutput.h5'
    try:
        os.remove(GSpath+FileName)
    except OSError:
        pass
    
    ParticleStatsWv, ParticleTime,ParticleTimesWv= LoadGrayScaleStats(GSpath,ExpDate)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(0,0,ParticleStatsWv,1,ThresholdFlag)
    D_KorolevCorr= KorolevCorrectedD(FilledArea, VoidArea,DiameterBG)
    D_ReuterBakan= ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
        
    #StartDate=datetime.datetime.combine(ParticleTime[0].date(),datetime.time(second=0))
    #IAT= GetIAT_vector(ParticleTime,ParticleTimesWv)    
    #ParticleTime_sinceMidnight= ParticleTimesWv[:,0]+(ParticleTimesWv[:,1]/1000000000)
        
    
    
    h5f = h5py.File(GSpath+FileName, 'w')
    
    #h5f.create_dataset('ParticleStatsWv', data=ParticleStatsWv)
    #h5f.create_dataset('ParticleTime', data=ParticleTime)
    h5f.create_dataset('ParticleTimesWv', data=ParticleTimesWv)
    #h5f.create_dataset('StartDate', data=StartDate)
    #h5f.create_dataset('IAT', data=IAT)
    #h5f.create_dataset('ParticleTime_sinceMidnight', data=ParticleTime_sinceMidnight)
    #h5f.create_dataset('AreaFraction0', data=AreaFraction0)
    #h5f.create_dataset('AreaFraction1', data=AreaFraction1)
    #h5f.create_dataset('AreaFraction2', data=AreaFraction2)
    #h5f.create_dataset('AreaRatio_2_1', data=AreaRatio_2_1)
    #h5f.create_dataset('AreaRatio_1_0', data=AreaRatio_1_0)
    #h5f.create_dataset('DiameterLevel0', data=DiameterLevel0)
    #h5f.create_dataset('DiameterLevel1', data=DiameterLevel1)
    #h5f.create_dataset('DiameterLevel2', data=DiameterLevel2)
    #h5f.create_dataset('FilledArea', data=FilledArea)
    #h5f.create_dataset('VoidArea', data=VoidArea)
    h5f.create_dataset('DiameterBG', data=DiameterBG)
    h5f.create_dataset('Zd_fromRatios', data=Zd_fromRatios)
    h5f.create_dataset('D_ReuterBakan', data=D_ReuterBakan)
    h5f.create_dataset('D_KorolevCorr', data=D_KorolevCorr)
    
    h5f.close()
    

#_______________________________________________________________________________________

# Only makes psd between StartAvg and EndAvg


def ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,Diameter,Zd_fromRatios,Zd_max,DOF_c,MaxParticlesPerImage, ArmSep, DLow_D0_Zd_Table):
    
    IATmin=1E-5
    #IATmin=1E-6
    #IATmin=0
    #Zd_max=6 # maximum Zd to include in PSD
    #DOF_c= 8 # c value for DOF calculation, normally should be same as Zd_max
    TAS=100
    Lambda = 0.658	# in um 
    Lambda_mm=  Lambda / 1000
    ProbeRes= 15
    #ArmSep = 70
    ArrayElements = 64
    
    #NParticles=np.array(ParticleStatsWv[:,0])
    
    StartDate=datetime.datetime.combine(ParticleTime[0].date(),datetime.time(second=0))
    IAT= GetIAT_vector(ParticleTime,ParticleTimesWv)    
    ParticleTime_sinceMidnight= ParticleTimesWv[:,0]+(ParticleTimesWv[:,1]/1000000000)
    
    
    Startime=int(datetime_2_s_sincemidnight(StartAvg))
    Endtime=int(datetime_2_s_sincemidnight(EndAvg))
    
    Ns=Endtime-Startime
    
    TimeBins_s=np.linspace(Startime+0.5,Endtime+0.5,num=Ns+1)    
    
    SizeBins=np.linspace((ProbeRes*1.5),(ProbeRes*ArrayElements+ProbeRes*0.5),num=ArrayElements)
    PSD_SizeMid=SizeBins-ProbeRes*0.5
    
    Time_idx = np.digitize(ParticleTime_sinceMidnight, TimeBins_s)
    Size_idx = np.digitize(Diameter, SizeBins)
    
    # Calculate concentration of individual particle
    DLow_D0_model, Zd_model=Get_DLow_D0_model(DLow_D0_Zd_Table)
    ConcentrationL= np.zeros(len(ParticleTime_sinceMidnight))
    for i in range(len(ParticleTime_sinceMidnight)):
        if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime):
            #ConcentrationL[i]=1/Allin_SVol_lawson(Diameter[i],DiameterLevel0x[i],Lambda,DOF_c,TAS,ProbeRes,ArmSep, ArrayElements)
            ConcentrationL[i]=1/sVol_EAW_INT(TAS,Diameter[i],DOF_c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model)
    
    # Filter and bin data
    Counts_PSD = np.zeros((len(TimeBins_s), len(SizeBins)))
    PSD_L = np.zeros((len(TimeBins_s), len(SizeBins))) 
    
    if MaxParticlesPerImage==0 : 
        for i in range(len(ParticleTime_sinceMidnight)-1) : 
            
            if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (Diameter[i] > (0.5*ProbeRes)) & (Diameter[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) & (Zd_fromRatios[i] <Zd_max): 
                #if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (DiameterLevel1[i] > (0.5*ProbeRes)) & (DiameterLevel1[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) :
                Counts_PSD[Time_idx[i],Size_idx[i]]+=1
                PSD_L[Time_idx[i],Size_idx[i]]+= ConcentrationL[i]
    else :
        print ('Need to include max particles per image')
    #else : 
    #    for i in range(len(Diameter)-1) : 
    #        if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (Diameter[i] > (0.5*ProbeRes)) & (Diameter[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) & (Zd_fromRatios[i] <Zd_max) & (NParticles[i] <= MaxParticlesPerImage): 
    #            #if (ParticleTime_sinceMidnight[i] > Startime) & (ParticleTime_sinceMidnight[i] < Endtime) & (DiameterLevel1[i] > (0.5*ProbeRes)) & (DiameterLevel1[i] < (ProbeRes*ArrayElements+ProbeRes*0.5)) & (IAT[i] > IATmin) :
    #            Counts_PSD[Time_idx[i],Size_idx[i]]+=1    
    #            PSD_L[Time_idx[i],Size_idx[i]]+= ConcentrationL[i]
    
    PSD_MidTime=[StartDate + datetime.timedelta(seconds=int(TimeBins_s[x]-0.5)) for x in range(len(TimeBins_s))]
   
    return PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime#, SVol,DoF
    
    #Zd = GetZd(AreaRatio_1_0,AreaFraction2)
    #return Zd_fromRatios, IAT


#_______________________________________________________________________________________
    
## PICASSO flight C098  ********************************
    
def PICASSO_C078_grayscalePSD():
    DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D40_D0_Zd.hdf5'
    GSpath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C078/CIP15/'
    HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C078/CIP15/'
    HaloName='C078_forSeb_v1.h5'

    ExpDate= datetime.datetime(2018, 2, 7, 0, 0, 0)
    
    
    PltTitle='C078_1602'
    StartAvg=datetime.datetime(2018, 2, 7, 16, 2, 0)
    EndAvg=datetime.datetime(2018, 2, 7, 16, 10, 0)  
    #HaloPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C098/Run6/'   
    #HaloName='C098_forSeb_093935_Tree084800optimized.h5'
    
    #SaveGrayscaleOutput(GSpath,ExpDate,255075)
    
    Data_h5 = h5py.File(GSpath+'GrayscaleOutput.h5', 'r')              
       
    #NParticles=np.array(Data_h5['NParticles'])
    #h5f.create_dataset('ParticleTime', data=ParticleTime)
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    ParticleTime=[ExpDate + datetime.timedelta(seconds=int(ParticleTimesWv[x][0])) for x in range(len(ParticleTimesWv))]
    
    DiameterBG=np.array(Data_h5['DiameterBG'])
    Zd_fromRatios=np.array(Data_h5['Zd_fromRatios'])
    
    Data_h5.close()
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 16})

    # Grayscale Plots
    for i in range(4,9,1):
        PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,i,i,0,40,DLow_D0_Zd_Table)
        SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
        CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
        plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='Z$_{d}$ < '+str(i))        

    #No Zd filtering    
    #ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,Diameter,Zd_fromRatios,Zd_max,DOF_c,MaxParticlesPerImage, ArmSep, DLow_D0_Zd_Table)
    PSD_L, Counts_PSD, PSD_SizeMid, PSD_MidTime=ParticlesStats2PSD_v2(StartAvg, EndAvg, ParticleTime,ParticleTimesWv,DiameterBG,Zd_fromRatios,1E6,9,0,40,DLow_D0_Zd_Table)
    SizeWidth, PSD_avg_dNdDp, PSD_avg, Counts=Avg_dndDp_CountStats(StartAvg,EndAvg,PSD_L,Counts_PSD, PSD_MidTime, PSD_SizeMid)
    CountUncertainty=  PSD_avg_dNdDp / np.sqrt(Counts)
    plt.errorbar(PSD_SizeMid,PSD_avg_dNdDp,yerr=(CountUncertainty/2), label='All Z$_{d}$') 

    
    BinMids, HaloHolo_dNdDp_L_avg, CountUncertainty=BinHaloHolo_v2(0.01, StartAvg, EndAvg, HaloName,HaloPath)
    plt.plot(BinMids,HaloHolo_dNdDp_L_avg,'o',label='HALOHolo')
    #plt.errorbar(BinMids, HaloHolo_dNdDp_L_avg,'o',yerr=(CountUncertainty/2),label='HALOHolo')  
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Diameter, μm')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    plt.legend()
    #plt.title(PltTitle)
    plt.savefig(GSpath+PltTitle,dpi=200)
    
    
    #return CDP_BinCentre, CDP_dNdDp_L_avg
  
#_______________________________________________________________________________________

def CompareSA_curves():
    
    Lambda = 0.658	# in um 
    Lambda_mm=  Lambda / 1000
    D25_D0_model, Zd_model_25= Get_DLow_D0_model('C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5')   
    D50_D0_model, Zd_model_50= Get_DLow_D0_model('C:/Users/Admin TEMP/Documents/DropletGun/D50_D0_Zd.hdf5') 
    
    D_array=np.linspace(15,960,64)
    SA_25=np.zeros(64)
    SA_50=np.zeros(64)
    
    # 25 
    for i in range(len(D_array)):
        D0=D_array[i]
        SA_25[i]=SA_EAW_INT(D0,7,64,15,Lambda_mm, 70,D25_D0_model,Zd_model_25)
        SA_50[i]=SA_EAW_INT(D0,7,64,15,Lambda_mm, 70,D50_D0_model,Zd_model_50)
        
        
    plt.plot(D_array, SA_25)
    plt.plot(D_array, SA_50)

    plt.xscale('log')
    
    
    return D_array, SA_25, SA_50