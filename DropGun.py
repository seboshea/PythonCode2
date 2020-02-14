# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 16:20:46 2018

@author: Admin TEMP
"""

import numpy as np
import matplotlib.pyplot as plt  
import datetime
import h5py
from matplotlib.dates import DateFormatter
import bisect
import csv
from MyFunctions import is_number,LoadCoreCloud,loadFAAMCore,ChangeTimeBaseAvg, BinLinearV2, KorolevCorrectedD
from MyFunctions import Level2_D_Correction_vector, ReuterBakanDiameter, DimensionlessDistance
import os
from numpy import genfromtxt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math
import pandas as pd
from GrayscaleFunctions import LoadGrayScaleStats, GrayscaleRatios, sVol_EAW_INT,Get_DLow_D0_model
from scipy.stats.stats import pearsonr
import matplotlib.ticker as plticker

#ParticleTimesWv=[]
#PartGrayStatsWv=[]


#_______________________________________________________________________________________   
 
#Collate stats from experiments
   
def DropGun_Stats_All():    

    DropGun_Stats={}
    StageFile='StagePositions_exp001.csv'

    
    #60um_25_50_75
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75/20180327125851/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_25_50_75_run2'


    #60um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/20180327143000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['60um_25_50_75_run2']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #60um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60/20180327133536/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_40_50_60'

    #60um_40_50_60_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60_run2/20180327140438/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['60um_40_50_60_run2']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #90um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_40_50_60/20180326143716/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 26, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['90um_40_50_60']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #90um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/20180326150258/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 26, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['90um_25_50_75_run2']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #120um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/20180327100000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['120um_25_50_75_run2']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #120um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_40_50_60/20180327092820/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    #AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #DropGun_Stats['120um_40_50_60']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)



    #********************************* September 2018 ****************************************************
    #60um_255075
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/20180927111105/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0) 
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    DropGun_Stats['60um_25_50_75']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)
    
    #90um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/20180927134617/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    DropGun_Stats['90um_25_50_75']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    #90um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/Both/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    DropGun_Stats['120um_25_50_75']=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)


    #AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg, VoidRatioAvg

    return DropGun_Stats



#_______________________________________________________________________________________

# pLot Particle_X vs Particle_Y with AreaFraction and Diameter as colourscale


def PlotXY_D_AF(StagePath, PltTitle, Particle_X, Particle_Y, AreaFraction,Diameter):
    
    cm = plt.cm.get_cmap('RdYlBu')

    fig=plt.figure(figsize=(10,10))
    
    fig.add_subplot(2, 1, 1)
    plt.title(PltTitle)
    #plt.gca().set_ylim([0, 8000]) 
    AFC = plt.scatter(Particle_X, Particle_Y, c=AreaFraction,vmin=0, vmax=1,cmap=cm)
    #AFC.ylabel('Stage Y, ,m')
    cbar=plt.colorbar(AFC)
    cbar.set_label('Level 0 / Total area')
    
    #plt.tick_params(axis='245° radial', labelcolor='tab:red')
    fig.add_subplot(2, 1, 2)
    DC = plt.scatter(Particle_X, Particle_Y, c=Diameter,vmin=0, vmax=100,cmap=cm)
    cbar2=plt.colorbar(DC)
    cbar2.set_label('Diameter, μm')
    
    
    #plt.show()
    plt.savefig(StagePath+'XY_D_AF'+PltTitle+'.png',dpi=200)
    plt.close(fig)


#_______________________________________________________________________________________

# Checking the impact on TAS on the DMT CIP DOF

def CIP_TAS_DOF():
    
    Dof_D, Dof_Y_1=DOF_curve(1,22)
    Dof_D, Dof_Y_2=DOF_curve(2,22)
    Dof_D, Dof_Y_4=DOF_curve(4,22)
    Dof_D, Dof_Y_10=DOF_curve(10,22)
    Dof_D, Dof_Y_20=DOF_curve(20,22)
    Dof_D, Dof_Y_100=DOF_curve(100,22)
    
    plt.plot(Dof_D, Dof_Y_1,'o',markersize=1,color='blue') 
    plt.plot(Dof_D, Dof_Y_2,'o',markersize=1,color='black')
    plt.plot(Dof_D, Dof_Y_4,'o',markersize=1,color='black')
    plt.plot(Dof_D, Dof_Y_10,'o',markersize=1,color='black')
    plt.plot(Dof_D, Dof_Y_20,'o',markersize=1,color='black')
    plt.plot(Dof_D, Dof_Y_100,'o',markersize=1,color='red')
    
    
#_______________________________________________________________________________________
# Return the DOF for CIP @ Diameter
# Lamda is hard coded    
# Diameter in um

def CIP_DOF(Diameter,TAS,ArmSep):
   #Diameter=np.linspace(0,50,51)
	
   Taur=0.1						# deatime of CIP in ms
   Lambda = 0.000658			# laser wavelength in mm (658nm)	// wavelength of laser in mm
   f=100*np.minimum(0.5/(1-math.exp(-Diameter/(TAS*Taur))),1)
   z=30.46-(0.628*f)+(0.003246*f*f)
   radius = Diameter/2000						# radius in mm
   DoF = min(z*radius*radius/Lambda,ArmSep)		# dof in mm		
   return DoF

#_______________________________________________________________________________________

# Curve DOF vs diameter


def DOF_curve(TAS,CentreOfDOF):
    Diameter=np.linspace(1,150,150)
    DOF=np.empty(len(Diameter), dtype=float)
    #TAS=4
    for i in range(len(Diameter)):
        DOF[i]=CIP_DOF(Diameter[i],TAS,70)
    
    Dof_Y_plus=CentreOfDOF+(DOF/2)
    Dof_Y_minus=CentreOfDOF-(DOF/2)
    
    Dof_Y= np.append(Dof_Y_plus,Dof_Y_minus)
    Dof_D= np.append(Diameter,Diameter)
    
    #plt.plot(Dof_D, Dof_Y,'o')
    
    
    return Dof_D, Dof_Y
#_______________________________________________________________________________________

    #c=6
    
def DOF_lawson_curve(c):
    Diameter=np.linspace(1,150,150) # um
    DOF=np.empty(len(Diameter), dtype=float)
    
    Lambda = 0.658	# in um   
    #c=6
    
    for i in range(len(Diameter)):
        DOF[i]=Dof_lawson(Diameter[i],Lambda,c)

    Dof_Y_plus=22+(DOF/2)
    Dof_Y_minus=22-(DOF/2)
    
    Dof_lawson_Y= np.append(Dof_Y_plus,Dof_Y_minus)
    Dof_lawson_D= np.append(Diameter,Diameter)
    
    #plt.plot(Dof_D, Dof_Y,'o')
    
    
    return Dof_lawson_D, Dof_lawson_Y


#_______________________________________________________________________________________
## From lawson et al 2006  equation 1
# Lambda in um
#Diameter in um    

def Dof_lawson(Diameter,Lambda,c):
    radius=Diameter/2 # um
    DoF=(2*c*radius**2)/Lambda ## 2* because there is a +/- in the equation
    DoF/=1000 # return in mm
    return DoF

#_______________________________________________________________________________________

 #AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre


def DropGun_Stats_Plot(DropGun_Stats):
    

    ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['60um_40_50_60_run2','90um_40_50_60','120um_40_50_60'] 
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/'
    PltTitle='CIP Model'
    fig=plt.figure(figsize=(10,10))
    i=0
    for i in range(len(ExpList)):
        D_Level1_avg=DropGun_Stats[ExpList[i]][4][:]
        BinCentre=DropGun_Stats[ExpList[i]][6][:]
        AF_Level0_avg=DropGun_Stats[ExpList[i]][2][:] 
        #plt.plot(BinCentre, D_Level1_avg)
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(BinCentre, D_Level1_avg, c=AF_Level0_avg,vmin=0, vmax=1,cmap=cm) 
    Dof_D,Dof_Y=DOF_curve(4,22)
    plt.plot(Dof_Y,Dof_D,'o',markersize=1,color='black')    
    cbar=plt.colorbar(sc)
    cbar.set_label('Level 0 / Total area')
    plt.xlabel('Stage Y, mm')
    plt.xlim([-5,55])
    plt.ylabel('Diameter >=1, μm')
    plt.title(PltTitle)
    #plt.show()
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    plt.close(fig)
   
#_______________________________________________________________________________________

 #AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre


def DropGun_Stats_korolev_Plot(DropGun_Stats):
    
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['60um_40_50_60_run2','90um_40_50_60','120um_40_50_60'] 
    #FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/'
    
    
    PltTitle='120_90_60um_GS_25_50_75_korolev'
    fig=plt.figure(figsize=(10,10))
    i=0
    for i in range(len(ExpList)):
        D_Level1_avg=DropGun_Stats[ExpList[i]][7][:]
        BinCentre=DropGun_Stats[ExpList[i]][6][:]
        AF_Level0_avg=DropGun_Stats[ExpList[i]][2][:] 
        #plt.plot(BinCentre, D_Level1_avg)
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(BinCentre, D_Level1_avg, c=AF_Level0_avg,vmin=0, vmax=1,cmap=cm) 

    # DoF Lines
    Dof_D,Dof_Y=DOF_curve(4,22)
    plt.plot(Dof_Y,Dof_D,'o',markersize=1,color='black', label='OASIS CIP DOF')

    Dof_lawson_c6_D, Dof_lawson_c6_Y=DOF_lawson_curve(6)
    plt.plot(Dof_lawson_c6_Y,Dof_lawson_c6_D,'o',markersize=1,color='black', label='Lawson et al, c=6')
    
    Dof_lawson_c8_D, Dof_lawson_c8_Y=DOF_lawson_curve(8)
    plt.plot(Dof_lawson_c8_Y,Dof_lawson_c8_D,'o',markersize=1,color='black', label='Lawson et al, c=8')
    
    cbar=plt.colorbar(sc)
    cbar.set_label('Level 0 / Total area')
    plt.xlabel('Stage Y, mm')
    plt.xlim([-5,55])
    plt.ylabel('Diameter (Korolev corrected), μm')
    plt.title(PltTitle)
    plt.legend()
    #plt.show()
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    plt.close(fig)
   
#_______________________________________________________________________________________
# Plot diameter and area fraction for given drop gun size and grayscale setting


def PlotDiameterAF():

    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/Autumn2017/60um_25_50_75/20171018/20171018091312/OASIS_Out/'
    #ExpDate= datetime.datetime(2017, 10, 18, 0, 0, 0)
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/Autumn2017/60um_25_50_75/'

    #***************** March 2018 ********************

    #60um_25_50_75
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75/20180327125851/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_25_50_75_run2'


    #60um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/20180327143000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_25_50_75_run2'
    #ActualDiameter=65

    #60um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60/20180327133536/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_40_50_60'

    #60um_40_50_60_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60_run2/20180327140438/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_40_50_60_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_40_50_60_run2'
    #ActualDiameter=65

    #90um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_40_50_60/20180326143716/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 26, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='90um_40_50_60'
    #ActualDiameter=90
    
    #90um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/20180326150258/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 26, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='90um_25_50_75_run2'
    #ActualDiameter=85

    #120um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/20180327100000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP   
    #PltTitle='120um_25_50_75_run2'
    #ActualDiameter=105
    
    #120um_40_50_60
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_40_50_60/20180327092820/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_40_50_60/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='120um_40_50_60'
    #ActualDiameter=105

    
    #********************************* September 2018 ****************************************************
    #60um_255075
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/20180927111105/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/'
    #ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0) 
    #ActualDiameter = 55 
    
    #90um_25_50_75
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/20180927134617/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/'
    #ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    #ActualDiameter = 80
    
    #120um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/Both/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    PltTitle='120um_25_50_75'
    ActualDiameter = 90


    StageFile='StagePositions_exp001.csv'
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    #PlotNumberOfParticles(FilledArea,VoidArea,Particle_X,DiameterLevel0,PartGrayStatsWv,StagePath,'Nparticles'+PltTitle)
    #Figurename=StagePath+'DiameterPosition'
    #PlotGrayscaleDiameter(DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X,Figurename,PltTitle)
    #Figurename=StagePath+'AreaFractionPosition'
    #PlotAreaFraction(AreaFraction0,AreaFraction1,AreaFraction2,Particle_X,Figurename,PltTitle)
    #Figurename=StagePath+'DiameterTimeSeries'
    #PlotGSTimeSeries(Particle_X, Particle_Y, DiameterLevel0,DiameterLevel1,DiameterLevel2, ParticleTime,Figurename,PltTitle)
    #PlotXY_D_AF(StagePath, PltTitle,Particle_X, Particle_Y, AreaFraction0,DiameterLevel2)
    
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    #D_2stage=Level2_D_Correction_vector(D_KorolevCorr, AreaFraction2)
    D_ReuterBakan=ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
    
    D_ReuterBakan_median, D_ReuterBakan_25,D_ReuterBakan_75, \
        D_ReuterBakan_c6_median,D_ReuterBakan_c6_25,D_ReuterBakan_c6_75, \
        D_Korolev_median,D_Korolev_25,D_Korolev_75, \
        D_Korolev_c6_median, D_Korolev_c6_25, D_Korolev_c6_75, \
        D_median, D_25, D_75, \
        D_c6_median, D_c6_25, D_c6_75=PercentilesKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,ActualDiameter,29)
    
    
    print(D_ReuterBakan_c6_median)
    print(D_ReuterBakan_c6_75-D_ReuterBakan_c6_25)
    print(D_Korolev_c6_median)
    print(D_Korolev_c6_75-D_Korolev_25)
    print(D_c6_median)
    print(D_c6_75-D_c6_25)

    
    #Figurename=StagePath+'D_2stage_'+PltTitle
    #Plot2StageCorrD(D_2stage, D_KorolevCorr,DiameterLevel1, Particle_X,Figurename,PltTitle,ActualDiameter,22)
    
    
    #D_BG_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterBG)
    #Figurename=StagePath+'KorolevD_'+PltTitle
    #PlotKorCorrD(D_KorolevCorr,DiameterLevel1,Particle_X,Figurename,PltTitle,ActualDiameter)
    #Figurename=StagePath+'BG_particles_'+PltTitle
    #Plot_BG_diameter(D_KorolevCorr,DiameterLevel1, D_BG_KorolevCorr, DiameterBG, Particle_X,Figurename,PltTitle,ActualDiameter,22)

    #Zd=DimensionlessDistance(Particle_X,ActualDiameter,22)
    #Figurename=StagePath+'AreaFraction_Zd'
    #PlotAreaFractionZd(AreaFraction0,AreaFraction1,AreaFraction2,Zd,Figurename,PltTitle)
    
#_______________________________________________________________________________________

#VoidRatio=VoidAreaRatio(FilledArea,VoidArea,Particle_Y)

#def VoidAreaRatio(FilledArea,VoidArea,Particle_Y): 
    
#    AreaRatio = VoidArea / FilledArea
#    AreaRatio_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaRatio,Particle_Y,0,50,50)
#    mask = ~np.isnan(AreaRatio_bin) # remove nans
#    filtered_data = [d[m] for d, m in zip(AreaRatio_bin.T, mask.T)]
#    VoidRatio = [np.nanmedian(x) for x in filtered_data]
#    return VoidRatio







#_______________________________________________________________________________________
# PLot Counts vs position


def PlotNumberOfParticles(FilledArea,VoidArea,Particle_Y,DiameterLevel0,PartGrayStatsWv,SavePath,Title):
    
    NParticles=np.array(PartGrayStatsWv[:,0])  
    AreaRatio = VoidArea / FilledArea
    
    
    for x in range(len(NParticles)-1):
        if (DiameterLevel0[x]!=DiameterLevel0[x]):
            NParticles[x]=np.nan          
        
    #Number of particles **************************
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2,1,1) 
    NParticles_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(NParticles,Particle_Y,0,50,50)     
    mask = ~np.isnan(NParticles_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(NParticles_bin.T, mask.T)]
    means = [np.nanmean(x) for x in filtered_data]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Number of particles')
    plt.xlabel('Stage Y, mm')
   
    #AreaRatio **************************
    plt.subplot(2,1,2,sharex=ax1)
    AreaRatio_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaRatio,Particle_Y,0,50,50)
    mask = ~np.isnan(AreaRatio_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaRatio_bin.T, mask.T)]
    means = [np.nanmean(x) for x in filtered_data]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Void area / Filled area (Excludes NANs)')
    plt.xlabel('Stage Y, mm')
 
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    plt.savefig(SavePath+Title,dpi=200)
    plt.close(fig)


#_______________________________________________________________________________________

# Load grayscale data and stage positions


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

def LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile):

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
            if 'ParticleTimesWv' not in dir():
                ParticleTimesWv=TmpTimesWv
                #PartGrayStatsWv=TmpStats
                ParticleStatsWv=TmpStats
            else:    
                ParticleTimesWv= np.append(ParticleTimesWv,TmpTimesWv,axis=0)
                #PartGrayStatsWv= np.append(PartGrayStatsWv,TmpStats,axis=0)
                ParticleStatsWv= np.append(ParticleStatsWv,TmpStats,axis=0)
            
        
    ParticleTime=[ExpDate + datetime.timedelta(seconds=int(ParticleTimesWv[x][0])) for x in range(len(ParticleTimesWv))]          

    #Load stage position
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/60um_25_50_75/'
    #StageFile='StagePositions_exp001.csv'  
    StageY = np.loadtxt(StagePath+StageFile,skiprows=2,usecols=range(2,3),delimiter=',')      
    StageX = np.loadtxt(StagePath+StageFile,skiprows=2,usecols=range(1,2),delimiter=',') 
    #str2date = lambda x: datetime.datetime.strptime(x.decode("utf-8"), '%d/%m/%Y %H:%M:%S')
    str2date = lambda x: datetime.datetime.strptime(x.decode("utf-8"), '%d/%m/%y %H:%M:%S.%f')
    StageTimeTmp= np.genfromtxt(StagePath+StageFile,skip_header=2,usecols=range(0),delimiter=',',converters= {0: str2date})
    StageTime=[StageTimeTmp[x][0] for x in range(len(StageTimeTmp))]
    
    Match=np.zeros(len(ParticleTime))*np.nan
    Particle_X=np.zeros(len(ParticleTime))*np.nan
    Particle_Y=np.zeros(len(ParticleTime))*np.nan

    for i in range(len(ParticleTime)):
        Match[i]=bisect.bisect_left(StageTime, ParticleTime[i])
        if Match[i] in range(0,len(StageX)-1):
            #if np.absolute(((ParticleTime[i]-StageTime[int(Match[i])])).seconds)<2 :
            if np.absolute(((ParticleTime[i]-StageTime[int(Match[i])])).total_seconds())<2 :
                Particle_X[i]=StageX[int(Match[i])]
                Particle_Y[i]=StageY[int(Match[i])]
            else: 
                Match[i]=np.nan	
        else: 
            Match[i]=np.nan
    #plt.plot(StageTime,StageY,"o")
    #plt.plot(ParticleTime,Particle_Y,"o") 

    #ParticleStage_Y=ChangeTimeBaseAvg(StageTime,StageY,ParticleTime)

    #plt.plt(StageTime,StageY )
    return Particle_X, Particle_Y, ParticleStatsWv, ParticleTime,ParticleTimesWv



#_______________________________________________________________________________________  


def PlotGSTimeSeries(Particle_X, Particle_Y, DiameterLevel0,DiameterLevel1,DiameterLevel2, ParticleTime,Figurename,PltTitle):

    fig=plt.figure(figsize=(10,10))
    
    formatter = DateFormatter('%H:%M')
    ax1=fig.add_subplot(2, 1, 1)
    plt.title(PltTitle)
    
    #plt.gca().set_ylim([0, 8000]) 
    plt.ylabel('Stage X, mm', color='tab:red')  # we already handled the x-label with ax1
    plt.plot(ParticleTime, Particle_X,"o",color='tab:red')
    #plt.tick_params(axis='245° radial', labelcolor='tab:red')
    
    ax2 = ax1.twinx()
    ax2.plot(ParticleTime, Particle_Y,"o",color='tab:blue')
    ax2.set_ylabel('Stage, Y',color='tab:blue')
    
    plt.subplot(2, 1, 2, sharex=ax1)
    plt.plot(ParticleTime, DiameterLevel0,"o")
    plt.plot(ParticleTime, DiameterLevel1,"o")
    plt.plot(ParticleTime, DiameterLevel2,"o")
    #plt.gca().set_ylim([-45, 10])
    plt.ylabel('Diameter, μm')
    
    plt.gcf().axes[0].xaxis.set_major_formatter(formatter)
    plt.xlabel('time')
    #StartPt=FindFirstGreaterThan(ALT_GIN,150)
    #plt.gca().set_xlim([Time_Core[StartPt],Time_Core[-1]]) #C082
    #plt.gca().set_xlim([Time_Core[StartPt], datetime.datetime(2018, 2, 14,20,30,00)]) #C082

    plt.tight_layout()       
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)

    


#_______________________________________________________________________________________   
    


def Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_Y):
    
    StartBin=0
    EndBin=100
    NumBins=100
    
    #Area fraction 2
    AreaFraction2_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction2,Particle_Y,StartBin,EndBin,NumBins)
    BinCentre=BinsMid
    
    mask = ~np.isnan(AreaFraction2_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction2_bin.T, mask.T)]
    AF_Level2_avg = [np.nanmedian(x) for x in filtered_data]

    #Area fractoin 1
    AreaFraction1_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction1,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(AreaFraction1_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction1_bin.T, mask.T)]
    AF_Level1_avg = [np.nanmedian(x) for x in filtered_data]
    
    #Area fraction 0
    AreaFraction0_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction0,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(AreaFraction0_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction0_bin.T, mask.T)]
    AF_Level0_avg = [np.nanmedian(x) for x in filtered_data]
    
    #Diameter 2
    DiameterLevel2_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel2,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(DiameterLevel2_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel2_bin.T, mask.T)]
    D_Level2_avg = [np.nanmedian(x) for x in filtered_data]
    
    #Diameter 1
    DiameterLevel1_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel1,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(DiameterLevel1_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel1_bin.T, mask.T)]
    #D_Level1_avg = [np.nanmedian(x) for x in filtered_data]    
    D_Level1_avg = [np.nanmean(x) for x in filtered_data]   
    #D_Level1_avg = [np.nanpercentile(x,25) for x in filtered_data]
    
    
    
    
    #Diameter 0
    DiameterLevel0_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel0,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(DiameterLevel0_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel0_bin.T, mask.T)]
    D_Level0_avg = [np.nanmedian(x) for x in filtered_data]
    
    #Korolev corrected diameter
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,StartBin,EndBin,NumBins)
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    D_KorolevCorr_avg = [np.nanmedian(x) for x in filtered_data]
    
    AreaRatio = VoidArea / FilledArea
    AreaRatio_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaRatio,Particle_Y,0,100,100)
    mask = ~np.isnan(AreaRatio_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaRatio_bin.T, mask.T)]
    VoidRatioAvg = [np.nanmedian(x) for x in filtered_data]
#    return VoidRatio
    
    return AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg, VoidRatioAvg
    
    
#_______________________________________________________________________________________    


#Plot stage positions versus grayscale fraction box plots 

def PlotAreaFraction(AreaFraction0,AreaFraction1,AreaFraction2,Particle_Y,Figurename,PltTitle):

    #Area fraction 2
    AreaFraction2_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction2,Particle_Y,0,50,50)

    mask = ~np.isnan(AreaFraction2_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction2_bin.T, mask.T)]

    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(3, 1, 3)
    
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.gca().set_ylim([-0.05, 1.05])
    plt.gca().set_ylim([-0.05, 1.05])
    
    plt.ylabel('Level 2 / Total area')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_xlim([-0.05, 50.05])

    #Area fraction 1
    AreaFraction1_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction1,Particle_Y,0,50,50)
    mask = ~np.isnan(AreaFraction1_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction1_bin.T, mask.T)]
    plt.subplot(3, 1, 2, sharex=ax1)
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Level 1 / Total area')
    plt.gca().set_ylim([-0.05, 1.05])

    #Area fraction 0
    AreaFraction0_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction0,Particle_Y,0,50,50)
    mask = ~np.isnan(AreaFraction0_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction0_bin.T, mask.T)]
    plt.subplot(3, 1, 1, sharex=ax1)
    plt.title(PltTitle)
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Level 0 / Total area ')
    plt.gca().set_ylim([-0.05, 1.05])
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    #Figurename='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/60um_25_50_75/AreaFraction.png'
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)
    
    
#_______________________________________________________________________________________    


def PlotGrayscaleDiameter(DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_Y,Figurename,PltTitle):

    #Diameter 2
    DiameterLevel2_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel2,Particle_Y,0,50,50)
    mask = ~np.isnan(DiameterLevel2_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel2_bin.T, mask.T)]

    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(3, 1, 3)
    
    #plt.title('Grayscale 50%')
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.gca().set_ylim([0, 200])
    plt.ylabel('Diameter level 2, μm')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_xlim([-0.05, 50.05])

    #Diameter 1
    DiameterLevel1_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel1,Particle_Y,0,50,50)
    mask = ~np.isnan(DiameterLevel1_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel1_bin.T, mask.T)]
    plt.subplot(3, 1, 2, sharex=ax1)
    #plt.title('Grayscale 25%')
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Diameter level >=1, μm')
    plt.gca().set_ylim([0, 200])

    #Diameter 0
    DiameterLevel0_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel0,Particle_Y,0,50,50)
    mask = ~np.isnan(DiameterLevel0_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(DiameterLevel0_bin.T, mask.T)]
    plt.subplot(3, 1, 1, sharex=ax1)
    plt.title(PltTitle)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Diameter level >=0, μm')
    plt.gca().set_ylim([0, 200])

    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    #Figurename='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/60um_25_50_75/GrayscaleDiameter.png'
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)    
    
#_______________________________________________________________________________________   
def Plot_BG_diameter(D_KorolevCorr,Diameter, D_BG_KorolevCorr, DiameterBg, Particle_Y,Figurename,PltTitle,ActualDiameter,Y_centre):

    
    #Plot_BG_diameter(D_KorolevCorr,DiameterLevel1, D_BG_KorolevCorr, DiameterBG, Particle_Y,Figurename,PltTitle,60)
    
    # Calculate DOF**********     
    #Y_centre=22
    OASIS_CIP_DOF=CIP_DOF(ActualDiameter,4,70)
    CIP_DOF_st=Y_centre-(OASIS_CIP_DOF/2)
    CIP_DOF_end=Y_centre+(OASIS_CIP_DOF/2)
    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=Y_centre-(Dof_lawson_c6/2)
    Dof_lawson_c6_end=Y_centre+(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=Y_centre-(Dof_lawson_c8/2)
    Dof_lawson_c8_end=Y_centre+(Dof_lawson_c8/2)  
    # ******************************      
    
    # diameter of bounding box korolev corrected
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,0,50,50)
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(4, 1, 1)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-',label='OASIS CIP DOF')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='Lawson et al 2006, c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='Lawson et al 2006, c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.gca().set_ylim([0, 200])
    plt.title(PltTitle)
    plt.legend(loc='upper right')
    plt.ylabel('D (Korolev), μm')
    plt.gca().set_xlim([-0.05, 50.05])
        
    # Diameter of bounding box 
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Particle_Y,0,50,50)
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(4, 1, 2, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('D, μm')
    plt.gca().set_ylim([0, 200])
    
    #Diameter of largest particle Korolev corrected
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_BG_KorolevCorr,Particle_Y,0,50,50)
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(4, 1, 3, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('D BG (Korolev), μm')
    plt.gca().set_ylim([0, 200])
    
    #Diameter of largest particle Korolev corrected
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterBg,Particle_Y,0,50,50)
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(4, 1, 4, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('D BG, μm')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 200])
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)    

#_______________________________________________________________________________________   

#_______________________________________________________________________________________   
def PlotKorCorrD(D_KorolevCorr,Diameter, Particle_Y,Figurename,PltTitle,ActualDiameter):

    # Calculate DOF**********     
    Y_centre=22
    OASIS_CIP_DOF=CIP_DOF(ActualDiameter,4,70)
    CIP_DOF_st=Y_centre-(OASIS_CIP_DOF/2)
    CIP_DOF_end=Y_centre+(OASIS_CIP_DOF/2)
    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=Y_centre-(Dof_lawson_c6/2)
    Dof_lawson_c6_end=Y_centre+(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=Y_centre-(Dof_lawson_c8/2)
    Dof_lawson_c8_end=Y_centre+(Dof_lawson_c8/2)  
    # ******************************      
    
    
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,0,50,50)
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2, 1, 1)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-',label='OASIS CIP DOF')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='Lawson et al 2006, c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='Lawson et al 2006, c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    
    plt.gca().set_ylim([0, 150])
    plt.title(PltTitle)
    plt.legend()
    plt.ylabel('Diameter (Korolev correction), μm')
    
    plt.gca().set_xlim([-0.05, 50.05])
        
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Particle_Y,0,50,50)
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(2, 1, 2, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    
    plt.ylabel('Diameter level >=1, μm')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 150])
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)    




#_______________________________________________________________________________________    
    
# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg


def PlotAFvsD(DropGun_Stats):
    ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['60um_40_50_60_run2','90um_40_50_60','120um_40_50_60'] 
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/'
    PltTitle='120_90_60um_GS_25_50_75'
    CentreOfDOF=22
    fig=plt.figure(figsize=(10,10))
    i=0
    for i in range(len(ExpList)):
        D_Level1_avg=DropGun_Stats[ExpList[i]][4][:]
        BinCentre=DropGun_Stats[ExpList[i]][6][:]
        AF_Level0_avg=DropGun_Stats[ExpList[i]][2][:]
        DOP=np.absolute(BinCentre-CentreOfDOF)  # Distance from object plane 
        
        #plt.plot(BinCentre, D_Level1_avg)
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(D_Level1_avg, AF_Level0_avg, c=DOP, vmin=0, vmax=25,cmap=cm) 

    # DoF Lines
    
    cbar=plt.colorbar(sc)
    cbar.set_label('Distance from object plane, mm')
    plt.ylabel('Level 0/ Total area')
    plt.xlim([0,200])
    plt.xlabel('Diameter, μm')
    plt.title(PltTitle)
    plt.show()
    #plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    #plt.close(fig)
   


#_______________________________________________________________________________________    



#Plot stage positions versus grayscale fraction box plots 

def PlotAreaFractionZd(AreaFraction0,AreaFraction1,AreaFraction2,Zd,Figurename,PltTitle):

    #Area fraction 2
    AreaFraction2_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction2,Zd,-13,13,50)

    mask = ~np.isnan(AreaFraction2_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction2_bin.T, mask.T)]

    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(3, 1, 3)
    
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid, widths=(0.5*(BinsLower[1]-BinsLower[0])))
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.gca().set_ylim([-0.05, 1.05])
    
    plt.ylabel('Level 2 / Total area')
    plt.xlabel('Z$_{d}$')
    plt.gca().set_xlim([-6, 6])

    #Area fraction 1
    AreaFraction1_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction1,Zd,-13,13,50)
    mask = ~np.isnan(AreaFraction1_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction1_bin.T, mask.T)]
    plt.subplot(3, 1, 2, sharex=ax1)
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid, widths=(0.5*(BinsLower[1]-BinsLower[0])))
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Level 1 / Total area')
    plt.gca().set_ylim([-0.05, 1.05])


    #Area fraction 0
    AreaFraction0_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(AreaFraction0,Zd,-13,13,50)
    mask = ~np.isnan(AreaFraction0_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(AreaFraction0_bin.T, mask.T)]
    plt.subplot(3, 1, 1, sharex=ax1)
    plt.title(PltTitle)
    
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid, widths=(0.5*(BinsLower[1]-BinsLower[0])))
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.ylabel('Level 0 / Total area ')
    plt.gca().set_ylim([-0.05, 1.05])
    
    majorLocator = MultipleLocator(1)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(0.5)
    ax1.xaxis.set_major_locator(majorLocator)
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax1.xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
        
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)
    
#_______________________________________________________________________________________    
# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg 
    
    
    
def AllAreaFractionZd(DropGun_Stats,SaveFig):     
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['120um_40_50_60','90um_40_50_60','60um_40_50_60_run2']  
    #ActualDiameterList=[90,75,60] 
    #FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/'
    #CentreOfDOF=22
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    ExpListTmp=np.linspace(50,150,21)
    ExpList=[]
    ActualDiameterList=[]
    CentreOfDOF=20
    for i in range(len(ExpListTmp)-1):
        ExpList.append(str(int(ExpListTmp[i])))
        ActualDiameterList.append((ExpListTmp[i]))
        
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(3, 1, 1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction=DropGun_Stats[ExpList[i]][2][:]
        plt.plot(Zd,AreaFraction,marker='o',markersize=2,linestyle='None')
    plt.ylabel('Level 0 / Total area')
    plt.ylim((-0.05,1.05))
    plt.grid(True)
    
    fig.add_subplot(3, 1, 2,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction=DropGun_Stats[ExpList[i]][1][:]
        plt.plot(Zd,AreaFraction,marker='o',markersize=2,linestyle='None')
    
    plt.ylabel('Level 1 / Total area')
    plt.ylim((-0.05,1.05))
    plt.grid(True)
    
    fig.add_subplot(3, 1, 3,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction=DropGun_Stats[ExpList[i]][0][:]
        Exp=ExpList[i]+', D0 = '+str(ActualDiameterList[i])+' μm'
        #Size='D0 = '+str(ActualDiameterList[i])+' μm'
        plt.plot(Zd,AreaFraction, label=Exp,marker='o',markersize=2,linestyle='None')
    #plt.legend()
    plt.ylabel('Level 2 / Total area')
    plt.ylim((-0.05,1.05))
    #plt.xlim((1,5))
    plt.xlabel('Z$_{d}$')
    plt.grid(True)
    
    Figurename=FilePath+'AF_Zd'
    
    if SaveFig == 1 :
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)

#_______________________________________________________________________________________    
#
#C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/

# Load modelled CIP data

def LoadModelGrayScale(GSpath):

    # Load grayscale stats
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    for filena in os.listdir(GSpath):
        if filena.endswith(".h5"):
            print(filena)
            Data_h5 = h5py.File(GSpath+filena, 'r')              
            TmpTimesWv=np.array(Data_h5['ParticleTimesWv'])
            TmpStats=np.array(Data_h5['ParticleStatsWv'])
            Data_h5.close()  	        
            if 'ParticleTimesWv' not in dir():
                ParticleTimesWv=TmpTimesWv
                ParticleStatsWv=TmpStats
            else:    
                ParticleTimesWv= np.append(ParticleTimesWv,TmpTimesWv,axis=0)
                ParticleStatsWv= np.append(ParticleStatsWv,TmpStats,axis=0)
            
    ActualDiameterAll=ParticleTimesWv[:,3]
    StageXAll=ParticleTimesWv[:,5]     

    #Load stage position
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/60um_25_50_75/'
    #StageFile='StagePositions_exp001.csv'  
    
    return ParticleStatsWv,ActualDiameterAll,StageXAll

#_______________________________________________________________________________________    
#

#VoidRatio=VoidAreaRatio(FilledArea,VoidArea,Particle_Y) 

    
def ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll):

    SavePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    Model_Stats={}
    Model_D_Percentiles={}
    
    
    for D0 in range(50,155,5): 
        #D0=50
        ParticleStatsTmp=ParticleStatsWv[np.where(ActualDiameterAll==D0),:] 
        ParticleStatsExp=ParticleStatsTmp[0,:,:]
        Particle_X=StageXAll[np.where(ActualDiameterAll==D0)]
        AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_X, ParticleStatsExp,0,255075)
        #PlotNumberOfParticles(FilledArea,VoidArea,Particle_X,DiameterLevel0,ParticleStatsExp,SavePath,'Nparticles_VoidRatio_'+str(D0)+'um')
        D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
        #Particle_X_Flag=Particle_X[np.where(DiameterLevel1==DiameterLevel1)]
        #CountsVsStageX(Particle_X_Flag,SavePath,'Counts'+str(int(D0)), D0,20)
        #Figurename='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/Diameter50'
        #PlotGrayscaleDiameter(DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X,Figurename,'D0=50')
        #D_BG_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterBG)
        #Figurename=SavePath+'Korolev_Reuter_correction_D0_'+str(D0)
        #Plot_BG_diameter(D_KorolevCorr,DiameterLevel1, D_BG_KorolevCorr, DiameterBG, Particle_X,Figurename,'D0='+str(D0)+' μm',D0,20)
        
        #VoidRatio = VoidArea / FilledArea
        #D_2stage=Level2_D_Correction_vector(D_KorolevCorr, AreaFraction2, VoidRatio)
        #Plot2StageCorrD(D_2stage, D_KorolevCorr,DiameterLevel1, Particle_X,Figurename,'D0='+str(D0)+' μm',D0,20)
                
        #D_ReuterBakan=ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
        
        #PlotKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,Figurename,'D0='+str(D0)+' μm',D0,20)
        
        #Model_D_Percentiles[str(D0)]=PercentilesKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,D0,20)
        
        #Figurename=SavePath+'Diameter_Z0_'+str(D0)
        #PlotZ0hist(DiameterLevel1, Particle_X,50,Figurename)
        
        Model_Stats[str(D0)]=Avg_Areafraction_Diameter(VoidArea,FilledArea,D_KorolevCorr,AreaFraction0,AreaFraction1,AreaFraction2,DiameterLevel0,DiameterLevel1,DiameterLevel2,Particle_X)

    return Model_Stats

    #return Model_D_Percentiles

#_______________________________________________________________________________________    
#

#def Add_Zd_2_Model_Stats(Model_Stats):

    
#    for D0 in range(50,155,5):      
#        Model_Stats[str(D0)]['Zd'] = DimensionlessDistance(Model_Stats[str(D0)][6][:],SD0,CentreOfDOF_Model)
  



#_______________________________________________________________________________________    
#

def CountsVsStageX(Particle_X,FilePath,PltTitle, ActualDiameter,Y_centre):
    
    fig=plt.figure(figsize=(10,10))
    plt.hist(Particle_X,bins=50,range=(0,50))
    OASIS_CIP_DOF=CIP_DOF(ActualDiameter,4,70)
    CIP_DOF_st=Y_centre-(OASIS_CIP_DOF/2)
    CIP_DOF_end=Y_centre+(OASIS_CIP_DOF/2)
    plt.axvline(x=CIP_DOF_st, color='r', linestyle='-',label='OASIS CIP DOF')
    plt.axvline(x=CIP_DOF_end, color='r', linestyle='-')
    #plt.gca().set_ylim([0, 200])
    plt.title(str(ActualDiameter)+'μm')
    plt.xlabel('Stage Y, mm')
    plt.ylabel('Counts')
    plt.legend(loc='upper right')
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    plt.close(fig)



#_______________________________________________________________________________________    
#    
def ModelGS_Stats_Plot(Model_Stats):
    
    ExpListTmp=np.linspace(50,150,21)
    ExpList=[]
    for i in range(len(ExpListTmp)-1):
        ExpList.append(str(int(ExpListTmp[i])))
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['60um_40_50_60_run2','90um_40_50_60','120um_40_50_60'] 
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    PltTitle='CIP Model'
    fig=plt.figure(figsize=(10,10))
    i=0
    for i in range(len(ExpList)):
        D_Level1_avg=Model_Stats[ExpList[i]][4][:]
        BinCentre=Model_Stats[ExpList[i]][6][:]
        AF_Level0_avg=Model_Stats[ExpList[i]][2][:] 
        #plt.plot(BinCentre, D_Level1_avg)
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(BinCentre, D_Level1_avg, c=AF_Level0_avg,vmin=0, vmax=1,cmap=cm) 
    Dof_D,Dof_Y=DOF_curve(4,20)
    plt.plot(Dof_Y,Dof_D,'o',markersize=1,color='black')    
    cbar=plt.colorbar(sc)
    cbar.set_label('Level 0 / Total area')
    plt.xlabel('Stage Y, mm')
    plt.xlim([-5,55])
    plt.ylabel('Diameter >=1, μm')
    plt.title(PltTitle)
    #plt.show()
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    plt.close(fig)
   
        

#_______________________________________________________________________________________    
#   
    
# Plot korolev corrected size vs position for the model    
    
def ModelGS_Stats_korolev_Plot(Model_Stats):
    
    ExpListTmp=np.linspace(50,150,21)
    ExpList=[]
    for i in range(len(ExpListTmp)-1):
        ExpList.append(str(int(ExpListTmp[i])))
    PltTitle='CIP Model Korolev'
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    fig=plt.figure(figsize=(10,10))
    i=0
    for i in range(len(ExpList)):
        D_Level1_avg=Model_Stats[ExpList[i]][7][:]
        BinCentre=Model_Stats[ExpList[i]][6][:]
        AF_Level0_avg=Model_Stats[ExpList[i]][2][:] 
        #plt.plot(BinCentre, D_Level1_avg)
        cm = plt.cm.get_cmap('RdYlBu')
        sc = plt.scatter(BinCentre, D_Level1_avg, c=AF_Level0_avg,vmin=0, vmax=1,cmap=cm) 

    # DoF Lines
    Dof_D,Dof_Y=DOF_curve(4,20)
    plt.plot(Dof_Y,Dof_D,'o',markersize=1,color='black', label='OASIS CIP DOF')

    #Dof_lawson_c6_D, Dof_lawson_c6_Y=DOF_lawson_curve(6)
    #plt.plot(Dof_lawson_c6_Y,Dof_lawson_c6_D,'o',markersize=1,color='black', label='Lawson et al, c=6')
    
    #Dof_lawson_c8_D, Dof_lawson_c8_Y=DOF_lawson_curve(8)
    #plt.plot(Dof_lawson_c8_Y,Dof_lawson_c8_D,'o',markersize=1,color='black', label='Lawson et al, c=8')
    
    cbar=plt.colorbar(sc)
    cbar.set_label('Level 0 / Total area')
    plt.xlabel('Stage Y, mm')
    plt.xlim([-5,55])
    plt.ylabel('Diameter (Korolev corrected), μm')
    plt.title(PltTitle)
    plt.savefig(FilePath+PltTitle+'.png',dpi=200)
    plt.close(fig)
    
#_______________________________________________________________________________________    

#PLot model and dropgun vs position btween arms


# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg    
    
def CompareModelDropgun(DropGun_Stats,Model_Stats):
    
     
    ExpList='120um_25_50_75_run2'
    DropgunLabel = '120 μm dropgun head'
    ModelSizes=['100','105','110']
    
    #ExpList='90um_25_50_75_run2'   
    #DropgunLabel = '90 μm dropgun head'
    #ModelSizes=['85','90','95']
    
    #ExpList='60um_25_50_75_run2'   
    #DropgunLabel = '60 μm dropgun head'
    #ModelSizes=['70','75','80']
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    CentreOfDOF_Model=20
    CentreOfDOF_DG=22
        
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(3, 1, 3) 
    i=0
    for i in range(len(ModelSizes)):
        D_Level2_avg=Model_Stats[ModelSizes[i]][3][:]
        BinCentre=Model_Stats[ModelSizes[i]][6][:]
        AF_Level2_avg=Model_Stats[ModelSizes[i]][0][:] 
        DOP=BinCentre-CentreOfDOF_Model  # Distance from object plane
        plt.plot(DOP, D_Level2_avg,color='black', label=ModelSizes[i]+' μm')
       
    D_Level2_avg=DropGun_Stats[ExpList][3][:]
    BinCentre=DropGun_Stats[ExpList][6][:]
    #AF_Level2_avg=DropGun_Stats[ExpList][0][:] 
    DOP=BinCentre-CentreOfDOF_DG  # Distance from object plane
    plt.plot(DOP, D_Level2_avg,color='Red', label=DropgunLabel)
    
    plt.legend(loc='upper right')
    plt.ylabel('Diameter level 0, μm')
    plt.ylim((0,200))
    plt.xlabel('Distance from object plane, mm')
    
    
    plt.subplot(3, 1, 2, sharex=ax1)
    for i in range(len(ModelSizes)):
        D_Level1_avg=Model_Stats[ModelSizes[i]][4][:]
        BinCentre=Model_Stats[ModelSizes[i]][6][:]
        AF_Level1_avg=Model_Stats[ModelSizes[i]][1][:] 
        DOP=BinCentre-CentreOfDOF_Model  # Distance from object plane
        plt.plot(DOP, D_Level1_avg,color='black')
    
    D_Level1_avg=DropGun_Stats[ExpList][4][:]
    BinCentre=DropGun_Stats[ExpList][6][:]
    #AF_Level1_avg=DropGun_Stats[ExpList][1][:] 
    DOP=BinCentre-CentreOfDOF_DG  # Distance from object plane
    plt.plot(DOP, D_Level1_avg,color='red')
    
    plt.ylabel('Diameter level >=1, μm')
    plt.ylim((0,200))    
    
    plt.subplot(3, 1, 1, sharex=ax1)
    for i in range(len(ModelSizes)):
        D_Level0_avg=Model_Stats[ModelSizes[i]][5][:]
        BinCentre=Model_Stats[ModelSizes[i]][6][:]
        AF_Level0_avg=Model_Stats[ModelSizes[i]][2][:] 
        DOP=BinCentre-CentreOfDOF_Model  # Distance from object plane
        plt.plot(DOP, D_Level0_avg,color='black')
    
    D_Level0_avg=DropGun_Stats[ExpList][5][:]
    BinCentre=DropGun_Stats[ExpList][6][:]
    #AF_Level0_avg=DropGun_Stats[ExpList][2][:] 
    DOP=BinCentre-CentreOfDOF_DG  # Distance from object plane
    plt.plot(DOP, D_Level0_avg,color='red')
    
    plt.ylabel('Diameter level 2, μm')
    plt.ylim((0,200))    
    plt.xlim((-20, 20))
    plt.title(DropgunLabel)
    
    #plt.show()
    plt.savefig(FilePath+'Model_'+ExpList+'.png',dpi=200)
    


#_______________________________________________________________________________________  

def PLot_Areafraction_6Panel_255075():
    
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
    Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
    
    DropGun_Stats=DropGun_Stats_All()
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 14})
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']
    ExpList=['120um_25_50_75','90um_25_50_75','60um_25_50_75']     
    #ActualDiameterList=[90,80,55] 
    ActualDiameterList=[90,80,55] 
    CentreOfDOF_DG=29
    
    
    DropGunModel_AF_Zd(DropGun_Stats,Model_Stats,1, ExpList, ActualDiameterList, CentreOfDOF_DG, '25', '50', '75')
    DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats,1, ExpList, ActualDiameterList, CentreOfDOF_DG, '25', '50', '75')
    Figurename=FilePath+'AreaFraction_ModelDG'
    plt.savefig(Figurename,dpi=200)
    plt.close(fig) 

#_______________________________________________________________________________________  

def PLot_Areafraction_6Panel_305070():
    
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_305070/'
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
    Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
    
    DropGun_Stats=DropGun_Stats_All()
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']
    ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']     
    ActualDiameterList=[100,85,65] 
    CentreOfDOF_DG=22
    
    
    DropGunModel_AF_Zd(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '30', '50', '70')
    DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '30', '50', '70')
    Figurename=FilePath+'AreaFraction_ModelDG'
    plt.savefig(Figurename,dpi=200)
    plt.close(fig) 

#_______________________________________________________________________________________  

def PLot_Areafraction_6Panel_405070():
    
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_405070/'
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
    Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
    
    DropGun_Stats=DropGun_Stats_All()
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']
    ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']     
    ActualDiameterList=[100,85,65] 
    CentreOfDOF_DG=22
    
    
    DropGunModel_AF_Zd(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '40', '50', '70')
    DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '40', '50', '70')
    Figurename=FilePath+'AreaFraction_ModelDG'
    plt.savefig(Figurename,dpi=200)
    plt.close(fig) 

#_______________________________________________________________________________________  

def PLot_Areafraction_6Panel_255067():
    
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255067/'
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
    Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
    
    DropGun_Stats=DropGun_Stats_All()
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']
    ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']     
    ActualDiameterList=[100,85,65] 
    CentreOfDOF_DG=22
    
    
    DropGunModel_AF_Zd(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '25', '50', '67')
    DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats,0, ExpList, ActualDiameterList, CentreOfDOF_DG, '25', '50', '67')
    Figurename=FilePath+'AreaFraction_ModelDG'
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)



#_______________________________________________________________________________________    

# Plot Area fraction vs Zd for both model and experiments


# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg 
    
    
    
def DropGunModel_AF_Zd(DropGun_Stats,Model_Stats,IncludeDropGun,ExpList, ActualDiameterList, CentreOfDOF_DG, Low, Medium, High):     
        
    
    ModelDiameterList=np.linspace(50,150,21)
    ModelList=[]
    CentreOfDOF_Model=50
    for i in range(len(ModelDiameterList)-1):
        ModelList.append(str(int(ModelDiameterList[i])))    
        
    # Level 0 **************

    plt.subplot(3, 2, 1)
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction=Model_Stats[ModelList[i]][2][:]
        plt.plot(Zd,AreaFraction,marker='o',markersize=2,linestyle='None',color='gray')
    
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction=DropGun_Stats[ExpList[i]][2][:]
            plt.plot(Zd,AreaFraction,marker='+',markersize=5,linestyle='None',color='red')
    
    
    plt.ylabel('A$_{'+Low+' - '+ Medium+'}$ / A$_{'+Low+' - 100}$')
    plt.ylim((-0.05,1.05))
    plt.xlim((-10,10))
    plt.grid(True)
    
    # Level 1 **************
    plt.subplot(3, 2, 3)
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction=Model_Stats[ModelList[i]][1][:]
        plt.plot(Zd,AreaFraction,marker='o',markersize=2,linestyle='None',color='gray')
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction=DropGun_Stats[ExpList[i]][1][:]
            plt.plot(Zd,AreaFraction,marker='+',markersize=5,linestyle='None',color='red')
    
    
    plt.ylabel('A$_{'+Medium+' - '+ High+'}$ / A$_{'+Low+' - 100}$')
    plt.ylim((-0.05,1.05))
    plt.xlim((-10,10))
    plt.grid(True)
    
    # Level 2 **************
    plt.subplot(3, 2, 5)
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction=Model_Stats[ModelList[i]][0][:]
        Exp=ModelList[i]+', D0 = '+str(ModelDiameterList[i])+' μm'
        #Size='D0 = '+str(ActualDiameterList[i])+' μm'
        plt.plot(Zd,AreaFraction, label=Exp,marker='o',markersize=2,linestyle='None',color='gray')
    
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction=DropGun_Stats[ExpList[i]][0][:]
            plt.plot(Zd,AreaFraction,marker='+',markersize=5,linestyle='None',color='red')

    #plt.legend()
    plt.ylabel('A$_{'+High+' - 100}$ / A$_{'+Low+' - 100}$')
    plt.ylim((-0.05,1.05))
    plt.xlim((-10,10))
    plt.xlabel('Z$_{d}$')
    plt.grid(True)


# 
    
#_______________________________________________________________________________________    

# Plot Area fraction vs Zd for both model and experiments


# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg 
    
    
    
def DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats,IncludeDropGun,ExpList, ActualDiameterList, CentreOfDOF_DG, Low, Medium, High):     
        
    ModelDiameterList=np.linspace(50,150,21)
    ModelList=[]
    CentreOfDOF_Model=50
    for i in range(len(ModelDiameterList)-1):
        ModelList.append(str(int(ModelDiameterList[i])))    
        
    # level 2 / level 1 **************
    #fig=plt.figure(figsize=(10,10))
    plt.subplot(3, 2, 2) # level 2 / level 1
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction2=np.array(Model_Stats[ModelList[i]][0][:])
        AreaFraction1=np.array(Model_Stats[ModelList[i]][1][:])
        AreaFraction21=AreaFraction2/AreaFraction1
        plt.plot(Zd,AreaFraction21,marker='o',markersize=2,linestyle='None',color='gray')
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction2=np.array(DropGun_Stats[ExpList[i]][0][:])
            AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
            AreaFraction21=AreaFraction2/AreaFraction1
            plt.plot(Zd,AreaFraction21,marker='+',markersize=5,linestyle='None',color='red')
    
    #plt.ylabel('Area '+High+'% / Area '+Medium+'%')
    plt.ylabel('A$_{'+High+' - 100}$ / A$_{'+Medium+' - '+ High+'}$')
    plt.xlim((-10,10))
    plt.grid(True)
    
    # Level 1 / level 0 **************
    plt.subplot(3, 2, 4)
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction1=np.array(Model_Stats[ModelList[i]][1][:])
        AreaFraction0=np.array(Model_Stats[ModelList[i]][2][:])
        AreaFraction10=AreaFraction1/AreaFraction0
        plt.plot(Zd,AreaFraction10,marker='o',markersize=2,linestyle='None',color='gray')
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
            AreaFraction0=np.array(DropGun_Stats[ExpList[i]][2][:])
            AreaFraction10=AreaFraction1/AreaFraction0
            plt.plot(Zd,AreaFraction10,marker='+',markersize=5,linestyle='None',color='red')
    
    
    #plt.ylabel('Area '+Medium+'% / Area '+Low+'%')
    plt.ylabel('A$_{'+Medium+' - '+ High+'}$ / A$_{'+Low+' - '+ Medium+'}$')
    plt.xlim((-10,10))
    plt.grid(True)
    
    # Level 2 / level 0 **************
    plt.subplot(3, 2, 6)
    i=0
    for i in range(len(ModelList)):
        ActualDiameter=ModelDiameterList[i]
        Zd=DimensionlessDistance(Model_Stats[ModelList[i]][6][:],ActualDiameter,CentreOfDOF_Model)
        AreaFraction2=np.array(Model_Stats[ModelList[i]][0][:])
        AreaFraction0=np.array(Model_Stats[ModelList[i]][2][:])
        AreaFraction20=AreaFraction2/AreaFraction0      
        Exp=ModelList[i]+', D0 = '+str(ModelDiameterList[i])+' μm'
        #Size='D0 = '+str(ActualDiameterList[i])+' μm'
        plt.plot(Zd,AreaFraction20, label=Exp,marker='o',markersize=2,linestyle='None',color='gray')
    
    if IncludeDropGun == 1 :
        for i in range(len(ExpList)):
            ActualDiameter=ActualDiameterList[i]
            Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF_DG)
            AreaFraction2=np.array(DropGun_Stats[ExpList[i]][0][:])
            AreaFraction0=np.array(DropGun_Stats[ExpList[i]][2][:])
            AreaFraction20= AreaFraction2 / AreaFraction0
            plt.plot(Zd,AreaFraction20,marker='+',markersize=5,linestyle='None',color='red')

    #plt.legend()
    #plt.ylabel('Area '+High+'% / Area '+Low+'%')
    plt.ylabel('A$_{'+High+' - 100}$ / A$_{'+Low+' - '+ Medium+'}$')
    plt.xlim((-10,10))
    plt.xlabel('Z$_{d}$')
    plt.grid(True)
    #Figurename=FilePath+'AF_Zd_DG_2'
    #plt.savefig(Figurename,dpi=200)
    #plt.close(fig) 






















    
#_______________________________________________________________________________________    

# Plot ratios Level2/Level1, level 1/ level0, and level2/level0     



# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg 
    

  
  
def AllAreaFractionZd_2(DropGun_Stats,SaveFig):     
    #ExpList=['120um_25_50_75_run2','90um_25_50_75_run2','60um_25_50_75_run2']    
    #ExpList=['120um_40_50_60','90um_40_50_60','60um_40_50_60_run2']  
    #ActualDiameterList=[90,75,60] 
    #FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/'
    #CentreOfDOF=22
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    ExpListTmp=np.linspace(50,150,21)
    ExpList=[]
    ActualDiameterList=[]
    CentreOfDOF=20
    for i in range(len(ExpListTmp)-1):
        ExpList.append(str(int(ExpListTmp[i])))
        ActualDiameterList.append((ExpListTmp[i]))
        
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(4, 1, 1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction2=np.array(DropGun_Stats[ExpList[i]][0][:])
        AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
        Ratio=AreaFraction2/AreaFraction1
        plt.plot(Zd,Ratio,marker='o',markersize=2,linestyle='None')
    plt.ylabel('Level 2 / Level 1')
    #plt.ylim((-0.05,1.05))
    plt.grid(True)
    
    fig.add_subplot(4, 1, 2,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
        AreaFraction0=np.array(DropGun_Stats[ExpList[i]][2][:])
        Ratio=AreaFraction1/AreaFraction0
        plt.plot(Zd,Ratio,marker='o',markersize=2,linestyle='None')
    
    plt.ylabel('Level 1 / Level 0')
    #plt.ylim((-0.05,1.05))
    plt.grid(True)
    
    fig.add_subplot(4, 1, 3,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction2=np.array(DropGun_Stats[ExpList[i]][0][:])
        AreaFraction0=np.array(DropGun_Stats[ExpList[i]][2][:])
        Ratio=AreaFraction2/AreaFraction0
        Exp=ExpList[i]+', D0 = '+str(ActualDiameterList[i])+' μm'
        plt.plot(Zd,Ratio, label=Exp,marker='o',markersize=2,linestyle='None')
    #plt.legend()
    plt.ylabel('Level 2 / Level 0')
    #plt.ylim((-0.05,1.05))
    plt.grid(True)
    
    
    fig.add_subplot(4, 1, 4,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        AreaFraction2=np.array(DropGun_Stats[ExpList[i]][0][:])
        AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
        AreaFraction0=np.array(DropGun_Stats[ExpList[i]][2][:])
        Ratio=(AreaFraction2+AreaFraction1)/AreaFraction0
        Exp=ExpList[i]+', D0 = '+str(ActualDiameterList[i])+' μm'
        plt.plot(Zd,Ratio, label=Exp,marker='o',markersize=2,linestyle='None')
    #plt.legend()
    plt.ylabel('(Level 2 + Level 1) / Level 0')
    #plt.ylim((-0.05,1.05))

    
    plt.xlabel('Z$_{d}$')
    #plt.xlim((1,5))
    plt.grid(True)
    
    Figurename=FilePath+'AF_Zd_MoreRatios'
    
    if SaveFig==1:    
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)
   
    
    

#_______________________________________________________________________________________    

# Plot korolev correction vs Zd
    
# DropGun_Stats columns: AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre,D_KorolevCorr_avg 



def KorolevCorr_vs_Zd(DropGun_Stats,SaveFlag):
    
    FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
    ExpListTmp=np.linspace(50,150,21)
    ExpList=[]
    ActualDiameterList=[]
    CentreOfDOF=20
    for i in range(len(ExpListTmp)-1):
        ExpList.append(str(int(ExpListTmp[i])))
        ActualDiameterList.append((ExpListTmp[i]))
        
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2, 1, 1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        D_Level1=np.array(DropGun_Stats[ExpList[i]][4][:])
        D_KorolevCorr=np.array(DropGun_Stats[ExpList[i]][7][:])
        
        #D_KorolevCorr=Diameter*CorrFactArray
        CorrFact=D_KorolevCorr/D_Level1
        
        #AreaFraction1=np.array(DropGun_Stats[ExpList[i]][1][:])
        #Ratio=AreaFraction2/AreaFraction1
        plt.plot(Zd,CorrFact,marker='o',markersize=2,linestyle='None')
    plt.ylabel('Korolev correction factor')
    
    fig.add_subplot(2, 1, 2,sharex=ax1)
    i=0
    for i in range(len(ExpList)):
        ActualDiameter=ActualDiameterList[i]
        Zd=DimensionlessDistance(DropGun_Stats[ExpList[i]][6][:],ActualDiameter,CentreOfDOF)
        D_Level1=np.array(DropGun_Stats[ExpList[i]][4][:])
        plt.plot(Zd,D_Level1,marker='o',markersize=2,linestyle='None')
    
    plt.ylabel('Diameter level >=1, μm')
    plt.xlabel('Z$_{d}$')
    #plt.ylim((-0.05,1.05))
    plt.grid(True)

    if SaveFlag == 1 : 
        Figurename=FilePath+'KorolevCorr_Zd'  
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)

    
    
#_______________________________________________________________________________________    
    
# Use area fraction to calculate Zd

def AF_to_Zd(Model_Stats,FilePath,SaveFlag):     
    
    Zd = []
    Ratio_2_0 = []  
    Ratio_1_0 = []
    
    CentreOfDOF_Model=50
    for D0 in range(50,155,5):      
        TmpZd = DimensionlessDistance(Model_Stats[str(D0)][6][:],D0,CentreOfDOF_Model)
        AreaFraction2=np.array(Model_Stats[str(D0)][0][:])
        AreaFraction1=np.array(Model_Stats[str(D0)][1][:])
        AreaFraction0=np.array(Model_Stats[str(D0)][2][:])
        TmpRatio_2_0=AreaFraction2/AreaFraction0
        TmpRatio_1_0=AreaFraction1/AreaFraction0
    
        Zd = np.append(Zd,TmpZd)
        Ratio_2_0 = np.append(Ratio_2_0,TmpRatio_2_0)
        Ratio_1_0 = np.append(Ratio_1_0,TmpRatio_1_0)

    #Ratio_2_0=(np.where((Zd<0) | (Zd>5), np.nan, Ratio_2_0))
    #Ratio_1_0=(np.where((Zd<3.5) , np.nan, Ratio_1_0))
    
    Ratio_1_0=(np.where( (Zd<2) | (Ratio_2_0>0) | (Ratio_1_0 == 0) , np.nan, Ratio_1_0))
    
    
    #Zd = (np.where((Zd<0) | (Zd>5), np.nan, Zd))
    #idx = np.isfinite(Zd) & np.isfinite(Ratio_2_0)
    #p = np.polyfit(Zd[idx], Ratio_2_0[idx], 3)
    #Fit_Ratio_2_0 = p[0] * Zd**3 + p[1] * Zd**2 + p[2] * Zd + p[3] # 3rd degree polynomial
    
    #plt.plot(Zd[idx],Ratio_2_0[idx],'o')
    #plt.xscale('log')
    #plt.plot(Zd[idx],Fit_Ratio_2_0[idx],'o',color='black')
    
    idx=np.isfinite(Zd) & np.isfinite(Ratio_1_0)
    p = np.polyfit(Zd[idx], Ratio_1_0[idx], 1)
    Fit_Ratio_1_0 = p[0] * Zd + p[1] # linear fit 
    
    print('pearsonr')
    print(pearsonr(Zd[idx], Ratio_1_0[idx]))
    
    #plt.rcParams.update({'font.size': 12})
    
    fig, ax = plt.subplots()
    plt.plot(Zd[idx],Ratio_1_0[idx],'o')
    #plt.xscale('log')
    plt.plot(Zd[idx],Fit_Ratio_1_0[idx],color='black')
    plt.grid(True)
    plt.ylabel('A$_{50 - 75}$ / $A_{25 - 50}$')
    plt.xlabel('Z$_{d}$')   
    textstr = 'y = '+ str(round(p[0],4))+'x + '+str(round(p[1],4))
    props = dict(boxstyle='square', facecolor='white', alpha=0.5)
    #props = dict(boxstyle='square', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14,
        verticalalignment='top',horizontalalignment='right', bbox=props)

    if SaveFlag == 1 : 
        #FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
        Figurename=FilePath+'Ratio_1_0_Zd'  
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)

    
    
    
    
    #return Ratio_2_0, Zd
    return p  # Coefficients to the fit Zd vs Level1/Level0

#_______________________________________________________________________________________    

# Use level 2 ratio to find Zd when less than 0.5

def Level2_2_Zd(Model_Stats,SaveFlag):
    
    #ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale()
    #Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)

    Zd = []
    Ratio_2_0 = []  
    Ratio_1_0 = []
    AreaFraction2 = []  
    VoidRatio = []
    CentreOfDOF_Model=20

    for D0 in range(50,155,5):      
        TmpVoidRatio=np.array(Model_Stats[str(D0)][8][:])
        TmpZd = DimensionlessDistance(Model_Stats[str(D0)][6][:],D0,CentreOfDOF_Model)
        TmpAreaFraction2=np.array(Model_Stats[str(D0)][0][:])
        TmpAreaFraction1=np.array(Model_Stats[str(D0)][1][:])
        TmpAreaFraction0=np.array(Model_Stats[str(D0)][2][:])
        TmpRatio_2_0=TmpAreaFraction2/TmpAreaFraction0
        TmpRatio_1_0=TmpAreaFraction1/TmpAreaFraction0
    
        Zd = np.append(Zd,TmpZd)
        Ratio_2_0 = np.append(Ratio_2_0,TmpRatio_2_0)
        Ratio_1_0 = np.append(Ratio_1_0,TmpRatio_1_0)
        AreaFraction2 = np.append(AreaFraction2,TmpAreaFraction2)
        VoidRatio = np.append(VoidRatio,TmpVoidRatio)
    
    #TmpAreaFraction2=(np.where( (Zd>0.5) | (Zd<0) , np.nan, AreaFraction2))
    TmpAreaFraction2=(np.where( (AreaFraction2<0.55) | (AreaFraction2>0.8) | (Zd<0), np.nan, AreaFraction2))   

    idx=np.isfinite(Zd) & np.isfinite(TmpAreaFraction2)
    p = np.polyfit(Zd[idx], TmpAreaFraction2[idx], 2)
    
    SortZd=np.sort(Zd[idx]) 
    Fit_AreaFraction2 = p[0] *(SortZd**2) + p[1] * SortZd + p[2] # Fit curve

    #p = np.polyfit(Zd[idx], TmpAreaFraction2[idx], 1)
    #Fit_AreaFraction2 =  p[0] * Zd + p[1] # Fit curve
   
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2, 1, 1)

    plt.plot(Zd,TmpAreaFraction2,'o')
    plt.plot(SortZd,Fit_AreaFraction2)
    plt.ylabel('Level 2 / Total Area')

    textstr = 'y = '+ str(round(p[0],4))+'$x^2$ + '+str(round(p[1],4))+'x + '+str(round(p[2],4))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #props = dict(boxstyle='round', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=ax1.transAxes, fontsize=14,
        verticalalignment='top',horizontalalignment='right', bbox=props)

    plt.subplot(2, 1, 2, sharex=ax1) 
    plt.plot(Zd,VoidRatio,'o')
    plt.ylabel('Void Area / Filled Area')
    plt.xlabel('Z$_{d}$')  

    plt.ylim((0,0.2))
    plt.xlim((0,2))

    if SaveFlag == 1 : 
        FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
        Figurename=FilePath+'Level2_Zd'  
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)

    return p


#_______________________________________________________________________________________   

# Plot diameter with korolev and areafraction2 corrections

def Plot2StageCorrD(D_2stage, D_KorolevCorr,Diameter, Particle_Y,Figurename,PltTitle,ActualDiameter,Y_centre):

    # Calculate DOF**********     
    #Y_centre=22
    #OASIS_CIP_DOF=CIP_DOF(ActualDiameter,4,70)
    #CIP_DOF_st=Y_centre-(OASIS_CIP_DOF/2)
    #CIP_DOF_end=Y_centre+(OASIS_CIP_DOF/2)
    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=Y_centre-(Dof_lawson_c6/2)
    Dof_lawson_c6_end=Y_centre+(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=Y_centre-(Dof_lawson_c8/2)
    Dof_lawson_c8_end=Y_centre+(Dof_lawson_c8/2)  
    Dof_lawson_c05= Dof_lawson(ActualDiameter,Lambda,0.5)
    Dof_lawson_c05_st=Y_centre-(Dof_lawson_c05/2)
    Dof_lawson_c05_end=Y_centre+(Dof_lawson_c05/2)     
    
    # ******************************      
    
    

    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    
    
    D_2stage_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_2stage,Particle_Y,0,50,50)
    mask = ~np.isnan(D_2stage_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_2stage_bin.T, mask.T)]
    ax1=fig.add_subplot(3, 1, 1)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c05_st, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c05_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('Diameter (Both corrections), μm')
    #plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 200])
    plt.title(PltTitle)
    
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,0,50,50)
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    plt.subplot(3, 1, 2, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c05_st, color='r', linestyle='-',label='Lawson et al 2006, c=0.5')
    plt.axvline(x=Dof_lawson_c05_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='Lawson et al 2006, c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='Lawson et al 2006, c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.gca().set_ylim([0, 200])
    plt.legend()
    plt.ylabel('Diameter (Korolev), μm')
    plt.gca().set_xlim([-0.05, 50.05])    
    
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Particle_Y,0,50,50)
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(3, 1, 3, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c05_st, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c05_end, color='r', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('Diameter, μm')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 200])
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)    

#_______________________________________________________________________________________    

# Use level 2 ratio to find Zd when less than 0.5

def Level2_level1_2_Zd(Model_Stats,SaveFlag):
    
    #ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale()
    #Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)

    Zd = []
    Ratio_2_0 = []  
    Ratio_1_0 = []
    AreaFraction2 = []  
    AreaFraction1 = [] 
    VoidRatio = []
    CentreOfDOF_Model=20

    for D0 in range(50,155,5):      
        TmpVoidRatio=np.array(Model_Stats[str(D0)][8][:])
        TmpZd = DimensionlessDistance(Model_Stats[str(D0)][6][:],D0,CentreOfDOF_Model)
        TmpAreaFraction2=np.array(Model_Stats[str(D0)][0][:])
        TmpAreaFraction1=np.array(Model_Stats[str(D0)][1][:])
        TmpAreaFraction0=np.array(Model_Stats[str(D0)][2][:])
        TmpRatio_2_0=TmpAreaFraction2/TmpAreaFraction0
        TmpRatio_1_0=TmpAreaFraction1/TmpAreaFraction0
    
        Zd = np.append(Zd,TmpZd)
        Ratio_2_0 = np.append(Ratio_2_0,TmpRatio_2_0)
        Ratio_1_0 = np.append(Ratio_1_0,TmpRatio_1_0)
        AreaFraction2 = np.append(AreaFraction2,TmpAreaFraction2)
        AreaFraction1 = np.append(AreaFraction1,TmpAreaFraction1)
        VoidRatio = np.append(VoidRatio,TmpVoidRatio)
    
    Ratio_2_1= AreaFraction2 / AreaFraction1
    
    #TmpAreaFraction2=(np.where( (Zd>0.5) | (Zd<0) , np.nan, AreaFraction2))
    TmpRatio_2_1=(np.where( (Ratio_2_1 < 1.5) | (Ratio_2_1>6) | (Zd<0), np.nan, Ratio_2_1))   

    idx=np.isfinite(Zd) & np.isfinite(TmpRatio_2_1)
    p = np.polyfit(Zd[idx], TmpRatio_2_1[idx], 2)
    
    SortZd=np.sort(Zd[idx]) 
    Fit_Ratio_2_1 = p[0] *(SortZd**2) + p[1] * SortZd + p[2] # Fit curve

    #p = np.polyfit(Zd[idx], TmpAreaFraction2[idx], 1)
    #Fit_AreaFraction2 =  p[0] * Zd + p[1] # Fit curve
   
    fig=plt.figure(figsize=(10,10))
    ax1=fig.add_subplot(2, 1, 1)

    plt.plot(Zd,TmpRatio_2_1,'o')
    plt.plot(SortZd,Fit_Ratio_2_1)
    plt.ylabel('Level 2 / Level 1')

    textstr = 'y = '+ str(round(p[0],4))+'$x^2$ + '+str(round(p[1],4))+'x + '+str(round(p[2],4))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #props = dict(boxstyle='round', alpha=0.5)
    plt.text(0.95, 0.95, textstr, transform=ax1.transAxes, fontsize=14,
        verticalalignment='top',horizontalalignment='right', bbox=props)

    plt.subplot(2, 1, 2, sharex=ax1) 
    plt.plot(Zd,VoidRatio,'o')
    plt.ylabel('Void Area / Filled Area')
    plt.xlabel('Zd')  

    plt.ylim((0,0.2))
    plt.xlim((0,2))

    if SaveFlag == 1 : 
        FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP/'
        Figurename=FilePath+'Level2_level1_Zd'  
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)

    return p

#_______________________________________________________________________________________  
    
# Reuter and Bakan Fig 3
 #AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre

def DiameterVSAreaFraction1(Model_Stats):

    #Diameter='90'
    for D0 in range(50,155,25):  
        D_Level1_avg=Model_Stats[str(D0)][4][:]
        AF_Level0_avg=Model_Stats[str(D0)][2][:]
        plt.plot(AF_Level0_avg,D_Level1_avg, 'o')
    
    
    
    plt.ylabel('Diameter level1')
    plt.xlabel('Level 0/ Total area')
    
    


#_______________________________________________________________________________________



def Model_Diameter_Percentiles(ParticleStatsWv,ActualDiameterAll,StageXAll):

    #Model_D_Percentiles={}
    
    Model_D= np.linspace(50,150,21)
    #Output=[]
    
    #for D0 in range(50,155,5):
    for D0 in Model_D : 
        #D0=50
        ParticleStatsTmp=ParticleStatsWv[np.where(ActualDiameterAll==D0),:] 
        ParticleStatsExp=ParticleStatsTmp[0,:,:]
        Particle_X=StageXAll[np.where(ActualDiameterAll==D0)]
        AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_X, ParticleStatsExp,0,255075)
        D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
                
        D_ReuterBakan=ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
        
        #PlotKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,Figurename,'D0='+str(D0)+' μm',D0,20)
        
        #Model_D_Percentiles[str(D0)]=PercentilesKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,D0,20)
        Tmp=np.asarray(PercentilesKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,D0,50))
        #Output= np.append(Output,Tmp)
        
        if D0 == 50:
                Model_D_Percentiles=Tmp[np.newaxis]
                #print(D0)
        else:    
                Model_D_Percentiles= np.append(Model_D_Percentiles,Tmp[np.newaxis],axis=0)
                #print(D0)
    

    



            
    return Model_D_Percentiles, Model_D


#_______________________________________________________________________________________   

#Bin diameter based on poisionts 


def PercentilesKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,Diameter, Particle_Y,ActualDiameter,Y_centre):

    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6) / 2 # Note divide by 2
    
    D_ReuterBakan_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_ReuterBakan,Particle_Y,0,100,100)
    D_ReuterBakan_binMedian=np.nanmedian(D_ReuterBakan_bin,axis=0)
    
    D_ReuterBakan_median=np.nanmedian(D_ReuterBakan_binMedian)
    D_ReuterBakan_25=np.nanpercentile(D_ReuterBakan_binMedian,25)
    D_ReuterBakan_75=np.nanpercentile(D_ReuterBakan_binMedian,75)

    D_ReuterBakan_bin_c6=(np.where((BinsMid < (Y_centre + Dof_lawson_c6)) & (BinsMid > (Y_centre - Dof_lawson_c6)),D_ReuterBakan_binMedian, np.nan))
    D_ReuterBakan_c6_median=np.nanmedian(D_ReuterBakan_bin_c6)
    D_ReuterBakan_c6_25=np.nanpercentile(D_ReuterBakan_bin_c6,25)
    D_ReuterBakan_c6_75=np.nanpercentile(D_ReuterBakan_bin_c6,75)


    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,0,100,100)
    D_KorolevCorr_binMedian=np.nanmedian(D_KorolevCorr_bin,axis=0)
    D_Korolev_median=np.nanmedian(D_KorolevCorr_binMedian)
    D_Korolev_25=np.nanpercentile(D_KorolevCorr_binMedian,25)
    D_Korolev_75=np.nanpercentile(D_KorolevCorr_binMedian,75)
    
    D_KorolevCorr_bin_c6=(np.where((BinsMid < (Y_centre + Dof_lawson_c6)) & (BinsMid > (Y_centre - Dof_lawson_c6)),D_KorolevCorr_binMedian, np.nan))
    D_Korolev_c6_median=np.nanmedian(D_KorolevCorr_bin_c6)
    D_Korolev_c6_25=np.nanpercentile(D_KorolevCorr_bin_c6,25)
    D_Korolev_c6_75=np.nanpercentile(D_KorolevCorr_bin_c6,75)

    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Particle_Y,0,100,100)
    Diameter_binMedian=np.nanmedian(Diameter_bin,axis=0)
    D_median=np.nanmedian(Diameter_binMedian)
    D_25=np.nanpercentile(Diameter_binMedian,25)
    D_75=np.nanpercentile(Diameter_binMedian,75)
    
    
    Diameter_bin_c6=(np.where((BinsMid < (Y_centre + Dof_lawson_c6)) & (BinsMid > (Y_centre - Dof_lawson_c6)),Diameter_binMedian, np.nan))
    D_c6_median=np.nanmedian(Diameter_bin_c6)
    D_c6_25=np.nanpercentile(Diameter_bin_c6,25)
    D_c6_75=np.nanpercentile(Diameter_bin_c6,75)


    #plt.hist(Diameter_bin_c6, bins=20, range=(90,160))

    return D_ReuterBakan_median, D_ReuterBakan_25,D_ReuterBakan_75, \
        D_ReuterBakan_c6_median,D_ReuterBakan_c6_25,D_ReuterBakan_c6_75, \
        D_Korolev_median,D_Korolev_25,D_Korolev_75, \
        D_Korolev_c6_median, D_Korolev_c6_25, D_Korolev_c6_75, \
        D_median, D_25, D_75, \
        D_c6_median, D_c6_25, D_c6_75
    
    

#_______________________________________________________________________________________   

# Plot diameter with korolev and areafraction2 corrections

def PlotKorolevReuterDiameter(D_ReuterBakan, D_KorolevCorr,Diameter, Particle_Y,Figurename,PltTitle,ActualDiameter,Y_centre):

    # Calculate DOF**********     
    #Y_centre=22
    #OASIS_CIP_DOF=CIP_DOF(ActualDiameter,4,70)
    #CIP_DOF_st=Y_centre-(OASIS_CIP_DOF/2)
    #CIP_DOF_end=Y_centre+(OASIS_CIP_DOF/2)
    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=Y_centre-(Dof_lawson_c6/2)
    Dof_lawson_c6_end=Y_centre+(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=Y_centre-(Dof_lawson_c8/2)
    Dof_lawson_c8_end=Y_centre+(Dof_lawson_c8/2)      
    
    # ******************************      
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})   
    
    D_ReuterBakan_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_ReuterBakan,Particle_Y,0,50,50)
    
    #D_ReuterBakan_median=np.nanmedian(D_ReuterBakan_bin)
    #D_ReuterBakan_25=np.nanpercentile(D_ReuterBakan_bin,25)
    #D_ReuterBakan_75=np.nanpercentile(D_ReuterBakan_bin,75)
    
    mask = ~np.isnan(D_ReuterBakan_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_ReuterBakan_bin.T, mask.T)]
    ax1=fig.add_subplot(3, 1, 1)
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('Diameter (Reuter & Bakan), μm')
    #plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 200])
    plt.title(PltTitle)
    
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_Y,0,50,50)
    
    #D_Korolev_median=np.nanmedian(D_KorolevCorr_bin)
    #D_Korolev_25=np.nanpercentile(D_KorolevCorr_bin,25)
    #D_Korolev_75=np.nanpercentile(D_KorolevCorr_bin,75)
       
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    plt.subplot(3, 1, 2, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='Lawson et al 2006, c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='Lawson et al 2006, c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.gca().set_ylim([0, 200])
    plt.legend()
    plt.ylabel('Diameter (Korolev), μm')
    plt.gca().set_xlim([-0.05, 50.05])    
    
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Particle_Y,0,50,50)
    
    #D_median=np.nanmedian(Diameter_bin)
    #D_25=np.nanpercentile(Diameter_bin,25)
    #D_75=np.nanpercentile(Diameter_bin,75)
       
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.subplot(3, 1, 3, sharex=ax1) 
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    plt.plot(BinsMid, means,"o", color='black')
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.ylabel('Diameter, μm')
    plt.xlabel('Stage Y, mm')
    plt.gca().set_ylim([0, 200])
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    plt.savefig(Figurename,dpi=200)
    plt.close(fig)   

#_______________________________________________________________________________________   

def Model_Dropgun_Diameter_6panel():
    
 

    #***************** March 2018 ********************

    #60um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/20180327143000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/60um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='60um_25_50_75_run2'
    #ActualDiameter=65

    #90um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/20180326150258/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/90um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 26, 1, 0, 0) #1 hour time difference between stage and CIP
    #PltTitle='90um_25_50_75_run2'
    #ActualDiameter=85

    #120um_25_50_75_run2
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/20180327100000/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/March2018/120um_25_50_75_run2/'
    #ExpDate= datetime.datetime(2018, 3, 27, 1, 0, 0) #1 hour time difference between stage and CIP   
    #PltTitle='120um_25_50_75_run2'
    #ActualDiameter=100

    #***************** September 2018 ********************

    #60um_25_50_75
#    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/20180927111105/Output/'
#    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/'
#    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
#    PltTitle='60um_25_50_75'
#    ActualDiameter=55

    #90um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/20180927134617/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    PltTitle='90um_25_50_75'
    ActualDiameter=80    

    #120um_25_50_75
    #GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/Both/Output/'
    #StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/'
    #ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    #PltTitle='120um_25_50_75'
    #ActualDiameter=90



    StageFile='StagePositions_exp001.csv'       
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    D_ReuterBakan=ReuterBakanDiameter(DiameterLevel0,AreaFraction0)
    fig=plt.figure(figsize=(10,8))
    plt.rcParams.update({'font.size': 12}) 
    #PlotDiameter6Panel(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,ActualDiameter,29.5, 1)
    PlotDiameter4Panel(D_ReuterBakan, D_KorolevCorr,DiameterLevel1, Particle_X,ActualDiameter,29.5, 1)


    #***************** Model ********************
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    ModelStatsWv,ActualDiameterAll,ModelXAll=LoadModelGrayScale(GSpath)
    ModelStatsTmp=ModelStatsWv[np.where(ActualDiameterAll==ActualDiameter),:] 
    ModelStatsExp=ModelStatsTmp[0,:,:]
    Model_X=ModelXAll[np.where(ActualDiameterAll==ActualDiameter)]
    ModelAreaFraction0,ModelAreaFraction1,ModelAreaFraction2,ModelAreaRatio_2_1,ModelAreaRatio_1_0,ModelDiameterLevel0,ModelDiameterLevel1,ModelDiameterLevel2, ModelFilledArea, ModelVoidArea,ModelDiameterBG,ModelZd_fromRatios=GrayscaleRatios(Model_X, Model_X, ModelStatsExp,0,255075)
    ModelD_KorolevCorr=KorolevCorrectedD(ModelFilledArea, ModelVoidArea,ModelDiameterLevel1)
    ModelD_ReuterBakan=ReuterBakanDiameter(ModelDiameterLevel0,ModelAreaFraction0)
    #PlotDiameter6Panel(ModelD_ReuterBakan, ModelD_KorolevCorr,ModelDiameterLevel1, Model_X,ActualDiameter,50, 2)
    PlotDiameter4Panel(ModelD_ReuterBakan, ModelD_KorolevCorr,ModelDiameterLevel1, Model_X,ActualDiameter,50, 2)


    Figurename='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'+PltTitle+'_Dcorrections4panels'

    plt.savefig(Figurename,dpi=200)
    plt.close(fig)  


#_______________________________________________________________________________________   

# Plot diameter with korolev and areafraction2 corrections

def PlotDiameter6Panel(D_ReuterBakan, D_KorolevCorr,Diameter, Particle_Y,ActualDiameter,Y_centre, Strip):

    # Calculate DOF**********

    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=-1*(Dof_lawson_c6/2)
    Dof_lawson_c6_end=(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=-1*(Dof_lawson_c8/2)
    Dof_lawson_c8_end=(Dof_lawson_c8/2)      
    
    # ******************************      
    
    Z=Particle_Y#-Y_centre

    ax1=plt.subplot(3, 2, Strip)
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Z,0,100,100)       
    BinsMid-=Y_centre
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
    
    plt.ylabel('Diameter, μm')    
    plt.gca().yaxis.grid(True)
    plt.gca().set_ylim([0, 200])
    
    if Strip == 1 : 
        plt.title('Drop generator')
        #plt.axhspan(ymin=0, ymax=200, xmin=-50, xmax=-28.9, facecolor='gray')
        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
    else:
        plt.title('Model')
        plt.legend()
    
    plt.subplot(3, 2, Strip+2, sharex=ax1) 
    D_ReuterBakan_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_ReuterBakan,Z,0,100,100)
    BinsMid-=Y_centre
    mask = ~np.isnan(D_ReuterBakan_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_ReuterBakan_bin.T, mask.T)]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
    plt.gca().set_ylim([0, 200])
    plt.gca().yaxis.grid(True)
    #plt.legend()
    plt.ylabel('Diameter (RB98), μm')
    if Strip == 1 : 
        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
    
  
    
    plt.subplot(3, 2, Strip+4, sharex=ax1) 
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Z,0,100,100)       
    BinsMid-=Y_centre
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
    plt.ylabel('Diameter (K07), μm')
    plt.xlabel('Distance from object plane, mm')
    plt.gca().set_ylim([0, 200])
    plt.gca().set_xlim([-39, 39]) 
    plt.gca().yaxis.grid(True)
    if Strip == 1 : 
        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    #plt.savefig(Figurename,dpi=200)
    #plt.close(fig)   

#_______________________________________________________________________________________   

# Plot diameter with korolev 

def PlotDiameter4Panel(D_ReuterBakan, D_KorolevCorr,Diameter, Particle_Y,ActualDiameter,Y_centre, Strip):

    # Calculate DOF**********

    Lambda = 0.658
    Dof_lawson_c6= Dof_lawson(ActualDiameter,Lambda,6)
    Dof_lawson_c6_st=-1*(Dof_lawson_c6/2)
    Dof_lawson_c6_end=(Dof_lawson_c6/2)   
    Dof_lawson_c8= Dof_lawson(ActualDiameter,Lambda,8)
    Dof_lawson_c8_st=-1*(Dof_lawson_c8/2)
    Dof_lawson_c8_end=(Dof_lawson_c8/2)      
    
    # ******************************      
    
    Z=Particle_Y#-Y_centre

    ax1=plt.subplot(2, 2, Strip)
    Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(Diameter,Z,0,100,100)       
    BinsMid-=Y_centre
    mask = ~np.isnan(Diameter_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='c=6')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-',label='c=8')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
    
    plt.ylabel('Diameter, μm')    
    plt.gca().yaxis.grid(True)
    plt.gca().set_ylim([0, 200])
    
    if Strip == 1 : 
        plt.title('Drop generator')
        #plt.axhspan(ymin=0, ymax=200, xmin=-50, xmax=-28.9, facecolor='gray')
        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
    else:
        plt.title('Model')
        plt.legend()
    
#    plt.subplot(3, 2, Strip+2, sharex=ax1) 
#    D_ReuterBakan_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_ReuterBakan,Z,0,100,100)
#    BinsMid-=Y_centre
#    mask = ~np.isnan(D_ReuterBakan_bin) # remove nans
#    filtered_data = [d[m] for d, m in zip(D_ReuterBakan_bin.T, mask.T)]
#    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
#    means = [np.nanmean(x) for x in filtered_data]
#    #Sdev=[np.std(x) for x in filtered_data]
#    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
#    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
#    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
#    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
#    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
#    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
#    plt.gca().set_ylim([0, 150])
#    plt.gca().yaxis.grid(True)
#    #plt.legend()
#    plt.ylabel('Diameter (RB98), μm')
#    if Strip == 1 : 
#        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
#        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
#    
  
    
    plt.subplot(2, 2, Strip+2, sharex=ax1) 
    D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Z,0,100,100)       
    BinsMid-=Y_centre
    mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
    filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
    plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
    means = [np.nanmean(x) for x in filtered_data]
    #Sdev=[np.std(x) for x in filtered_data]
    #plt.plot(BinsMid, means,"o", color='black', markersize=3)
    plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_st, color='g', linestyle='-')
    plt.axvline(x=Dof_lawson_c8_end, color='g', linestyle='-')
    plt.axhline(y=ActualDiameter, color='r', linestyle=':')
    plt.ylabel('Diameter (K07), μm')
    plt.xlabel('Distance from object plane, mm')
    plt.gca().set_ylim([0, 200])
    plt.gca().set_xlim([-39, 39]) 
    plt.gca().yaxis.grid(True)
    if Strip == 1 : 
        plt.axvspan(-50, -29.5, alpha=0.5, color='grey')
        plt.axvspan(15.5, 50, alpha=0.5, color='grey')
    
    majorLocator = MultipleLocator(10)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    plt.gca().xaxis.set_major_locator(majorLocator)
    plt.gca().xaxis.set_major_formatter(majorFormatter)
    plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
    
    #plt.savefig(Figurename,dpi=200)
    #plt.close(fig)   



#_______________________________________________________________________________________
 
    
 #AF_Level2_avg, AF_Level1_avg, AF_Level0_avg, D_Level2_avg, D_Level1_avg, D_Level0_avg, BinCentre    
    
def D25_D0(Model_Stats):
    

    i=0
    CentreOfDOF_Model = 50 
    for i in range(50,155,5):
        ActualDiameter=i
        Zd=DimensionlessDistance(Model_Stats[str(i)][6][:],ActualDiameter,CentreOfDOF_Model)
        D25=np.array(Model_Stats[str(i)][5][:])    
        D25_D0_ratio = D25 / ActualDiameter
        plt.plot(Zd,D25_D0_ratio,'o',markersize=2,linestyle='None',color='black')

#_________________________________________________________________________________________
    
# Remove particles without 67% pixels

    
def Model_67grayscale():
    
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255067/'
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(GSpath)
        
    
    Y_centre=50
    
    
    
    Model_D= np.linspace(50,150,21)
    #Output=[]
    
    #for D0 in range(50,155,5):
    for D0 in Model_D : 
        
        Lambda = 0.658
        Dof_lawson_c6= Dof_lawson(D0,Lambda,6)
        Dof_lawson_c6_st=-1*(Dof_lawson_c6/2)
        Dof_lawson_c6_end=(Dof_lawson_c6/2)   
        Dof_lawson_c48= Dof_lawson(D0,Lambda,4.8)
        Dof_lawson_c48_st=-1*(Dof_lawson_c48/2)
        Dof_lawson_c48_end=(Dof_lawson_c48/2) 
                
        #D0=50
        ParticleStatsTmp=ParticleStatsWv[np.where(ActualDiameterAll==D0),:] 
        ParticleStatsExp=ParticleStatsTmp[0,:,:]
        Particle_X=StageXAll[np.where(ActualDiameterAll==D0)]
        AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_X, ParticleStatsExp,0,255067)
        D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
        
        D_KorolevCorr=np.where(AreaFraction2>0,D_KorolevCorr,np.nan)
        DiameterLevel1=np.where(AreaFraction2>0,DiameterLevel1,np.nan)
    
    
        fig=plt.figure(figsize=(10,10))
        plt.rcParams.update({'font.size': 14})
        ax1=plt.subplot(2, 1, 1)
        Diameter_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(DiameterLevel1,Particle_X,0,100,100)       
        BinsMid-=Y_centre
        mask = ~np.isnan(Diameter_bin) # remove nans
        filtered_data = [d[m] for d, m in zip(Diameter_bin.T, mask.T)]
        plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
        means = [np.nanmean(x) for x in filtered_data]
        #plt.plot(BinsMid, means,"o", color='black', markersize=3)
        plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-',label='c=6')
        plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
        plt.axvline(x=Dof_lawson_c48_st, color='g', linestyle='-',label='c=4.8')
        plt.axvline(x=Dof_lawson_c48_end, color='g', linestyle='-')
        plt.ylabel('Diameter, μm')    
        plt.gca().yaxis.grid(True)
        plt.gca().set_ylim([0, 200])
        plt.legend()
        
        
        plt.subplot(2, 1, 2, sharex=ax1) 
        D_KorolevCorr_bin, BinsLower, BinsUpper, BinsMid=BinLinearV2(D_KorolevCorr,Particle_X,0,100,100)
        BinsMid-=Y_centre
        mask = ~np.isnan(D_KorolevCorr_bin) # remove nans
        filtered_data = [d[m] for d, m in zip(D_KorolevCorr_bin.T, mask.T)]
        plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=BinsMid)
        means = [np.nanmean(x) for x in filtered_data]
        #plt.plot(BinsMid, means,"o", color='black', markersize=3)
        plt.axvline(x=Dof_lawson_c6_st, color='b', linestyle='-')
        plt.axvline(x=Dof_lawson_c6_end, color='b', linestyle='-')
        plt.axvline(x=Dof_lawson_c48_st, color='g', linestyle='-')
        plt.axvline(x=Dof_lawson_c48_end, color='g', linestyle='-')
        plt.gca().set_ylim([0, 200])
        plt.gca().yaxis.grid(True)
        #plt.legend()
        plt.ylabel('Diameter (Korolev), μm')
        plt.xlabel('Distance from object plane, mm')
        
        plt.gca().set_xlim([-50, 50]) 
        plt.gca().yaxis.grid(True)
    
        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(5)
        plt.gca().xaxis.set_major_locator(majorLocator)
        plt.gca().xaxis.set_major_formatter(majorFormatter)
        plt.gca().xaxis.set_minor_locator(minorLocator)# for the minor ticks, use no labels; default NullFormatter
        
        
        Figurename= GSpath+str(int(D0))+'um.png'
        plt.savefig(Figurename,dpi=200)
        plt.close(fig) 


#_________________________________________________________________________________________

# Make size distribution of dropgun scan with and without 67%threshold

def Dropgun_67_pdf():
    StageFile='StagePositions_exp001.csv'
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255067/20180927120501/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255067/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D40_D0_Zd.hdf5' # just use this one not too important
    DLow_D0_model, Zd_model=Get_DLow_D0_model(DLow_D0_Zd_Table) 
    SizeBinsEdge=np.linspace(7.5,967.5,num=129)
    SizeBinsMid= (SizeBinsEdge[0:-1:1]+SizeBinsEdge[1::1])/2
    Lambda = 0.658	# in um 
    Lambda_mm=  Lambda / 1000
    
    Svol_noThresh=np.zeros(len(SizeBinsMid))
    Svol_67=np.zeros(len(SizeBinsMid))
    for i in range(len(SizeBinsMid)):
        #sVol_EAW_INT(TAS,D0,c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model)    
        Svol_noThresh[i]=sVol_EAW_INT(1,SizeBinsMid[i],8,64,15,Lambda_mm, 70, DLow_D0_model, Zd_model)
        Svol_67[i]=sVol_EAW_INT(1,SizeBinsMid[i],4.8,64,15,Lambda_mm, 70, DLow_D0_model, Zd_model)
    
    
    D_noNan = D_KorolevCorr[~np.isnan(D_KorolevCorr)] # remove nans
    #D_noNan = DiameterLevel1[~np.isnan(DiameterLevel1)] # remove nans
    
    D_hist, notused=np.histogram(D_noNan, bins=SizeBinsEdge) 
    #D_hist= D_hist/ np.nansum(D_hist)
    
    D_hist= D_hist/(Svol_noThresh*7.5*398) # 398s is length of scan
    #plt.plot(SizeBinsMid, D_hist,'o-')
    plt.bar(SizeBinsMid, D_hist, width=7.5, label='No threshold', color='grey')
    
    D_67= np.where(AreaFraction2 > 0, D_KorolevCorr, np.nan)
    #D_67= np.where(AreaFraction2 > 0, DiameterLevel1, np.nan)
    
    D_67_noNan = D_67[~np.isnan(D_67)] 
    D_hist, notused=np.histogram(D_67_noNan, bins=SizeBinsEdge) 
    #D_hist= D_hist/ np.nansum(D_hist)
    D_hist= D_hist/(Svol_67*7.5*398) # 398s is length of scan
    plt.plot(SizeBinsMid, D_hist,'ro-',label='67% threshold')
    plt.legend()
    
    plt.xlim([0,200])
    plt.xlabel('Diameter (K07), μm')
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    Figurename= StagePath+'Pdf_67threshold.png'
    plt.savefig(Figurename,dpi=200)
    #plt.close(fig) 
    
#_________________________________________________________________________________________

def BatchPlot_Dropgun_PSD():
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    #FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    #ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
    #Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
    
    
    fig=plt.figure(figsize=(7,10))
    plt.rcParams.update({'font.size': 12})
    #plt.subplot(2, 2, 1)
    #AF_to_Zd(Model_Stats,FilePath,0)
    
    
    plt.subplot(3, 1, 1)
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/20180927111105/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/'
    Dropgun_PSD(ExpDate,GSpath, StagePath, 1,'60 μm printhead')
    
    plt.subplot(3, 1, 2)
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/20180927134617/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/'
    Dropgun_PSD(ExpDate,GSpath, StagePath, 0,'90 μm printhead')
    
    plt.subplot(3, 1, 3)
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/Both/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/'
    Dropgun_PSD(ExpDate,GSpath, StagePath, 2,'120 μm printhead')
    
    Figurename= 'C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/DropgunPSD_255075.png'
    plt.tight_layout() 
    plt.savefig(Figurename,dpi=200)

#_________________________________________________________________________________________


# Make size distribution of dropgun scan

def Dropgun_PSD(ExpDate,GSpath, StagePath, LabelFlag, PltTitle):
    StageFile='StagePositions_exp001.csv'   
    
    
    Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
    D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)
    DLow_D0_Zd_Table= 'C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5' # 
    DLow_D0_model, Zd_model=Get_DLow_D0_model(DLow_D0_Zd_Table) 
    SizeBinsEdge=np.linspace(7.5,967.5,num=129)
    SizeBinsMid= (SizeBinsEdge[0:-1:1]+SizeBinsEdge[1::1])/2
    Lambda = 0.658	# in um 
    Lambda_mm=  Lambda / 1000
    
    #Svol_noThresh=np.zeros(len(SizeBinsMid))
    Svol=np.zeros(len(SizeBinsMid))
    #for i in range(len(SizeBinsMid)):
        #sVol_EAW_INT(TAS,D0,c,ArrayElements,ProbeRes,Lambda_mm, ArmSep, DLow_D0_model, Zd_model)    
   #     Svol_noThresh[i]=sVol_EAW_INT(1,SizeBinsMid[i],8,64,15,Lambda_mm, 70, DLow_D0_model, Zd_model)
         
     #means = [np.nanmean(x) for x in filtered_data]
    
    for i in range(4,9,1) :
        D= np.where(Zd_fromRatios < i, D_KorolevCorr, np.nan)
    #D_67= np.where(AreaFraction2 > 0, DiameterLevel1, np.nan)
        D_noNan = D[~np.isnan(D)] 
        D_hist, notused=np.histogram(D_noNan, bins=SizeBinsEdge) 
    #D_hist= D_hist/ np.nansum(D_hist)
        Svol = np.array([sVol_EAW_INT(1,SizeBinsMid[x],i,64,15,Lambda_mm, 70, DLow_D0_model, Zd_model) for x in range(len(SizeBinsMid)) ])
        D_hist= D_hist/(Svol*7.5*398)
        plt.plot(SizeBinsMid, D_hist,'-',label='Z$_{d}$ < '+str(i))
    
    
    D_noNan = D_KorolevCorr[~np.isnan(D_KorolevCorr)] # remove nans
    #D_noNan = DiameterLevel1[~np.isnan(DiameterLevel1)] # remove nans
    D_hist, notused=np.histogram(D_noNan, bins=SizeBinsEdge) 
    #D_hist= D_hist/ np.nansum(D_hist)
    Svol = np.array([sVol_EAW_INT(1,SizeBinsMid[x],8,64,15,Lambda_mm, 70, DLow_D0_model, Zd_model) for x in range(len(SizeBinsMid)) ])
    D_hist= D_hist/(Svol*7.5*398)
    #plt.plot(SizeBinsMid, D_hist,'-', label='No threshold', color='grey')
    plt.bar(SizeBinsMid, D_hist, width=7.5, label='No threshold', color='grey')
    
    
    if LabelFlag==1:
        plt.legend()
    plt.xlim([0,200])
    
    plt.ylabel('dN dDp$^{-1}$, L$^{-1}$ μm$^{-1}$')
    if LabelFlag == 2 :
        plt.xlabel('Diameter (K07), μm')
    
    plt.title(PltTitle)
    
    #plt.yscale('log')
    
    #Figurename= StagePath+'Dropgun_PSD.png'
    #plt.savefig(Figurename,dpi=200)
    #plt.close(fig) 
    
    
 #_________________________________________________________________________________________
# Do a histogram of D at Z= 0. Compare model and dropgun

def D_at_Z0_model_dg():

    StageFile='StagePositions_exp001.csv'  
    
    #***************** September 2018 ********************

    #60um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/20180927111105/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X_60, Particle_Y_60, PartGrayStatsWv_60, ParticleTime_60,ParticleTimesWv_60=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0_60,AreaFraction1_60,AreaFraction2_60,AreaRatio_2_1_60,AreaRatio_1_0_60,DiameterLevel0_60,DiameterLevel1_60,DiameterLevel2_60, FilledArea_60,VoidArea_60,DiameterBG_60,Zd_fromRatios_60=GrayscaleRatios(Particle_X_60, Particle_Y_60, PartGrayStatsWv_60,0,255075)


    #90um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/20180927134617/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/90um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X_90, Particle_Y_90, PartGrayStatsWv_90, ParticleTime_90,ParticleTimesWv_90=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0_90,AreaFraction1_90,AreaFraction2_90,AreaRatio_2_1_90,AreaRatio_1_0_90,DiameterLevel0_90,DiameterLevel1_90,DiameterLevel2_90, FilledArea_90,VoidArea_90,DiameterBG_90,Zd_fromRatios_90=GrayscaleRatios(Particle_X_90, Particle_Y_90, PartGrayStatsWv_90,0,255075)
 

    #120um_25_50_75
    GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/Both/Output/'
    StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/120um_255075/'
    ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
    Particle_X_120, Particle_Y_120, PartGrayStatsWv_120, ParticleTime_120,ParticleTimesWv_120=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
    AreaFraction0_120,AreaFraction1_120,AreaFraction2_120,AreaRatio_2_1_120,AreaRatio_1_0_120,DiameterLevel0_120,DiameterLevel1_120,DiameterLevel2_120, FilledArea_120,VoidArea_120,DiameterBG_120,Zd_fromRatios_120=GrayscaleRatios(Particle_X_120, Particle_Y_120, PartGrayStatsWv_120,0,255075)

    SavePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'
    
    
    ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(SavePath)
    
    for D0 in range(50,155,5): 
        #D0=50
        ParticleStatsTmp=ParticleStatsWv[np.where(ActualDiameterAll==D0),:] 
        ParticleStatsExp=ParticleStatsTmp[0,:,:]
        Particle_X_model=StageXAll[np.where(ActualDiameterAll==D0)]
        AreaFraction0_model,AreaFraction1_model,AreaFraction2_model,AreaRatio_2_1_model,AreaRatio_1_0_model,DiameterLevel0_model,DiameterLevel1_model,DiameterLevel2_model, FilledArea_model, VoidArea_model,DiameterBG_model,Zd_fromRatios_model=GrayscaleRatios(Particle_X_model, Particle_X_model, ParticleStatsExp,0,255075)
        Figurename=SavePath+'Diameter_Z0_'+str(D0)
        
        fig=plt.figure(figsize=(7,7))
        PlotZ0hist(DiameterLevel1_model, Particle_X_model,50,'Model '+str(D0)+ ' μm' , 1)
        PlotZ0hist(DiameterLevel1_60, Particle_X_60,29,'Dropgun 60 μm', 0)
        PlotZ0hist(DiameterLevel1_90, Particle_X_90,29,'Dropgun 90 μm', 0)
        PlotZ0hist(DiameterLevel1_120, Particle_X_120,29,'Dropgun 120 μm', 0)
        
        plt.xlabel('Diameter, μm')
        plt.ylabel('Normalised frequency')
        plt.legend()
        plt.xlim([0,180])
        loc = plticker.MultipleLocator(base=15)
        plt.gca().xaxis.set_major_locator(loc)
        plt.grid(which='major', axis='x', linestyle='-')
        
        Figurename=SavePath+'Diameter_Z0_'+str(D0)+'.png'
        plt.savefig(Figurename,dpi=200)
        plt.close(fig) 

 #_________________________________________________________________________________________


def PlotZ0hist(DiameterLevel1, Particle_X,X_centre, Label, Bar):
    
    #DiameterLevel1_Z0=DiameterLevel1[Particle_X<(X_centre-0.5) or Particle_X>(X_centre+0.5)]
    DiameterLevel1_Z0=DiameterLevel1[np.logical_and(Particle_X>(X_centre-0.1),Particle_X<(X_centre+0.1))]
    DiameterLevel1_noNan = DiameterLevel1_Z0[~np.isnan(DiameterLevel1_Z0)]
    SizeBinsEdge=np.linspace(7.5,967.5,num=65)
    SizeBinsMid= (SizeBinsEdge[0:-1:1]+SizeBinsEdge[1::1])/2
    #print(np.nansum(DiameterLevel1_Z0))
    
    if  not np.isnan(np.sum(DiameterLevel1_noNan)) :
        #print('test')
        D_hist, notused=np.histogram(DiameterLevel1_noNan, bins=SizeBinsEdge) 
        D_hist= D_hist / np.sum(D_hist)
        
        #fig=plt.figure(figsize=(7,7))
        if Bar == 1: 
            plt.bar(SizeBinsMid, D_hist, width=7.5, color='grey', label=Label)
        else :
            plt.plot(SizeBinsMid, D_hist, '-D', label=Label)    
        
        #plt.xlabel('Diameter, μm')
        #plt.ylabel('Counts')
        #plt.legend()
        #plt.xlim([0,180])
        
        #if SavePath != 0 :
        
        #Spacing between each line
        #intervals = float(sys.argv[1])

        #loc = plticker.MultipleLocator(base=15)
        #plt.gca().xaxis.set_major_locator(loc)
        #ax.yaxis.set_major_locator(loc)
        #plt.grid(which='major', axis='x', linestyle='-')
        

    
    