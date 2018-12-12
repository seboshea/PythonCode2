# -*- coding: utf-8 -*-
"""
Created on Wed May 16 10:00:46 2018

@author: Admin TEMP
"""
import numpy as np
import matplotlib.pyplot as plt  
from DropGun import AllAreaFractionZd_2,LoadModelGrayScale, ModelGS_Stats, ModelGS_Stats_Plot,DropGun_Stats_All
from DropGun import AF_to_Zd,CompareModelDropgun,PlotDiameterAF,AllAreaFractionZd,DropGunModel_AF_Zd,LoadGrayScaleStats_Stage
from DropGun import Model_Dropgun_Diameter_6panel,KorolevCorr_vs_Zd,Level2_level1_2_Zd, Dof_lawson,ReuterBakanDiameter,Model_Diameter_Percentiles
from DropGun import DropGunModel_AF_Zd_2, PLot_Areafraction_6Panel_255075, D25_D0

from GrayscaleFunctions import GrayscaleRatios, LoadGrayScaleStats, ParticlesStats2PSD, DLow_D0_Zd_Table

from PICASSO_functions import Avg_dndDp

from MyFunctions import KorolevCorrectedD

import datetime

#DropGun_Stats=DropGun_Stats_All()

#PlotDiameterAF()

#AF_2_Zd(Model_Stats,1)


#FilePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/ModelCIP_255075/'

#ParticleStatsWv,ActualDiameterAll,StageXAll=LoadModelGrayScale(FilePath)
#Model_Stats=ModelGS_Stats(ParticleStatsWv,ActualDiameterAll,StageXAll)
#DLow_D0_Zd_Table(Model_Stats, 1, 'D40_D0_Zd.hdf5')


#Model_D_Percentiles,Model_D=Model_Diameter_Percentiles(ParticleStatsWv,ActualDiameterAll,StageXAll)

#PLot_Areafraction_6Panel(DropGun_Stats,Model_Stats)
#p=Level2_level1_2_Zd(Model_Stats,1)

#DropGunModel_AF_Zd(DropGun_Stats,Model_Stats)

#DropGunModel_AF_Zd_2(DropGun_Stats,Model_Stats)


#ParticleStatsWv,ActualDiameterAll,StageXAll=ModelGS_Stats_Plot(Model_Stats)

#CompareModelDropgun(DropGun_Stats,Model_Stats)
#PlotDiameterAF()

#AllAreaFractionZd(Model_Stats,1)
#AllAreaFractionZd_2(Model_Stats,1)

#p= AF_to_Zd(Model_Stats,FilePath,1)

#DropGunModel_AF_Zd(DropGun_Stats,Model_Stats)

#KorolevCorr_vs_Zd(Model_Stats,0)


#Model_Dropgun_Diameter_6panel()

#D25_D0_means,  BinsMid= D25_D0_Zd_Table(Model_Stats, 1)


#StageFile='StagePositions_exp001.csv'
#GSpath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255067/20180927120501/Output/'
#StagePath='C:/Users/Admin TEMP/Documents/DropletGun/CIPscan/September2018/60um_255067/'
#ExpDate= datetime.datetime(2018, 9, 27, 0, 0, 0)
#Particle_X, Particle_Y, PartGrayStatsWv, ParticleTime,ParticleTimesWv=LoadGrayScaleStats_Stage(GSpath,ExpDate,StagePath,StageFile)
#AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2,FilledArea,VoidArea,DiameterBG,Zd_fromRatios=GrayscaleRatios(Particle_X, Particle_Y, PartGrayStatsWv,0,255075)
#D_KorolevCorr=KorolevCorrectedD(FilledArea, VoidArea,DiameterLevel1)

    
SizeBinsEdge=np.linspace(7.5,967.5,num=130)
SizeBinsMid= (SizeBinsEdge[0:-1:1]+SizeBinsEdge[1::1])/2

D_noNan = D_KorolevCorr[~np.isnan(D_KorolevCorr)] # remove nans
#D_noNan = DiameterLevel1[~np.isnan(DiameterLevel1)] # remove nans

D_hist, notused=np.histogram(D_noNan, bins=SizeBinsEdge) 
D_hist= D_hist/ np.nansum(D_noNan)
#plt.plot(SizeBinsMid, D_hist,'o-')
plt.bar(SizeBinsMid, D_hist, width=7.5, label='No threhold')

D_67= np.where(AreaFraction2 > 0, D_KorolevCorr, np.nan)
#D_67= np.where(AreaFraction2 > 0, DiameterLevel1, np.nan)

D_67_noNan = D_67[~np.isnan(D_67)] 
D_hist, notused=np.histogram(D_67_noNan, bins=SizeBinsEdge) 
D_hist= D_hist/ np.nansum(D_67_noNan)
plt.plot(SizeBinsMid, D_hist,'ro-',label='67% threhold')
plt.legend()

plt.xlim([0,200])