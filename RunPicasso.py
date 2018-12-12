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
from MyFunctions import is_number, BinLinearV2
import os
from numpy import genfromtxt
import bisect
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import datetime

from MyFunctions import LoadCoreCloud,loadFAAMCore,LoadNevzorov, Average_nPts_2D
from MyFunctions import Average_nPts,Average_nPts_datetime,DateTime2IgorTime,ChangeTimeBaseAvg,ChangeTimeBase2DAvg

from PICASSO_functions import PlotCrossSections2CHB,PlotCrossSections, LoadFlightData2Dict, DataPaths,PlotFlightSummary, AverageFlightDict,CompareIWC,PlotFlight,BatchPlotCDPCore,TrackPlot
from PICASSO_functions import AddRadar2FlightDict,Plot_TvsHI, LoadFlightSummaryCSV, PSD_runs_dict,CompositePSD, SingleRunPSD, CalculateMVD, LoadRadar, Colocate, PlotRadar, BinAverage2D
 
from Mass_Dimension import a_b_MonteCarlo_hist,Find_a_b_MonteCarlo,Find_a_b,BrownFrancis

PathDict=DataPaths()

#FlightList=['C072','C073','C074','C075', 'C076', 'C077','C078','C079','C080','C081','C082','C097','C098']
#FlightList=['C074','C075', 'C076', 'C077','C078','C079','C080','C081','C082','C097','C098']
#FlightList=['C079','C080','C081','C082','C097','C098']

#FlightList=['C082']
#PlotFlight(PathDict,FlightList)


FlightDict_C081=LoadFlightData2Dict(1,0,1,0,'C081',PathDict)
#FlightDict_C098=LoadFlightData2Dict(0,1,0,1,'C098',PathDict)
#AddRadar2FlightDict(FlightDict_C081)

#CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c098-apr-24/'
#CSVName='flight-sum_faam_20180424_r0_c098.csv'
#FlightStr='C098'

#CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
#CSVName='flight-sum_faam_20180213_r0_c081.csv'
#FlightStr='C081'

#PSD_runs_dict(PathDict, FlightStr,CSVPath,CSVName)

#PlotCrossSections2CHB(FlightDict_C081,PathDict['C081','FullPath'],'C081')

#BatchPlotCDPCore()

#TrackPlot(FlightDict_C081,PathDict['C081','FullPath'], 'CO81')






#TrackPlot(FlightDict_C081,PathDict['C081','FullPath'], 'CO81')


#MVD=CalculateMVD(NC_All_accept_2DS_psd,NC_2DS_size)
#PSD_runs_dict(PathDict)
#CompareIWC(FlightDict_C082) 
#FlightDict_C081_10=AverageFlightDict(FlightDict_C081,10)
#PlotRadar(ZED_H_1D,Distance_1D,Altitude_1D)

#SavePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'

#CompareIWC(FlightDict_C081_60,'C081',CSVPath)

#NC_All_accept_total,NC_HI_total,NC_MI_total,NC_LI_total,NC_all_x, NC_all_y, NC_all_z, NC_S_z, NC_MI_z, NC_HI_z, NC_All_accept_CH0_z, NC_All_accept_CH1_z, NC_DateTime=BrownFrancis(PathDict,'C082')
#Time_Core=FlightDict_C082['Time_Core']
#NC_All_accept_CH0_core=ChangeTimeBase2DAvg(NC_DateTime,NC_All_accept_CH0_z,Time_Core, 1)


#FilePath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
#FileName='flight-sum_faam_20180213_r0_c081.csv'
#StartTime=LoadFlightSummaryCSV(FilePath,FileName)

#Plot_TvsHI(FlightDict_C081,PathDict['C081','FilePathOAP'],'C081_HI_IWC')




