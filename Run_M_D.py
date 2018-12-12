# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 11:12:42 2018

@author: Admin TEMP
"""

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
FlightDict_C081=LoadFlightData2Dict(1,1,1,'C081',PathDict)
AddRadar2FlightDict(FlightDict_C081)
CSVPath='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/c081-feb-13/'
CSVName='flight-sum_faam_20180213_r0_c081.csv'


#Filename='20180213_0552_CAMRA_RHI.nc'
#Filename='20180213_0652_CAMRA_RHI.nc'
Filename='20180213_0734_CAMRA_RHI.nc'
RadarData, ZED_H_1D,Distance_1D,Altitude_1D, Range, Elevation=LoadRadar(Filename)
Distance2CHB, FlightTrack_ZED_H, LAT, LON, ALT,Time_Core,ZED_H_grid, DistanceBins, AltitudeBins = Colocate(FlightDict_C081,ZED_H_1D,Distance_1D,Altitude_1D,Filename.replace('.nc','.png'))

IWC_avg_colocate, IWC_avg, Zed_Alt_avg, Zed_1D_avg, Zed_colocate_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr= SingleRunPSD(7,FlightDict_C081,FlightTrack_ZED_H,AltitudeBins,ZED_H_grid,CSVPath,CSVName)
#Find_a_b(CSVPath,IWC_avg_colocate, Zed_1D_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr)

aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random =Find_a_b_MonteCarlo(CSVPath,IWC_avg_colocate, Zed_1D_avg,CompositeSize, Composite_dN, Composite_dNdDp,RunStr)
a_b_MonteCarlo_hist(aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random,RunStr)