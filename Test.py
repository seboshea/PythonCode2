# -*- coding: utf-8 -*-
"""
Created on Thu May 24 10:45:55 2018

@author: Admin TEMP
"""

import datetime
import numpy as np
from netCDF4 import Dataset
import math
import bisect
import h5py
import matplotlib.pyplot as plt 
import scipy
from scipy.optimize import fsolve
import scipy.optimize as opt
from scipy.spatial import distance
from scipy.spatial.distance import pdist, squareform
from MyFunctions import LoadOAP_nc, loadmat, TimeSeriesPlot, Matlab2PythonTimeArray, Matlab2PythonTime, Time2Secs,Load_HaloHolo_PSD

from Mass_Dimension import Minimise_Fb

from GrayscaleFunctions import LoadGrayScaleStats,GetIAT_vector, ParticlesStats2PSD, GrayscaleRatios, DvsSvol

from PICASSO_functions import Avg_dndDp_CountStats

from DropGun import DimensionlessDistance, LoadModelGrayScale, ModelGS_Stats, LoadGrayScaleStats_Stage
import pandas as pd
from dateutil import parser

#SEA_path='C:/Users/Admin TEMP/Documents/PICASSO/Flights/FAAM_Data/SEA/SEAmerge_c081/'
#SEA_filename='seaprobe_core_processed_20180213_c081_MO-Nev.csv'

#df = pd.read_csv(SEA_path+SEA_filename)
#TimeStartStr=df['TIMESTAMP']
#SEA_Time=[]

#SEA_Time=[parser.parse(TimeStartStr[x]) for x in range(len(TimeStartStr))]
#SEA_Time= [datetime.datetime.strptime(TimeStartStr[x], "%Y-%m-%d %H:%M:%S.%f") for x in range(len(TimeStartStr))]

#021_LWC=df['021_LWC']
#021_LWC=df['021_IWC']
#021_LWC=df['021_TWC']

#083_LWC=df['083_LWC']
#083_LWC=df['083_IWC']
#083_LWC=df['083_TWC']


DLow_D0_Zd_Table='C:/Users/Admin TEMP/Documents/DropletGun/D25_D0_Zd.hdf5'
DvsSvol(DLow_D0_Zd_Table)
