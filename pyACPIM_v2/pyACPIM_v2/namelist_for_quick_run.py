#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# ==============================================================================
# This file contains values of variables that will be used as input to the model
# ==============================================================================
# NOTES:
# 1. DO NOT change variable names
# 2. DO NOT change the format of variables



nbins	= 60 # number of size bins for particles
nmodes	= 1 # number of aerosol modes
n_sv = 1 # number of semi volatile compounds, this part of the code is still under development leave as 1

NAER 	=  [1.6623e9, 	55e6, 		9200e6, 	0, 0] # number of aerosol in each mode
D_AER 	=  [40.59e-9, 1500e-9,	222e-9, 	0, 0] # median diameter of aerosol in each mode
sigma	    =  [0.4904,	0.4,		0.43, 		0, 0] # width of the aerosol mode (ln(geometeric standard deviation))

# mass fractional composition of aerosol in each mode
Mass_frac = {'ammonium sulphate': [1.0,   0.2],
            'sea salt':           [0.0,   0.0],
            'sulphuric acid':     [0.0,   0.0],
            'fulvic acid':        [0.0,   0.0],
            'Kaolinite':          [0.0,   0.0],
            'Montmorinillite':    [0.0,   0.0],
            'Feldspar':           [0.0,   0.0],
            'Illite' :            [0.0,   0.0],
            'Bio' :               [0.0,   0.0],
            'test' :              [0.0,   0.0]}

######## this part of code is under developement #####
semi_vols = ['SV01','SV02']                          #               
SV_MF = 1.0											 #
SV_MR = 2.58e-10									 #
SV_flag = False # leave as False                     #
######################################################

RH = 0.955741 # initial relative humidity
T = 287.3133 # initial temperature in kelvin
P = 97677 # initial pressure in Pa
w = 5 # updraft velocity for parcel simulations, this is a constant value (ms-1)
runtime = 20 # total run time of the model simulation in seconds

# choose type of simulation
simulation_type = 'chamber' # parcel or chamber

# choose criteria for heteogeneous freezing
freezing_crit = 'RH>1' # activation, RH>1 or critical mass of water
alpha = 70 # for critical mass of water for freezing criteria, see Simpson et al (2018) ACP

Dlow=10e-9 # size of smallest aerosol size bin in meters

dt = 1 # model time step in seconds

#constants for chamber temperature and pressure profiles
PRESS1 =414.22
PRESS2 = 0.0015178
Temp1 = 3.1033
Temp2 = 0.012514

output_file_name = 'output.nc'
