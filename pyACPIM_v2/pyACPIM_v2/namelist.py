#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

nbins	= 60
nmodes	= 1
ncomps = 9
n_sv = 1

NAER 	=  [1.6623e9, 	55e6, 		9200e6, 	0, 0]
D_AER 	=  [40.59e-9, 1500e-9,	222e-9, 	0, 0]
sigma	    =  [0.4904,	0.4,		0.43, 		0, 0]

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

semi_vols = ['SV01','SV02']                                         
SV_MF = 1.0
SV_MR = 2.58e-10
SV_flag = False

RH = 0.955741
T = 287.3133
P = 97677
w = 5
runtime = 20

# choose type of simulation
simulation_type = 'chamber' # parcel or chamber

# choose criteria for heteogeneous freezing
freezing_crit = 'RH>1' # activation, RH or critical mass of water
alpha = 70 # for critical mass of water for freezing criteria

Dlow=10e-9
rkm=1
dt = 1

#constants for chamber T and P fits
PRESS1 =414.22
PRESS2 = 0.0015178
Twall = T
Therm_chamber = 0.0015
Temp1 = 3.1033
Temp2 = 0.012514

output_file_name = 'output_sv1_run1_0107.nc'
