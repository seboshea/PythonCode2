#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:56:34 2018

@author: mbexkes3
"""
from numpy import zeros, zeros_like, reshape, float32
from netCDF4 import Dataset
import os
import subprocess
import pickle
import copy

#import namelist as n
import constants as c
import setup_grids as s

import importlib # NOTE: need to be in the correct directory for this to work

def run():
  #  importlib.reload(n)
 #   importlib.reload(c)
    #nmodes = n.nmodes
    #nbins = n.nbins
  #  from namelist import nbins, nmodes, Mass_frac, P,output_filename, Simulation_type, Dlow, rkm, dt, PRESS1, PRESS2, Temp1, Temp2, Heterogeneous_freezing_criteria, alpha_crit, runtime, T, RH, sig, NAER, D_AER, w
    #import namelist as n
 #   import namelist_0107_run1_sv1 as n

    with open('test.pkl','rb') as pk:
        [nbins, nmodes, NAER, D_AER, sig, Mass_frac, RH, T, P,
        output_filename, simulation_type, runtime, Dlow,
        dt, PRESS1, PRESS2, Temp1, Temp2, Heterogeneous_freezing_criteria, alpha_crit,
         w, SV_flag, SV_MR, SV_MF, n_sv, semi_vols] = pickle.load(pk)
   
 #   SV_MF = sv_mass_frac
    # -------------------- create index values -----------------------------------
    IND1 = nbins*nmodes # index for mass
    IND2 = IND1*2 # index for number
    IND3 = IND1*3
    INDSV1 = IND3
    INDSV2 = INDSV1 + IND1*n_sv
    IPRESS = -3
    ITEMP = -2
    IRH = -1

    IPRESS_ICE = -3
    ITEMP_ICE = -2
    IRH_ICE = -1


    # ----------------------------------------------------------------------------

    # --------------------- initailizes arrays ------------------------------------ 
    t_final = runtime

    Y = zeros(nbins*nmodes*3+nbins*nmodes*n_sv+3+n_sv)
    Y_AER = zeros([int(t_final/dt),IND1]) # mass of aerosol bin
    YICE_AER = zeros_like(Y_AER) # mass of aerosol in ice bins

    Y[IPRESS] = P
    Y[ITEMP] = T
    Y[IRH] = RH

 #   nmodes = nmodes
  #  nbins = nbins
   # dt = dt

    rhobin = zeros(nmodes*nbins)
    kappabin = zeros(nmodes*nbins)   
    molwbin = zeros(nmodes*nbins) 

    k = zeros(nmodes)
    rhoa = zeros(nmodes)
    molw_aer = zeros(nmodes)

    Kappa = zeros([nmodes,nbins])# use this one in all calculations except setup grids

    output = zeros([int(t_final/dt),len(Y)])
    output_ice = zeros([int(t_final/dt),len(Y)])

    CDP_CONC_liq = zeros([int(t_final/dt),30])
    CDP_CONC_ice = zeros([int(t_final/dt),30])
    CDP_CONC_total = zeros([int(t_final/dt),30])

    dummy  = zeros_like(Y)
    dummy2 = zeros_like(Y)
    dummy4 = zeros_like(Y)

    ice_aer = zeros([int(t_final/dt),IND1])

    ACT_DROPS = zeros([int(t_final/dt),IND1])
    ncomps = len(Mass_frac)
    
    Mass_frac_aer = copy.deepcopy(Mass_frac)
    
    # make dictionaries of mass fractions #######################
    if SV_flag:
        sv_mass_frac = {}
        for i in semi_vols[:n_sv]:
   #         Mass_frac[i] = [SV_MF]*nmodes
            sv_mass_frac[i] = [SV_MF]*nmodes
    ############################################################
    
    # Add in new aerosol types from user ------------------------------------------
    if ncomps > len(c.aerosol_dict):
        for key in Mass_frac:
            if key not in c.aerosol_dict:
                if sum(n.Mass_frac[key][:]) > 0.0:
                    print('Aerosol type', key, 'not found')
                    print('please enter values for new aerosol type ')
                    new_aerosol_type = list(range(4))
                    print('Molecular weight in kg per mole ')
                    new_aerosol_type[0] = input()
                    print('Density kg/m^3 ')
                    new_aerosol_type[1] = input()
                    print('Kappa')
                    new_aerosol_type[3] = input()
                    #create new temporary dictionary entry 
                    c.aerosol_dict[key] = list(map(float, new_aerosol_type))
                    # add IN properties to new aerosol type
                    print('Is aerosol type',key,'an INP?')
                   
                    INP = input('yes / no   ')
                    if INP[0] == 'y' : 
                        is_INP = True
                    else:
                        is_INP = False
                    if is_INP:
                        new_INP_type = zeros(3)
                        print('Please provide ns values for ',key)
                        new_INP_type[0] = 2 # just use feldspar type at the moment
                        print('ns value 1')
                        new_INP_type[1] = input()
                        print('ns value 2')
                        new_INP_type[2] = input()
                        c.ns_dict[key] = list(map(float, new_INP_type))
        
                    print('do you want to save new aerosol type ',key,' for future use')
                    to_save = input('yes / no   ')
                    if to_save[0] == 'y': 
                        to_save = True
                    else:
                        to_save = False
                    if to_save:

                        with open("constants.py", "r") as in_file:
                            buf = in_file.readlines()
                     
                        with open("constants.py", "w") as out_file:
                           for line in buf:
                             if line == "aerosol_dict = {'ammonium sulphate': [132.14e-3, 1770, 3, 0.61],\n":
                                 line = line + "\'"+key+"\': ["+str(new_aerosol_type[0])+','+str(new_aerosol_type[1])+','+str(new_aerosol_type[2])+','+str(new_aerosol_type[3])+"],"
                             if line == "ns_dict = {'Kaolinite':[0, -0.8881, 12.9445, -0.8802, 231.383],\n":
                                 if is_INP:
                                     line = line + "\'"+key+"\': ["+str(new_INP_type[0])+','+str(new_INP_type[1])+','+str(new_INP_type[2])+"],"
                             
                             out_file.write(line)
        # -----------------------------------------------------------------------------
    # --------- catch error if sum of mass fractions is greater than 1 ------------
    total = zeros(nmodes)

  #  for mode in range(nmodes):
   #     for key in list(n.Mass_frac.keys())[:]:
    #        total[mode] += n.Mass_frac[key][mode]
     #   if total[mode] != 1:
      #      print('ERROR - Simulation Stopped')
       #     print('sum of mass fractions in mode', mode+1, 'does not equal 1')
        #    print(n.Mass_frac)
         #   ERROR_FLAG = True
          #  break
        #else:
         #   ERROR_FLAG = False
          #  continue

    ERROR_FLAG = False        
    # -----------------------------------------------------------------------------
    nu = zeros_like(k)
    nubin = zeros_like(rhobin)
 #   rho_sv = zeros_like(rhobin)
     # calculate average values for each mode
    print(Mass_frac)
    for mode in range(nmodes):
        total_moles = sum([Mass_frac[key][mode]*100/c.aerosol_dict[key][0] for key in Mass_frac.keys()]) # mole weighted kappa
       
        k[mode] = sum([(Mass_frac[key][mode]*100/c.aerosol_dict[key][0]/total_moles)*c.aerosol_dict[key][3] for key in Mass_frac.keys()]) # mole weighted kappa
        nu[mode] = sum([(Mass_frac[key][mode]*100/c.aerosol_dict[key][0]/total_moles)*c.aerosol_dict[key][2] for key in Mass_frac.keys()]) # mole weighted kappa

   #     k[mode] = sum([n.Mass_frac[key][mode]*c.aerosol_dict[key][3] for key in c.aerosol_dict.keys()])
        rhoa[mode] = 1/sum([Mass_frac[key][mode]/c.aerosol_dict[key][1] for key in Mass_frac.keys()])
        molw_aer[mode] = sum([Mass_frac[key][mode]*c.aerosol_dict[key][0] for key in Mass_frac.keys()])
        #rho_sv[mode] = 1/sum([n.sv_mass_frac[key][mode]*c.semi_vol_dict[key][1] for key in n.sv_mass_frac.keys()])
        
                                        
    for i in range(nmodes):
        rhobin[i*nbins:nbins+nbins*i] = rhoa[i]
        kappabin[i*nbins:nbins+nbins*i] = k[i]
        molwbin[i*nbins:nbins+nbins*i] = molw_aer[i]
        nubin[i*nbins:nbins+nbins*i] = nu[i]
        
    ############################### set-up grid ###################################
    rkm = 1
    GRID = s.setup_grid(rhobin, kappabin, rkm, Dlow, nbins, 
                        nmodes, RH, sig, NAER, D_AER, T, rhoa, k)

    Y[nbins*nmodes:IND2] = GRID[0] # number of aerosol in each bin
    Y[:IND1] = GRID[1]               # water mass in each bin
    Y_AER[0,:] = GRID[2]                  # mass of aerosol in each bin
    YICE_AER[0,:] = GRID[2]                 # mass of aerosol in ice bins 

    for mode in range(nmodes):
        ind1 = nbins*mode
        ind2 = ind1+nbins
        Kappa[mode] = sum([Mass_frac[key][mode]*Y_AER[0,ind1:ind2]/
                             c.aerosol_dict[key][1]*c.aerosol_dict[key][3] for key
                             in Mass_frac.keys()])/sum([Mass_frac[key][mode]
                             *Y_AER[0,ind1:ind2]/c.aerosol_dict[key][1] for key
                             in Mass_frac.keys()])

    Kappa = reshape(Kappa,[nbins*nmodes])

    # setup output file
    def setup_output():
        var_names = ['ice_total','liq_total','drop_total','temperature','pressure','RH','liquid_water_content','ice_water_content']
        bin_var_names = ['ice_number','liq_number','ice_mass','liq_mass','activated_drops']
        
        position = [pos for pos, char in enumerate(output_filename) if char == '/']
        try:
            short_filename = output_filename[position[-1]+1:]
        except:
            short_filename = output_filename
        if os.path.isfile(short_filename):
        
            subprocess.call(['rm '+short_filename],shell=True)
            
        output_file = Dataset(short_filename,'w',format='NETCDF4_CLASSIC')
        #output_file = Dataset('output13.nc','w',format='NETCDF4_CLASSIC')
        
        output_file.createDimension('time',len(output))
        output_file.createDimension('bins',nbins)
        output_file.createDimension('modes', nmodes)
        output_file.createDimension('CDP_bins',30)
        
        for variable_name in var_names:
            output_file.createVariable(variable_name,float32,dimensions=('time',))
        for variable_name in bin_var_names:
            output_file.createVariable(variable_name,float32,dimensions=('time','modes','bins'))
        
        output_file.createVariable('CDP_CONC_total',float32,dimensions=('time','CDP_bins'))
        return output_file
    return nbins, nmodes, simulation_type, PRESS1, PRESS2, w, rhobin, molwbin, nubin, Kappa,Temp1, Temp2, Y, ERROR_FLAG, Y_AER, CDP_CONC_liq, CDP_CONC_ice, CDP_CONC_total,YICE_AER, alpha_crit, ncomps, output, output_ice, dummy2, ACT_DROPS, setup_output, Heterogeneous_freezing_criteria, Mass_frac, Mass_frac_aer, SV_flag, SV_MF, SV_MR, n_sv


if __name__ == "__main__":
    run()
           











