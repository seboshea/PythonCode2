# -*- coding: utf-8 -*-
"""
Created on Fri Jul 13 10:12:19 2018

@author: Admin TEMP
"""


import numpy as np
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
from matplotlib.mlab import bivariate_normal
from matplotlib.colors import LogNorm
import math
from netCDF4 import Dataset
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import scipy.io as sio
from scipy import optimize
from PICASSO_functions import SingleRunPSD,Colocate,LoadRadar, DataPaths,LoadFlightData2Dict, AddRadar2FlightDict

#_____________________________________________________________________________________________________
#mass diameter relationship from brown and francis 1995

def BrownFrancis(PSD,MidSize):
    	
    MidMass=np.zeros(len(MidSize))
    for i in range(len(MidSize)):
        if(MidSize[i]<100):
            prefix = 4.82
            exponent = -13
            power = 3.0
        else:
            prefix = 7.38
            exponent = -11
            power = 1.9
        MidMass[i] = prefix * (MidSize[i]**power) * (10**exponent)
	
    PSD_mass= PSD * np.transpose(MidMass)      
    
    TotalMass = np.sum(PSD_mass,axis=1)
    
    return PSD_mass, TotalMass  



#_____________________________________________________________________________________________________

#Find a and b using fit and optimise methods

def Find_a_b(CSVPath,IWC_avg, Zed_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr):   
    
    
    b_array=np.linspace(1.5,2.5,num=20)
    
    # Find a b using optimize*********************************
    aOptimise, bOptimise =Minimise_Fb(b_array,IWC_avg,Zed_avg,CompositeSize, Composite_dN)

    print('Optimise a = '+str(aOptimise))
    print('Optimise b = '+str(bOptimise))
    
    
    # Find a b using fit to b vs log(a) *********************************
    aFit, bFit = Find_a_b_fit(b_array,CSVPath,IWC_avg, Zed_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr)
    
    print('Fit a = '+str(aFit))
    print('Fit b = '+str(bFit)) 
    
#_____________________________________________________________________________________________________

# Do a monte carlo simulation varying IWC and Z calculating a and b using fitting and optimize methods


def Find_a_b_MonteCarlo(CSVPath,IWC_avg, Zed_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr):
    
    b_array=np.linspace(1.5,2.5,num=20)
    
    #Zed_avg= Zed_Alt_avg # choose which Z to use
    
    IWCsigma=IWC_avg*0.05
    Zsigma=1

    Npts=1000

    IWC_random = np.random.normal(loc=IWC_avg, scale=IWCsigma, size=Npts)
    Z_random = np.random.normal(loc=Zed_avg, scale=Zsigma, size=Npts)
    aFit_array = np.zeros(Npts)*np.nan
    bFit_array = np.zeros(Npts)*np.nan
    aOpt_array = np.zeros(Npts)*np.nan
    bOpt_array = np.zeros(Npts)*np.nan
        
    for i in range(0,Npts,1):
        aFit_array[i], bFit_array[i] = Find_a_b_fit(b_array,CSVPath,IWC_random[i], Z_random[i], CompositeSize, Composite_dN, Composite_dNdDp,RunStr)        
        aOpt_array[i], bOpt_array[i] = Minimise_Fb(b_array,IWC_random[i],Z_random[i],CompositeSize, Composite_dN)
      
    return aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random  
        
    print('test')


#_____________________________________________________________________________________________________


def a_b_MonteCarlo_hist(aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random,RunStr):
    
    Figurename=RunStr+'_a_b_sensitivity_1000'
    SavePath= 'C:/Users/Admin TEMP/Documents/PICASSO/Flights/ProcessedData/C081/'
    
    b_bins=np.linspace(-0.1,5,num=100)
    a_bins=np.logspace(-20,-5,num=100,base=10)

    fig=plt.figure(figsize=(10,10))
    plt.subplot(3, 2, 1)
    plt.hist(IWC_random, bins=100)
    plt.xlabel('IWC, g m$^{-3}$')
    plt.ylabel('Counts')
    #plt.text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)

    plt.subplot(3, 2, 2)
    plt.hist(Z_random, bins=100)
    plt.xlabel('Reflectivity, dBZ')
    plt.ylabel('Counts')

    ax1=plt.subplot(3, 2, 3)
    plt.hist(aOpt_array, bins=a_bins)
    plt.xlabel('a Optmisation')
    plt.xscale('log')
    plt.ylabel('Counts')
    #plt.title('a='+str(np.nanmean(aOpt_array)))
    plt.text(0.5, 0.9,'a='+str(np.nanmean(aOpt_array)), ha='center', va='center', transform=ax1.transAxes)
    
    
    ax2=plt.subplot(3, 2, 4)
    plt.hist(aFit_array, bins=a_bins)
    plt.xlabel('a Fit')
    plt.xscale('log')
    plt.ylabel('Counts')
    #plt.title('a='+str(np.nanmean(aFit_array)))
    plt.text(0.5, 0.9,'a='+str(np.nanmean(aFit_array)), ha='center', va='center', transform=ax2.transAxes)
    
    ax3=plt.subplot(3, 2, 5)
    plt.hist(bOpt_array, bins=b_bins)
    plt.xlabel('b Optmisation')
    plt.ylabel('Counts')
    #plt.title('b='+str(np.nanmean(bOpt_array)))
    plt.text(0.5, 0.9,'b='+str(np.nanmean(bOpt_array)), ha='center', va='center', transform=ax3.transAxes)
    
    ax4=plt.subplot(3, 2, 6)
    plt.hist(bFit_array, bins=b_bins)
    plt.xlabel('b Fit')
    plt.ylabel('Counts')
    #plt.title('b='+str(np.nanmean(bFit_array)))
    plt.text(0.5, 0.9,'b='+str(np.nanmean(bFit_array)), ha='center', va='center', transform=ax4.transAxes)
    
    plt.savefig(SavePath+Figurename,dpi=200)
    plt.close(fig)
    

#_____________________________________________________________________________________________________

# Use a fit to b vs log(a) for IWC and Z to estimate a and b that satify IWC and Z equations


def Find_a_b_fit(b_array,CSVPath,IWC_avg, Zed_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr):
    
    # Find a b using fit to b vs log(a)*********************************
    
    #b_array=np.linspace(1.5,2.5,num=20)
    #a_calcIWC,a_calcZ=Find_a_given_b(b_array,0.054926273720117796,-0.1149671037856154,CompositeSize, Composite_dN) # C081 run 2 with brown and francis
    #a_calcIWC,a_calcZ=Find_a_given_b(b_array,IWC_avg,Radar_ref,CompositeSize, Composite_dN) 
    
    a_calcIWC,a_calcZ=Find_a_given_b(b_array,IWC_avg,Zed_avg,CompositeSize, Composite_dN) 
    a_calcIWC_log=np.log10(a_calcIWC)
    mIWC,cIWC = np.polyfit(a_calcIWC_log, b_array, 1)
    a_calcIWC_fit=10**((b_array-cIWC)/mIWC)
    a_calcZ_log=np.log10(a_calcZ)
    mZ,cZ = np.polyfit(a_calcZ_log, b_array, 1)
    a_calcZ_fit=10**((b_array-cZ)/mZ)
    aResult=10**((cZ-cIWC)/(mIWC-mZ))
    bResult=mZ*math.log(aResult,10)+ cZ
    

       
    NumberofPlots=0
    
    if NumberofPlots>0 :
    
        b_array,a_array, IWC_sensitivity, Z_sensitivity=IWC_Z_sensitivity_2_ab(CompositeSize, Composite_dN)
        
        plt.figure(figsize=(10,10)) 
        plt.rcParams.update({'font.size': 16})
        plt.subplot(NumberofPlots, 1, 1)
        plt.title(RunStr)
    
        plt.pcolor(a_array,b_array,IWC_sensitivity, cmap='jet', norm=LogNorm())
        plt.plot(a_calcIWC,b_array,'+',color='black')
        plt.plot(a_calcIWC_fit,b_array,color='black')
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlim([1E-11, 1E-9])
        plt.ylabel('b')
        plt.xlabel('a')
        cbar=plt.colorbar(orientation='vertical')
        cbar.set_label('IWC, g m$^{-3}$')
    
        plt.subplot(NumberofPlots, 1, 2)
        plt.pcolor(a_array,b_array, Z_sensitivity, cmap='jet')
        plt.plot(a_calcZ,b_array,'+',color='black')
        plt.plot(a_calcZ_fit,b_array,color='black')
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlim([1E-11, 1E-9])
        plt.ylabel('b')
        plt.xlabel('a')
        cbar=plt.colorbar(orientation='vertical')
        cbar.set_label('Reflectivity, dBZ')
        
        if NumberofPlots==3 : 
            plt.subplot(NumberofPlots, 1, 3)
            plt.plot(a_calcIWC_fit,b_array,color='black')
            plt.plot(a_calcZ_fit,b_array,color='black')
            plt.xscale('log')
            #plt.yscale('log')
            plt.xlim([1E-11, 1E-9])
            plt.ylabel('b')
            plt.xlabel('a')
            ResultStr='a='+str(aResult)+', b='+str(bResult)
            plt.text(1E-9,2.5,ResultStr,fontsize=12,ha='right',va='top')

        plt.savefig(CSVPath+RunStr,dpi=200)
        
        plt.show()
    
    return aResult, bResult
    

    #return b_array,a_array, IWC_sensitivity, Z_sensitivity,a_calcIWC,a_calcZ, CompositeSize, Composite_dN

#_____________________________________________________________________________________________________

#Use PSD to create 2D array of IWC and Z using different a and b values.
    
    
 # brown and francis =  7.38E-11, 1. 
 
 

def IWC_Z_sensitivity_2_ab(CompositeSize, Composite_dN):    

    #Vary b = 1 to 3, increment 0.1
    # a =5E-12 to 5E-10, increment 0.25E-12

    #b_array=np.linspace(1,3,num=20)
    b_array=np.linspace(1.5,2.5,num=20)
    a_array=np.logspace(-11,-9,num=40,base=10)


    IWC_sensitivity=np.zeros((len(b_array),len(a_array)))
    Z_sensitivity=np.zeros((len(b_array),len(a_array)))

    for i in range(len(b_array)):
        for j in range (len(a_array)):
            Radar_ref, IWC = CalculateReflectivityIce(CompositeSize, Composite_dN, a_array[j], b_array[i])
            IWC_sensitivity[i,j]=IWC
            Z_sensitivity[i,j]=Radar_ref

    return b_array,a_array, IWC_sensitivity, Z_sensitivity

#_____________________________________________________________________________________________________
    
#_____________________________________________________________________________________________________

# Calculate ice reflectivity


# brown and francis a= 7.38E-11, b = 1.9, D in um, M in g
#Eq 2 Hogan et al., 2006  


def CalculateReflectivityIce(Size, dN_L, a, b):
    # dN in L-1
    dN=dN_L*1000 # m-3
    density_ice= 0.9167/1000			# g mm^-2    
    radar_scaleF = ( (0.174/0.93) * (36/(math.pi*math.pi*density_ice*density_ice)) )  # g^-2 mm^6	   
    #Size/= 1E6 # m
    #SecondMoment= (dN)*Size**2
    #ForthMoment= (dN)*Size**4
    #Size*=1E6 #um
    massVal= np.zeros(len(Size))
    for i in range(len(Size)):		
        if(Size[i] >= 100) : 
            massVal[i]=a*(Size[i]**b)			# g
        else : 
            massVal[i]=(4.82E-13)*(Size[i]**3)		# g
    Radar_ref_array= dN * massVal * massVal		# g^2 m^-3     
    Radar_ref= 10* math.log (radar_scaleF * np.sum(Radar_ref_array),10)		# Log to the base 10
    IWC= np.sum(dN* massVal) # g m3
    #print(Radar_ref) 
    #print(IWC)    
    return Radar_ref, IWC


#_____________________________________________________________________________________________________

# Calculate reflecivtiy using drop PSD


#Mason et al 2017 eq7     

def CalculateReflectivityLiquid(Size_um, dN_L, a, b):
    # dN in L-1
    dN=dN_L*1000 # m-3
    Size= Size_um/1000 # mm
   
    
    #Size/= 1E6 # m
    #SecondMoment= (dN)*Size**2
    #ForthMoment= (dN)*Size**4
    #Size*=1E6 #um
    
    MRratio=1 #Mie–Rayleigh backscatter ratio at the radar frequency
    Radar_ref_array= np.zeros(len(Size))
    for i in range(len(Size)):		
            Radar_ref_array[i]=dN*(Size[i]**6)*MRratio  # mm6/m3
    
    Radar_ref= 10* math.log (np.sum(Radar_ref_array),10)		# Log to the base 10
    
    #print(Radar_ref) 
    #print(IWC)    
    return Radar_ref



#_____________________________________________________________________________________________________

# Calculate a using b values in b_array


def Find_a_given_b(b_array,IWC,Z,CompositeSize, Composite_dN):
    
    #plt.plot(CompositeSize, Composite_dN)
    #plt.xscale('log')
    #plt.yscale('log')
    
    dN=Composite_dN*1000 # m-3
    density_ice= 0.9167/1000			# g mm^-2    
    radar_scaleF = ( (0.174/0.93) * (36/(math.pi*math.pi*density_ice*density_ice)) )  # g^-2 mm^6
    
    Reflectivity=10**(Z/10)
    
    lessThan100=np.zeros(len(CompositeSize))
    moreThan100=np.zeros(len(CompositeSize))
    
    a_calcIWC=np.zeros(len(b_array))
    a_calcZ=np.zeros(len(b_array))
    
    for j in range(len(b_array)):
        b=b_array[j]
        #Tmp=0
        for i in range(len(CompositeSize)):    
            if(CompositeSize[i] >= 100) : 
                moreThan100[i]=(CompositeSize[i]**b)
            else : 
                lessThan100[i]=(4.82E-13)*(CompositeSize[i]**3)
                #massVal[i]=(4.82E-13)*(Size[i]**3)		
        a_calcIWC[j]=(IWC-np.sum(lessThan100*dN))/np.sum(moreThan100*dN)
        a_calcZ[j]=np.sqrt(((Reflectivity/radar_scaleF)-np.sum(dN*lessThan100*lessThan100))/np.sum(dN*moreThan100*moreThan100))
         
    return a_calcIWC, a_calcZ
    

#____________________________________________________________________________________________________

# Find minimum in F_b. 


def Minimise_Fb(b_array,IWC,Zed,CompositeSize, Composite_dN):

    # Do optimisation to find b
    F_b_array = np.zeros(len(b_array))   
    for i in range(len(b_array)): 
        F_b_array[i]=F_b(b_array[i],IWC,Zed,CompositeSize, Composite_dN)
    
    #plt.plot(b_array,F_b_array)
    #plt.xlabel('b')
    #plt.ylabel('F(b)')
    
    b=2
    b_opt=optimize.minimize(F_b, b, args=(IWC,Zed,CompositeSize, Composite_dN),method='Nelder-Mead')
    b=b_opt['x']
    #print('Optimise b = '+str(b))
    
    # Use b to find a
    dN_m3=Composite_dN*1000 # m-3
    TempD=(np.where(CompositeSize>100,(CompositeSize**b)*dN_m3, np.nan))
    Part1=np.nansum(TempD)
    TempD =(np.where(CompositeSize<100, (4E-13*CompositeSize**3)*dN_m3,np.nan))
    Part2 = np.nansum(TempD)
    a = (IWC - Part2)/ Part1
    #print('Optimise a = '+str(a))    
    
    return a, b 
    
    
#____________________________________________________________________________________________________

# Rearrange IWC and Z mass dimension equations. To the form F(b)=0


def F_b(b, IWC,Zed,CompositeSize, Composite_dN):
    
    Reflectivity=10**(Zed/10)
    dN_m3=Composite_dN*1000 # m-3
    density_ice= 0.9167/1000			# g mm^-2    
    radar_scaleF = ( (0.174/0.93) * (36/(math.pi*math.pi*density_ice*density_ice)) ) # g^-2 mm^6
         
    TempD=(np.where(CompositeSize>100,(CompositeSize**b)*dN_m3 , np.nan))
    Part1 = np.nansum(TempD)

    Part2 = Reflectivity / radar_scaleF
    
    TempD=(np.where(CompositeSize<100,(4E-13*CompositeSize**3)**2*dN_m3, np.nan ))
    Part3=np.nansum(TempD)
    
    TempD= (np.where(CompositeSize>100,(CompositeSize**(2*b))*dN_m3, np.nan))
    Part4=np.nansum(TempD)
    
    TempD=(np.where(CompositeSize<100,(4E-13*CompositeSize**3)*dN_m3,np.nan))
    
    Part5 = np.nansum(TempD) - IWC
 
    return np.absolute(Part1*math.sqrt((Part2-Part3)/Part4) + Part5)


#____________________________________________________________________________________________________

# Calculate M_D for all runs during flight
 
def Calculate_a_b_runs():  

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


    for i in range(1,14,1):    
        IWC_avg_colocate, IWC_avg, Zed_Alt_avg, Zed_1D_avg, Zed_colocate_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr= SingleRunPSD(i,FlightDict_C081,FlightTrack_ZED_H,AltitudeBins,ZED_H_grid,CSVPath,CSVName)
        
        
        if ((IWC_avg==IWC_avg) & (Zed_1D_avg==Zed_1D_avg)):
            #Find_a_b(CSVPath,IWC_avg_colocate, Zed_1D_avg, CompositeSize, Composite_dN, Composite_dNdDp,RunStr)
            aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random =Find_a_b_MonteCarlo(CSVPath,IWC_avg_colocate, Zed_1D_avg,CompositeSize, Composite_dN, Composite_dNdDp,RunStr)
            a_b_MonteCarlo_hist(aFit_array, bFit_array, aOpt_array, bOpt_array, IWC_random, Z_random,RunStr)


