# -*- coding: utf-8 -*-
"""
Created on Mon May 13 15:57:38 2019

@author: Admin TEMP
"""

import numpy as np
import matplotlib.pyplot as plt 
#from skimage.draw import ellipse
from skimage.measure import label, regionprops
#from skimage.transform import rotate
import h5py
import os
from DropGun import Dof_lawson
from GrayscaleFunctions import GetZd
from scipy import ndimage
from skimage import filters
from MyFunctions import BinLinearV2, KorolevCorrectedD, BinLinear_median
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
import bisect
import matplotlib.animation as animation
#import cv2


#_______________________________________________________________________________________

# Compute diffraction intensity (I) at a given Z

def compute_diffraction(Z, l, pixel_size, x, y, X, Y, M, filter_fmax):

    # Parameters
    k = 2*np.pi/l
    
    delta_fx = 1./(x.size*pixel_size) # In order to be able to use FFT, Eq. (4.4-12) of Ersoy (2006)
    delta_fy = 1./(y.size*pixel_size)
    fx_max = x.size/2*delta_fx # µm-1
    fy_max = y.size/2*delta_fy # µm-1
    fx = np.arange(-fx_max, fx_max-delta_fx+.0000001, delta_fx)
    fy = np.arange(-fy_max, fy_max-delta_fy+.0000001, delta_fy)

    FX, FY = np.meshgrid(fx, fy)

    # Field at Z = 0
    U0 = np.copy(M)

    # Angular spectrum at Z = 0 (Fourier tranform), Eq. (4.3-2) of Ersoy (2006) 
    A0 = np.fft.fftshift(np.fft.fft2(U0))
    # fft2: Compute the 2-dimensional discrete Fourier Transform
    # fftshift: Shift the zero-frequency component to the center of the spectrum

    # Transfer function, Eq. (4.3-14) of Ersoy (2006)
    non_homogeneous_waves = 4*np.pi**2*(FX**2+FY**2) > k**2 # Eq. (4.3-13) of Ersoy (2006)
    H = np.exp( 1j*Z*np.abs(np.sqrt(k**2-4*np.pi**2*(FX**2+FY**2))) )
    H[non_homogeneous_waves] = 0

    # Angular spectrum at Z = Z, Eq. (4.3-11) of Ersoy (2006)
    AZ = A0*H

    # Lowpass filter to remove noise
    if False:
        filtre = np.zeros((y.size, x.size))
        A = (np.abs(FY) <= filter_fmax)
        B = (np.abs(FX) <= filter_fmax)
        filtre[A&B] = 1
        #if False:
        #    plot_image(y_b, x_b, filtre)
        AZ = AZ*filtre

    # Field at Z = Z, Eq. (4.3-12) of Ersoy (2006)
    UZ = np.fft.ifft2(np.fft.fftshift(AZ))

    # Intensity at Z = Z
    I = np.abs(UZ)**2

    # Phase at Z = Z
    P = np.arctan(np.imag(UZ)/np.real(UZ))
    
    return I, A0, fx, fy 




#_______________________________________________________________________________________

def plot_diffraction(xOAP,yOAP, ImageLow, ImageMid, ImageHigh,x, y, imageZ0, imageDiffraction,A0, fx, fy, Z, Zd,SaveFlag, Figurename):
    
    
    imageGreyscale = (ImageLow - ImageMid )*25 + (ImageMid - ImageHigh )*50 + (ImageHigh )*75
       
    
    fig=plt.figure(figsize=(8,8)) 
    plt.rcParams.update({'font.size': 12})
#    plt.subplot(3, 1, 1)
#    plt.title('z=0')
#    plt1=plt.pcolormesh(x,y,imageZ0, cmap='Greys_r')
#    plt1.set_clim(vmin=0, vmax=1)
#    plt.ylabel('y, μm')
#    cbar=plt.colorbar(orientation='vertical')
#    plt.ylim([-200,200])
#    plt.xlim([-200,200])
#    
#    plt.subplot(3, 1, 2)
#    plt.title('z='+str(Z) + ' μm, Zd = '+str(np.around(Zd,2)))
#    plt2=plt.pcolormesh(x,y,imageDiffraction, cmap='Greys_r', vmin=0, vmax = 1.5)
#    plt.ylabel('y, μm')
#    cbar=plt.colorbar(orientation='vertical')
#    plt.ylim([-200,200])
#    plt.xlim([-200,200])
    
    #imageDiffractionThreshold= (np.where(imageDiffraction>0.5, 1, 0))
#    plt.subplot(3, 1, 3)
    plt.title('z='+str(Z/1000) + ' mm, Zd = '+str(np.around(Zd,2)))
    plt3= plt.pcolormesh(xOAP,yOAP,imageGreyscale, cmap='gist_stern_r')
    plt3.set_clim(vmin=0, vmax=75)
    plt.ylabel('y, μm')
    plt.xlabel('x, μm')
    plt.ylim([-200,200])
    plt.xlim([-200,200])
    
    if SaveFlag == 1 :
        #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
        plt.savefig(Figurename,dpi=200)
        plt.close(fig)
        
    
#    plt.subplot(2, 2, 4)
#    plt.title('A0')
#    plt.pcolormesh(fx, fy,np.abs(A0), cmap='Greys_r')
#    plt.xlabel('fx')
#    plt.ylabel('fy')
#    cbar=plt.colorbar(orientation='vertical')
    #plt.ylim([-2000,2000])
    #plt.xlim([-2000,2000])
    
    #cbar.set_label('IWC, g m$^{-3}$') 

#_______________________________________________________________________________________

#Average 1 um pixels to 10 um pixels. 

def Average2OAPpixels(I, x, y) : 

    I_subset=I[int(I.shape[0]/2)-500: int(I.shape[0]/2) + 500, int(I.shape[1]/2)-500: int(I.shape[1]/2) + 500]
    # Average to 2DS size bins
    #I_binnned= np.mean(I_subset.reshape((I_subset.shape[0]//10,I_subset.shape[1]//10,-1)),2)
    #small = big.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
    I_binnned = I_subset.reshape([100, 10, 100, 10]).mean(3).mean(1)

    x_bins = np.arange(x[int(I.shape[0]/2)-500], x[int(I.shape[0]/2) + 500], 10)
    
    I_binned_75 = (np.where(I_binnned<0.25, 1, 0))
    I_binned_50 = (np.where(I_binnned<0.5, 1, 0))
    I_binned_25 = (np.where(I_binnned<0.75, 1, 0))
#    plt.pcolormesh(x_bins,x_bins,I_binnned, cmap='jet')
    return x_bins, I_binnned, I_binned_75, I_binned_50, I_binned_25 


#_______________________________________________________________________________________

#Average AveragingFactor elements of 2D array ssuming that each dimension of the new shape is a factor of the corresponding dimension in the old one


def AverageFactorOAPpixels(I, x, y, AveragingFactor, OAP_PixelSize) : 

    I_subset=I[int(I.shape[0]/2)-1000: int(I.shape[0]/2) + 1000, int(I.shape[1]/2)-1000: int(I.shape[1]/2) + 1000]
   
    nsmallx = int(I_subset.shape[0] / AveragingFactor)
    nsmally = int(I_subset.shape[1] / AveragingFactor)
    AveragingFactor =int(AveragingFactor)
    
    # Average to 2DS size bins
    #I_binnned= np.mean(I_subset.reshape((I_subset.shape[0]//10,I_subset.shape[1]//10,-1)),2)
    #I_binnned = big.reshape([nsmall, nbig//nsmall, nsmall, nbig//nsmall]).mean(3).mean(1)
    I_binnned = I_subset.reshape([nsmallx, AveragingFactor, nsmally, AveragingFactor]).mean(3).mean(1)

    x_bins = np.arange(x[int(I.shape[0]/2)-1000], x[int(I.shape[0]/2) + 1000], OAP_PixelSize)
    y_bins = np.arange(x[int(I.shape[1]/2)-1000], x[int(I.shape[1]/2) + 1000], OAP_PixelSize)
    
    I_binned_75 = (np.where(I_binnned<0.25, 1, 0))
    I_binned_50 = (np.where(I_binnned<0.5, 1, 0))
    I_binned_25 = (np.where(I_binnned<0.75, 1, 0))
#    plt.pcolormesh(x_bins,x_bins,I_binnned, cmap='jet')
    return x_bins,y_bins, I_binnned, I_binned_75, I_binned_50, I_binned_25

#_______________________________________________________________________________________
#Calculate stats (diameter, area, etc) of an image

def ImageParticelStats(BinaryImage, OAPPixelSize):
    
    #label_img = label(BinaryImage)
    if np.sum(BinaryImage) > 0 :
        stats = regionprops(BinaryImage, cache=False)
        bbox = [r.bbox for r in stats]
        
        MeanXY = OAPPixelSize * (bbox[0][2] - bbox[0][0] + bbox[0][3] - bbox[0][1]) / 2
        #Area= np.sum(BinaryImage)
        Area= [r.area for r in stats]
        Area=Area[0]
    else: 
        MeanXY =0
        Area = 0
    return MeanXY, Area
    

#_______________________________________________________________________________________
#calculate stats (diameter, area fraction etc) for an image as a funciton of Z


#CP20080401_204059_387_1.jpg
# THE PLOTTING WILL BE WRONG IF NOT STARTED AT z=0


#ShapeVsZ(SavePath,'Rosette_CPI',1,CPIImagePath,1,10)

# ShapeFlag = 0 Create mask 
# ShapeFlag = 1 Use CPI image as mask
# ShapeFlag = 2 Use general binary image as mask   


def ShapeVsZ(SavePath,Prefix,ShapeFlag,SourceImage,PixelSizeInput,PixelSizeOutput):
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    #Prefix='Column20by200'
    #Prefix='Circle125'
    #Prefix= 'Rosette_CPI'
    SaveDataFlag = 1
    PlotDataFlag = 0 
    
    Lambda= 0.658 # CIP laser wavelength
    #Lambda= 0.785 # 2DS laser wavelength
    #Zarray = np.arange(0, 100000, 10000)
    Zarray = np.arange(0, 100000, 2000)
    #Zarray = np.arange(0, 100000, 2000)
    #Zarray = np.arange(0, 50000, 1000)


    # Parameters    
    
    # Set up shape
    
    if ShapeFlag == 0 : 
        
        pixel_size = PixelSizeInput #um
        x_min = -1024*pixel_size # (µm)
        x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        y_min = -1024*pixel_size # (µm)
        y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        x = np.arange(x_min, x_max+.001, pixel_size)
        y = np.arange(y_min, y_max+.001, pixel_size)
        X, Y = np.meshgrid(x, y)
        
        ### Circle
        #D= 150 # um
        #mask = X**2 + Y**2 < (D/2)**2 # Circle
        
        ###Rectangle
        #H=350
        #W=H/3
        #W=30
        #H=120
        W = 50
        H = 300
        mask = (abs(X) < (W/2)) & (abs(Y) < (H/2)) # rectangle
        M= np.ones((x.size,y.size))
        M[mask] = 0
        pixel_size = 1
    
    if ShapeFlag == 1 : # Use CPI image as mask
        pixel_size = PixelSizeInput #um
        x_min = -1024*pixel_size # (µm)
        x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        y_min = -1024*pixel_size # (µm)
        y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        x = np.arange(x_min, x_max+.001, pixel_size)
        y = np.arange(y_min, y_max+.001, pixel_size)
        X, Y = np.meshgrid(x, y)
        #CPIImagePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/Rosette/CP20080401_204059_387_1.jpg'
        M=LoadCpiImageMask(SourceImage)
        pixel_size = 1
    
    if ShapeFlag == 2 : # Use general binary image as mask
        pixel_size = PixelSizeInput #um
        #pixel_size = 0.89 #um
        #pixel_size = 1 #um
        #pixel_size = 0.57 #um
        x_min = -1024*pixel_size # (µm)
        x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        y_min = -1024*pixel_size # (µm)
        y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
        x = np.arange(x_min, x_max+.001, pixel_size)
        y = np.arange(y_min, y_max+.001, pixel_size)
        X, Y = np.meshgrid(x, y)
        M = SourceImage
        
    
    # These areas and diameter are for all pixels where the drop in light intestiy is greater than threshold.
    # Be aware of this when comparing with OShea AMT 2019. 
    DiameterLevel0_BG=np.zeros(len(Zarray))
    DiameterLevel1_BG=np.zeros(len(Zarray))
    DiameterLevel2_BG=np.zeros(len(Zarray))
    Area0_BG =np.zeros(len(Zarray)) 
    Area1_BG =np.zeros(len(Zarray))
    Area2_BG =np.zeros(len(Zarray))
    Area0_Filled =np.zeros(len(Zarray)) 
    Area1_Filled =np.zeros(len(Zarray))
    Area2_Filled =np.zeros(len(Zarray))
    Circularity0_filled =np.zeros(len(Zarray)) 
    Circularity1_filled =np.zeros(len(Zarray))
    Circularity2_filled =np.zeros(len(Zarray))
    BoxDiameterLevel0=np.zeros(len(Zarray))
    BoxDiameterLevel1=np.zeros(len(Zarray))
    BoxDiameterLevel2=np.zeros(len(Zarray))
    DiameterBGy1=np.zeros(len(Zarray))
    DiameterBGx1=np.zeros(len(Zarray))
    Area_BBoxLevel0=np.zeros(len(Zarray))
    Area_BBoxLevel1=np.zeros(len(Zarray))
    Area_BBoxLevel2=np.zeros(len(Zarray))
    Perimeter_BBox_FilledLevel1 =np.zeros(len(Zarray))
    Area_BBox_FilledLevel1 =np.zeros(len(Zarray))
    
    #Z= 10000 # distance from object plane um
    for i in range(len(Zarray)):
        Z= Zarray[i]
        I, A0, fx, fy = compute_diffraction(Z, Lambda, pixel_size, x, y, X, Y, M,0)
        #xOAP, I_binnned, I_binned_75, I_binned_50, I_binned_25= Average2OAPpixels(I, x, y)
        #OAP_PixelSize = 0.57
        #OAP_PixelSize = pixel_size
        AveragingFactor = PixelSizeOutput / pixel_size # This needs to be an integer
        xOAP, yOAP, I_binnned, I_binned_75, I_binned_50, I_binned_25= AverageFactorOAPpixels(I, x, y, AveragingFactor, PixelSizeOutput)
        DiameterLevel0_BG[i], Area0_BG[i], tmp, tmp, Circularity0_filled[i],BoxDiameterLevel0[i], Area_BBoxLevel0[i], Area0_Filled[i],tmp,tmp = ImageParticleMoreStats(I_binned_25, PixelSizeOutput)
        DiameterLevel1_BG[i], Area1_BG[i], DiameterBGx1[i],DiameterBGy1[i], Circularity1_filled[i],BoxDiameterLevel1[i],Area_BBoxLevel1[i], Area1_Filled[i],Perimeter_BBox_FilledLevel1[i], Area_BBox_FilledLevel1[i] = ImageParticleMoreStats(I_binned_50, PixelSizeOutput)
        DiameterLevel2_BG[i], Area2_BG[i], tmp, tmp, Circularity2_filled[i],BoxDiameterLevel2[i],Area_BBoxLevel2[i], Area2_Filled[i],tmp,tmp = ImageParticleMoreStats(I_binned_75, PixelSizeOutput)
        D0 = DiameterLevel1_BG[0]
        if PlotDataFlag == 1 :
            if i % 2 == 0 : 
                Figurename=SavePath+Prefix+'at'+str(Z)+'.png.'
                Zd_true = (4 * Lambda * Z) / ((D0)**2) 
                plot_diffraction(xOAP, yOAP, I_binned_25, I_binned_50,I_binned_75,x,y, M,I,A0, fx, fy, Z,Zd_true, 1, Figurename)

    
        #x_bins, I_binnned, I_binned_75, I_binned_50, I_binned_25= Average2OAPpixels(I, x,y)
    
    AreaFraction0=(Area0_BG-Area1_BG)/Area0_BG # Alow notation from OShea AMT 2019   
    AreaFraction1=(Area1_BG-Area2_BG)/Area0_BG # Amid
    AreaFraction2=Area2_BG/Area0_BG # Ahigh   
    
    if SaveDataFlag == 1 : 
        #FileName= 'Column50by200.h5'
        #FileName= 'Circle_125.h5'
        FileName= Prefix+'.h5'
        try:
            os.remove(SavePath+FileName)
        except OSError:
            pass   
        h5f = h5py.File(SavePath+FileName, 'w')
    
        h5f.create_dataset('Zarray', data=Zarray)
        h5f.create_dataset('DiameterLevel0', data=DiameterLevel0_BG)
        h5f.create_dataset('DiameterLevel1', data=DiameterLevel1_BG)
        h5f.create_dataset('DiameterLevel2', data=DiameterLevel2_BG)
        h5f.create_dataset('AreaFraction0', data=AreaFraction0)
        h5f.create_dataset('AreaFraction1', data=AreaFraction1)
        h5f.create_dataset('AreaFraction2', data=AreaFraction2)
        h5f.create_dataset('Circularity0_filled', data=Circularity0_filled)
        h5f.create_dataset('Circularity1_filled', data=Circularity1_filled)
        h5f.create_dataset('Circularity2_filled', data=Circularity2_filled)
        h5f.create_dataset('BoxDiameterLevel0',data=BoxDiameterLevel0)
        h5f.create_dataset('BoxDiameterLevel1',data=BoxDiameterLevel1)
        h5f.create_dataset('BoxDiameterLevel2',data=BoxDiameterLevel2)
        h5f.create_dataset('Area0_Filled',data=Area0_Filled)
        h5f.create_dataset('Area1_Filled',data=Area1_Filled)
        h5f.create_dataset('Area2_Filled',data=Area2_Filled)
        h5f.create_dataset('Area0_BG',data=Area0_BG)
        h5f.create_dataset('Area1_BG',data=Area1_BG)
        h5f.create_dataset('Area2_BG',data=Area2_BG)
        h5f.create_dataset('DiameterBGx1',data=DiameterBGx1)
        h5f.create_dataset('DiameterBGy1',data=DiameterBGy1)
        h5f.create_dataset('Area_BBoxLevel0',data=Area_BBoxLevel0)
        h5f.create_dataset('Area_BBoxLevel1',data=Area_BBoxLevel1)
        h5f.create_dataset('Area_BBoxLevel2',data=Area_BBoxLevel2)
        h5f.create_dataset('Perimeter_BBox_FilledLevel1',data=Perimeter_BBox_FilledLevel1)
        h5f.create_dataset('Area_BBox_FilledLevel1',data=Area_BBox_FilledLevel1)
      
        
        h5f.close()

    #return Zarray, DiameterLevel0, DiameterLevel1, DiameterLevel2,AreaFraction0, AreaFraction1, AreaFraction2

#_______________________________________________________________________________________

    
def PlotShapeStatsVsZ(SaveFlag):
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
#    FilePreFix='Column30by120'
    FilePreFix='Rosette_CPI'
#    FilePreFix='Circle150'
    #FilePreFix='Column20by200'
    #FilePreFix='Circle125'
    FileName= FilePreFix+'.h5'
#    D0 = 75 #Mean XY at Z=0
#    FileName= 'Column50by200.h5'
#    FileName= 'Circle_125.h5'
    D0 = 190 #Mean XY at Z=0
   
    Lambda= 0.658 # laser wavelength
    Lambda_mm = Lambda /1000
    Data_h5 = h5py.File(SavePath+FileName, 'r')              
    Zarray=np.array(Data_h5['Zarray'])
    DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
    DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
    DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
    AreaFraction0=np.array(Data_h5['AreaFraction0'])
    AreaFraction1=np.array(Data_h5['AreaFraction1'])
    AreaFraction2=np.array(Data_h5['AreaFraction2'])
    Circularity0=np.array(Data_h5['Circularity0'])
    Circularity1=np.array(Data_h5['Circularity1'])
    Circularity2=np.array(Data_h5['Circularity2'])
    BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        
    Zarray_mm= Zarray / 1000
    AreaRatio_1_0 = AreaFraction1 / AreaFraction0
    
    Zd_AMT2019 = GetZd(AreaRatio_1_0,AreaFraction2,AreaFraction1,AreaFraction0,255075)
    Zd8 = Dof_lawson(D0,Lambda,8) /2 # Divide by 2 because only looking at +ve part
    Zd6 = Dof_lawson(D0,Lambda,6) /2 
    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
    
    fig=plt.figure(figsize=(10,5)) 
    plt.rcParams.update({'font.size': 12})
    plt.subplot(1, 2, 1)
#    plt.title('z=0')
    plt.plot(Zarray_mm,DiameterLevel0, color='r', label='25%')
    plt.plot(Zarray_mm,DiameterLevel1, color='k',label='50%')
    plt.plot(Zarray_mm,DiameterLevel2, color='c',label='75%')
    plt.plot(Zarray_mm,BoxDiameterLevel1, color='m',label='50% (BBox)')
    plt.axvline(x=Zd8, color='b', linestyle='--',label='c=8')
    plt.axvline(x=Zd6, color='g', linestyle='--',label='c=6')
    plt.ylabel('Diameter, μm')
    plt.xlabel('Distance from object plane, mm')
    plt.legend()
        
    plt.subplot(1, 2, 2)
#    plt.plot(Zarray_mm,DiameterLevel1, color='k',label='50%')
#    plt.axvline(x=Zd8, color='b', linestyle='--',label='c=8')
#    plt.axvline(x=Zd6, color='g', linestyle='--',label='c=6')
#    plt.ylabel('Diameter, μm')
#    plt.xlabel('Distance from object plane, μm')
#    plt.legend()
    plt.plot(Zd_true, Circularity0, color='r', label='25%')
    plt.plot(Zd_true, Circularity1, color='k',label='50%')
    plt.plot(Zd_true, Circularity2, color='c',label='75%')
    plt.ylabel('Circularity')
    plt.xlabel('Zd')
    plt.xlim([0,9])
    
#    plt.subplot(3, 2, 3)
#    plt.plot(Zd_true, AreaFraction0)
#    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
#    plt.xlabel('Zd')
#    plt.xlim([0,9])

#    plt.subplot(3, 2, 4)
#    plt.plot(Zd_true, AreaFraction1)
#    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
#    plt.xlabel('Zd')
#    plt.xlim([0,9])

#    plt.subplot(3, 2, 5)
#    plt.plot(Zd_true, AreaFraction2)
#    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
#    plt.xlabel('Zd')
#    plt.xlim([0,9])
    
#    plt.subplot(3, 2, 6)
#    plt.plot(Zd_true, Zd_AMT2019,'o')
#    plt.ylabel('Zd SOS AMT 2019')
#    plt.xlabel('Zd true')
#    plt.ylim([3.5,8.5])
#    plt.xlim([3.5,8.5])
    
    if SaveFlag == 1 :
        #SaveName= 'GreyscaleCircle125.png'
        #SaveName= 'GreyscaleColumn50by200.png'
        SaveName= FilePreFix+'circularity2.png'
        plt.savefig(SavePath+SaveName,dpi=200)
        plt.close(fig)
#_______________________________________________________________________________________
  
def LoadCpiImageMask(ImagePath):
    
    MaskArray = np.ones((2048,2048)) 
    
    #ImagePath = 'C:/Users/Admin TEMP/Documents/Diffraction/CPI_images/Rosette/CP20080401_204059_387_1.jpg'
    import_image = (ndimage.imread(ImagePath,flatten=True).astype(float))
    val = filters.threshold_otsu(import_image)
    BinaryImage= (np.where(import_image<val, 0, 1))
    shape=np.shape(BinaryImage)
    xlen=shape[0]
    ylen=shape[1]
    COMx = 100
    COMy = 100
    
    MaskArray[1024-int(xlen/2) : 1024+int(xlen/2+0.5), 1024-int(ylen/2+0.5) :1024+int(ylen/2) ]= BinaryImage
    
                        
    return MaskArray
                        
#_______________________________________________________________________________________

#Get stats for biggest paricle in region
# used example from https://scikit-image.org/docs/dev/auto_examples/segmentation/plot_label.html#sphx-glr-auto-examples-segmentation-plot-label-py


def ImageParticleMoreStats(BinaryImage, OAPPixelSize):
    
    if np.sum(BinaryImage) > 0 :
        #all particles in image
        BoxStats = regionprops(BinaryImage, cache=False)
        Boxbbox = [r.bbox for r in BoxStats]
        MeanXY = OAPPixelSize * (Boxbbox[0][2] - Boxbbox[0][0] + Boxbbox[0][3] - Boxbbox[0][1]) / 2
        Area_BBox = BoxStats[0].area
        #fill internal voids
        FilledImage = ndimage.morphology.binary_fill_holes(BinaryImage).astype(int)
        FilledStats = regionprops(FilledImage, cache=False)
        Perimeter_BBox_Filled = FilledStats[0].perimeter
        Area_BBox_Filled = FilledStats[0].area
        
        #select largest particle in images
        labels_max= SelectLargestParticle(BinaryImage)
        stats = regionprops(labels_max, cache=False)
        bbox = [r.bbox for r in stats]                     
        MeanXY_BG = OAPPixelSize * (bbox[0][2] - bbox[0][0] + bbox[0][3] - bbox[0][1]) / 2
        DiameterBGx = OAPPixelSize * (bbox[0][3] - bbox[0][1])
        DiameterBGy = OAPPixelSize * (bbox[0][2] - bbox[0][0])
        #Area_BG= [r.area for r in stats]
        #Area_BG=Area_BG[0] # pixels     
        Area_BG = stats[0].area
        
        #fill internal voids largest particle in image
        FilledImageBG = ndimage.morphology.binary_fill_holes(labels_max).astype(int) # fill internal voids for circularity 
        FilledStatsBG = regionprops(FilledImageBG, cache=False) 
        Circularity_Filled = FilledStatsBG[0].perimeter** 2 / (4 * np.pi * FilledStatsBG[0].area )
        #Area_Filled = [r.area for r in stats]
        Area_Filled = FilledStatsBG[0].area
    else:
        
        MeanXY_BG =0
        Area_BG =0
        DiameterBGx =0
        DiameterBGy =0
        Circularity_Filled =np.nan
        MeanXY =0
        Area_Filled = 0
        Area_BBox = 0
        Perimeter_BBox_Filled =0
        Area_BBox_Filled=0
        
    return MeanXY_BG, Area_BG, DiameterBGx,DiameterBGy, Circularity_Filled, MeanXY,Area_BBox, Area_Filled,Perimeter_BBox_Filled, Area_BBox_Filled
#
#
#for region in regionprops(label_image):
#    # take regions with large enough areas
#    if region.area >= 100:
#        # draw rectangle around segmented coins
#        minr, minc, maxr, maxc = region.bbox
#        rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
#                                  fill=False, edgecolor='red', linewidth=2)
#        ax.add_patch(rect)
#
#ax.set_axis_off()
#plt.tight_layout()
#plt.show()

#_______________________________________________________________________________________


def SelectLargestParticle(segmentation):
    labels = label(segmentation)
    unique, counts = np.unique(labels, return_counts=True)
    list_seg=list(zip(unique, counts))[1:] # the 0 label is by default background so take the rest
    largest=max(list_seg, key=lambda x:x[1])[0]
    labels_max=(labels == largest).astype(int)
    return labels_max

#_______________________________________________________________________________________

#plot stats vs Zd for several particles in list

def BatchPlotGrey():

    SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    
    #H5FileList = ['Column20by200','Column30by200','Column40by200','Column50by200','Column60by200','Column70by200','Column80by200','Column90by200','Column100by200','Column110by200']
    #D0List = [110, 115, 120, 125, 130,135,140,145,150,155]
    
    #H5FileList = ['Circle125','Rosette_CPI','Column20by200','Column30by200','Column40by200','Column50by200','Column60by200']
    #D0List = [125,190, 110, 115, 120, 125, 130]
    
    H5FileList = ['ColumnAspect3Dxy67','ColumnAspect3Dxy100','ColumnAspect3Dxy133','ColumnAspect3Dxy166','ColumnAspect3Dxy200','ColumnAspect3Dxy233']
    D0List = [67,100, 133, 166, 200, 233]
    
    fig=plt.figure(figsize=(10,10)) 
    plt.rcParams.update({'font.size': 12})
    panel1= plt.subplot(2, 2, 1)    
    plt.ylabel('D / D$_{0}$')
    plt.xlabel('Zd')
    plt.xlim([0,9])
    
    
    panel2=plt.subplot(2, 2, 2)
    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,9])

    panel3=plt.subplot(2, 2, 3)
    
    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,9])

    panel4=plt.subplot(2, 2, 4)
    
    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,9])
    
    
    for i in range(len(H5FileList)) :
        Filename = H5FileList[i]
        Data_h5 = h5py.File(SavePath+Filename+'.h5', 'r')              
        Zarray_mm = np.array(Data_h5['Zarray']) / 1000
        DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
        DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
        DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
        AreaFraction0=np.array(Data_h5['AreaFraction0'])
        AreaFraction1=np.array(Data_h5['AreaFraction1'])
        AreaFraction2=np.array(Data_h5['AreaFraction2'])
        #Circularity0=np.array(Data_h5['Circularity0'])
        #Circularity1=np.array(Data_h5['Circularity1'])
        #Circularity2=np.array(Data_h5['Circularity2'])
        Data_h5.close()
        D0 = D0List[i]
        #D0 = DiameterLevel1[0]
        Lambda= 0.658 # laser wavelength
        Lambda_mm = Lambda /1000
        Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
        
        D_D0 = DiameterLevel1 / D0
        
        panel1.plot(Zd_true, D_D0,label=Filename)
        panel2.plot(Zd_true, AreaFraction0,label=Filename)
        panel3.plot(Zd_true, AreaFraction1,label=Filename)
        panel4.plot(Zd_true, AreaFraction2,label=Filename)

        plt.legend()

    SaveName= SavePath+'ColumnAspect3GreyRatios.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig)




#_______________________________________________________________________________________

#plot stats vs Zd for several particles in sub folders

def BatchCPI_ShapeVsZ():

    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'    
    for filena in os.listdir(CPIPath):
        #if filena.endswith(".h5"):
        print(filena)
        for ImageName in os.listdir(CPIPath+filena):
                if ImageName.endswith(".jpg"):
                    print(ImageName)    
                    Prefix=ImageName.replace('.jpg', '')
                    CPIImagePath=CPIPath+filena+'/'+ImageName
                    ShapeVsZ(CPIPath+filena+'/',Prefix,1,CPIImagePath,1,10)


#ShapeVsZ(SavePath,Prefix,ShapeFlag,SourceImage,PixelSizeInput,PixelSizeOutput)

#_______________________________________________________________________________________

#Plot the CPI greyscale ratios vs z
                    


def PlotALLCPI_ShapeVsZ_bg():

    
    fig=plt.figure(figsize=(10,10)) 
    plt.rcParams.update({'font.size': 14})
    panel1= plt.subplot(2, 2, 1)    
    plt.ylabel('D / D$_{0}$')
    #plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,2.5])
    
    panel2=plt.subplot(2, 2, 2)
    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
    #plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])

    panel3=plt.subplot(2, 2, 3)
    
    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    panel4=plt.subplot(2, 2, 4)
    
    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #CPIPath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/'
    #SavePath=C:\Users\Admin TEMP\Dropbox (The University of Manchester)\Diffraction\CPI_output
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    for filena in os.listdir(CPIPath):
        #if filena.endswith(".h5"):
        print(filena)
        j=0
        LW = 0.1
        for ImageName in os.listdir(CPIPath+filena):
                if ImageName.endswith(".h5"):
                    print(ImageName)    
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
                    DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
                    DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
                    AreaFraction0=np.array(Data_h5['AreaFraction0'])
                    AreaFraction1=np.array(Data_h5['AreaFraction1'])
                    AreaFraction2=np.array(Data_h5['AreaFraction2'])
                    Data_h5.close()
                    #D0 = D0List[i]
                    D0 = DiameterLevel1[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    D_D0 = DiameterLevel1 / D0
        
#                    if j == 0  :
#                        panel1.plot(Zd_true, D_D0,color=ColourList[i],label=filena,linewidth=LW)
#                        panel2.plot(Zd_true, AreaFraction0,color=ColourList[i],label=filena,linewidth=LW)
#                        panel3.plot(Zd_true, AreaFraction1,color=ColourList[i],label=filena,linewidth=LW)
#                        panel4.plot(Zd_true, AreaFraction2,color=ColourList[i],label=filena,linewidth=LW)
#                    else :
                        
                    panel1.plot(Zd_true, D_D0,color='silver',linewidth=LW)
                    panel2.plot(Zd_true, AreaFraction0,color='silver',linewidth=LW)
                    panel3.plot(Zd_true, AreaFraction1,color='silver',linewidth=LW)
                    panel4.plot(Zd_true, AreaFraction2,color='silver',linewidth=LW)               
                
                j+=1

        i+=1
    
    Habit_D0_D_median={}
    HabitList= {'Bullets','ColumnAgg','Columns', 'Plate', 'PlateAgg', 'QuasiSpherical', 'Rosette', 'RosetteAgg'}
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    for filena in HabitList :
        LW=4
        Hist,ZdBinsMid,ZdBins_edge,D_D0Bins_edge, D_D0Bins_mid, D_D0Median,D_D0_IQR, AreaFraction0Median, AreaFraction1Median, AreaFraction2Median = CPI_HabitDmeasD0PDF(CPIPath,filena,0)
        panel1.plot(ZdBinsMid, D_D0Median,linewidth=LW, label=filena)
        panel2.plot(ZdBinsMid, AreaFraction0Median,linewidth=LW)
        panel3.plot(ZdBinsMid, AreaFraction1Median,linewidth=LW)
        panel4.plot(ZdBinsMid, AreaFraction2Median,linewidth=LW, label=filena)
        Habit_D0_D_median[filena+'_median'] = D_D0Median
        Habit_D0_D_median[filena+'_IQR'] = D_D0_IQR
    
    Habit_D0_D_median['ZdBinsMid'] = ZdBinsMid
    panel4.legend()
    
    
    
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPIGreyscaleRatios2_BG.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig)

    return Habit_D0_D_median

#_______________________________________________________________________________________

#plot stats vs Zd for specific habit


def PlotHabitCPI_ShapeVsZ():

    
    fig=plt.figure(figsize=(10,10)) 
    plt.rcParams.update({'font.size': 12})
    panel1= plt.subplot(2, 2, 1)    
    plt.ylabel('D / D$_{0}$')
    plt.xlabel('Zd')
    plt.xlim([0,10])
    plt.ylim([0,2.5])
    
    panel2=plt.subplot(2, 2, 2)
    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,10])
    plt.ylim([0,1])

    panel3=plt.subplot(2, 2, 3)
    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    panel4=plt.subplot(2, 2, 4)
    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
    plt.xlabel('Zd')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    
    
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #CPIPath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/'
    #SavePath=C:\Users\Admin TEMP\Dropbox (The University of Manchester)\Diffraction\CPI_output
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    filena= 'PlateAgg'
    LW = 0.2
    for ImageName in os.listdir(CPIPath+filena):
        if ImageName.endswith(".h5"):
            print(ImageName)    
                    #load cpi h5 file
            Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
            Zarray_mm = np.array(Data_h5['Zarray']) / 1000
            DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
            DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
            DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
            AreaFraction0=np.array(Data_h5['AreaFraction0'])
            AreaFraction1=np.array(Data_h5['AreaFraction1'])
            AreaFraction2=np.array(Data_h5['AreaFraction2'])
            Data_h5.close()
                    #D0 = D0List[i]
            D0 = DiameterLevel1[0]
            Lambda= 0.658 # laser wavelength
            Lambda_mm = Lambda /1000
            Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
            D_D0 = DiameterLevel1 / D0
            panel1.plot(Zd_true, D_D0,color='k',label=filena,linewidth=LW)
            panel2.plot(Zd_true, AreaFraction0,color='k',label=filena,linewidth=LW)
            panel3.plot(Zd_true, AreaFraction1,color='k',label=filena,linewidth=LW)
            panel4.plot(Zd_true, AreaFraction2,color='k',label=filena,linewidth=LW)

        

    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_'+filena+'GreyRatios.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig)

#_______________________________________________________________________________________

# plot zd contours on d vs Z axis


def ZvsDvsZd():
    
    
    Lambda= 0.658 # laser wavelength
    Lambda_mm = Lambda / 1000 
    D_array = np.arange(50, 19200, 50)
    Z_mm_atZd8 =   8 * ((D_array/1000)**2) /  (4 * Lambda_mm ) 
    Z_mm_atZd7 =   7 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd6 =   6 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd5 =   5 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd4 =   4 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd3 =   3 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd2 =   2 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd1 =   1 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd05 =   0.5 * ((D_array/1000)**2) /  (4 * Lambda_mm )
    Z_mm_atZd01 =   0.1 * ((D_array/1000)**2) /  (4 * Lambda_mm )

    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    plt.plot(D_array,Z_mm_atZd8, label = 'Zd = 8')
    plt.plot(D_array,Z_mm_atZd7, label = 'Zd = 7')
    plt.plot(D_array,Z_mm_atZd6, label = 'Zd = 6')
    plt.plot(D_array,Z_mm_atZd5, label = 'Zd = 5')
    plt.plot(D_array,Z_mm_atZd4, label = 'Zd = 4')
    plt.plot(D_array,Z_mm_atZd3, label = 'Zd = 3')
    plt.plot(D_array,Z_mm_atZd2, label = 'Zd = 2')
    plt.plot(D_array,Z_mm_atZd1, label = 'Zd = 1')
    plt.plot(D_array,Z_mm_atZd05, label = 'Zd = 0.5')
    plt.plot(D_array,Z_mm_atZd01, label = 'Zd = 0.1')
    
    plt.legend()    
    plt.axhline(y=35, color='black', linestyle='--')
    plt.axhline(y=81, color='black', linestyle='--')
    #plt.xlim([0,])
    plt.xlim([50,10000])
    plt.ylim([1,100000])
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Diameter, μm')
    plt.ylabel('Distance from object plane, mm')
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/ZvsDvsZd.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig)
    
#_______________________________________________________________________________________

#D/D0 pdf for a habit and also box plots of D/Do vs Zd. Returns histogram 

def CPI_HabitDmeasD0PDF(CPIPath,filena,SavePlot):

    #CPIPath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/'
    #SavePath=C:\Users\Admin TEMP\Dropbox (The University of Manchester)\Diffraction\CPI_output

    
    Batch_Zd=[]
    Batch_D_D0=[]
    Batch_AreaFraction0=[]
    Batch_AreaFraction1=[]
    Batch_AreaFraction2=[]
    
    ZdBins_min = 0
    ZdBins_max = 10
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
           
    for ImageName in os.listdir(CPIPath+filena):
        if ImageName.endswith(".h5"):
            print(ImageName)    
            #load cpi h5 file
            Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
            Zarray_mm = np.array(Data_h5['Zarray']) / 1000
            DiameterLevel0=np.array(Data_h5['DiameterLevel0'])
            DiameterLevel1=np.array(Data_h5['DiameterLevel1'])
            DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
            AreaFraction0=np.array(Data_h5['AreaFraction0'])
            AreaFraction1=np.array(Data_h5['AreaFraction1'])
            AreaFraction2=np.array(Data_h5['AreaFraction2'])
            Data_h5.close()
            D0 = DiameterLevel1[0]
            Lambda= 0.658 # laser wavelength
            Lambda_mm = Lambda /1000
            
            Zd = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
            D_D0 = DiameterLevel1 / D0
            D_D0[np.isnan(D_D0)]=0       
            AreaFraction0[np.isnan(AreaFraction0)] = 0
            AreaFraction1[np.isnan(AreaFraction1)] = 0
            AreaFraction2[np.isnan(AreaFraction2)] = 0
            
            D_D0_interp = np.interp(ZdBins_mid,Zd,D_D0) # interpolate diameter ratio to common zd base
            AreaFraction0_interp = np.interp(ZdBins_mid,Zd,AreaFraction0)
            AreaFraction1_interp = np.interp(ZdBins_mid,Zd,AreaFraction1)
            AreaFraction2_interp = np.interp(ZdBins_mid,Zd,AreaFraction2)
            
            Batch_Zd= np.append(Batch_Zd,ZdBins_mid)
            Batch_D_D0= np.append(Batch_D_D0,D_D0_interp)
            Batch_AreaFraction0= np.append(Batch_AreaFraction0,AreaFraction0_interp)
            Batch_AreaFraction1= np.append(Batch_AreaFraction1,AreaFraction1_interp)
            Batch_AreaFraction2= np.append(Batch_AreaFraction2,AreaFraction2_interp)
    #Bin DataY based on DataX value. Now returns lower, upper and mid bins 

    D_D0Bins_min =-0.05
    D_D0Bins_max =2.5
    D_D0Bins_delta =0.1
    D_D0Bins_edge=np.arange(D_D0Bins_min, D_D0Bins_max, D_D0Bins_delta)
    Hist,xbins_edge,ybins_edge = np.histogram2d(Batch_Zd,Batch_D_D0 , bins=[ZdBins_edge,D_D0Bins_edge]) # Hist y=D/D0, x=Zd and Z= counts
    #normalise Zd    
    RowSum=np.sum(Hist,axis=1)
    PDF=Hist/RowSum[:,np.newaxis] # The sum of all D/D0 at a given Zd will be 1. Row = Zd, Column = D/D0
    #median for bin
    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    D_D0Median = np.nanmedian(Output, axis=0) # median D/D0 at a given Zd
    D_D0_IQR = np.nanpercentile(Output,75, axis=0) - np.nanpercentile(Output,25, axis=0)
    
    D_D0Bins_mid = (D_D0Bins_edge[:-1:1] + D_D0Bins_edge[1::1]) /2
    TotalPDF = np.sum(Hist/np.sum(Hist),axis=0) # D/D0 frequnecy across whole DOF
    
    Output_AreaFraction0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    AreaFraction0Median = np.nanmedian(Output_AreaFraction0, axis=0) # median at a given Zd
    
    Output_AreaFraction1, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction1,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    AreaFraction1Median = np.nanmedian(Output_AreaFraction1, axis=0) # median at a given Zd  
    
    Output_AreaFraction2, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction2,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    AreaFraction2Median = np.nanmedian(Output_AreaFraction2, axis=0) # median at a given Zd

    
    if SavePlot == 1: 
        fig=plt.figure(figsize=(8,8))
        plt.rcParams.update({'font.size': 12})
        plt.subplot(2, 1, 1)
        #plt.pcolormesh(ZdBins_edge,DBins_edge,PDF.transpose(), cmap='viridis') #The grid orientation follows the standard matrix convention: An array C with shape (nrows, ncolumns) is plotted with the column number as X and the row number as Y.
        mask = ~np.isnan(Output) # remove nans
        filtered_data = [d[m] for d, m in zip(Output.T, mask.T)]
        #means = [np.nanmean(x) for x in filtered_data]
        plt.boxplot(filtered_data,whis=[10,90],showfliers=False,positions=ZdBinsMid,widths=(ZdBins_delta*0.75))
        #plt.plot(ZdBinsMid, means,"o", color='black')
        plt.xlabel('Z$_{d}$')
        plt.ylabel('D / D$_{0}$')
    
        majorLocator = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(0.5)
        plt.gca().xaxis.set_major_locator(majorLocator)
        plt.gca().xaxis.set_major_formatter(majorFormatter)
        plt.gca().xaxis.set_minor_locator(minorLocator)# 
        plt.xlim([0,8])
        plt.ylim([0,2])
            
        plt.subplot(2, 1, 2)
        plt.plot(D_D0Bins_mid ,TotalPDF,'o' )
        plt.xlabel('D / D$_{0}$')
        plt.ylabel('Normalised Frequency')
        plt.axhline(y=0, color='r', linestyle='--')
        plt.axvline(x=1, color='r', linestyle='--')
        SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/D_D0_'+filena
        plt.savefig(SaveName,dpi=200)    
        plt.close(fig)


    return PDF,ZdBinsMid,ZdBins_edge,D_D0Bins_edge, D_D0Bins_mid, D_D0Median,D_D0_IQR, AreaFraction0Median, AreaFraction1Median, AreaFraction2Median  


#_______________________________________________________________________________________

# plot median D/D0 for each habit as a function of Zd.

def Compare_CPI_HabitSizePDF():

    HabitList= {'Bullets','ColumnAgg','Columns', 'Plate', 'PlateAgg', 'QuasiSpherical', 'Rosette', 'RosetteAgg'}
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    fig=plt.figure(figsize=(8,8))
    for filena in HabitList :
        Hist,ZdBinsMid,ZdBins_edge,D_D0Bins_edge, D_D0Bins_mid,D_D0Median,D_D0_IQR, AreaFraction0Median, AreaFraction1Median, AreaFraction2Median = CPI_HabitDmeasD0PDF(CPIPath,filena,0)
        plt.plot(ZdBinsMid, D_D0Median, label=filena)
        plt.xlabel('Z$_{d}$')
        plt.ylabel('D / D$_{0}$')

    plt.legend() 
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/D_D0_AllHabits'
    plt.savefig(SaveName,dpi=200)    
    plt.close(fig)


#_______________________________________________________________________________________

# D/D0 pdf vs Do for each habit

def DmeasD0VsD0_habit(filena,SaveFlag) :
#    filena='Rosette'
    Armwidth_mm = 70 
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    PDF,ZdBinsMid,ZdBins_edge,D_D0Bins_edge, D_D0Bins_mid, D_D0Median,D_D0_IQR,  AreaFraction0Median, AreaFraction1Median, AreaFraction2Median= CPI_HabitDmeasD0PDF(CPIPath,filena,0)
    #PDF Row = Zd, Column = D/D0
    D0array=np.arange(10, 1280, 10)
    D0arrayEdge= np.arange(5,1285, 10)


    Lambda= 0.658 # laser wavelength
    Lambda_mm = Lambda /1000        
    Zd_Max = (4 * Lambda_mm * (Armwidth_mm/2) ) / ((D0array/1000)**2) # the Zd of the armwidth as function of D
    Zd_Max_idx= np.zeros(len(Zd_Max))
    
    DmeasD0VsD0= np.zeros([len(D0array),len(D_D0Bins_mid)])
    
    for i in range(len(Zd_Max)):
        if Zd_Max[i] < np.max(ZdBinsMid) :
            Zd_Max_idx[i]=bisect.bisect_left(ZdBinsMid,Zd_Max[i])  
            tmp =np.sum(PDF[:int(Zd_Max_idx[i]),:],axis=0) 
            tmp[0]=np.nan # remove D/D0 = 0
            tmp /= np.nansum(tmp)
            DmeasD0VsD0[i,:]= np.transpose(tmp)
        else : 
            Zd_Max_idx[i] = len(Zd_Max)-1
            tmp = np.sum(PDF[:int(Zd_Max_idx[i]),:],axis=0) 
            tmp[0]=np.nan # remove D/D0 = 0
            tmp /= np.nansum(tmp)
            DmeasD0VsD0[i,:]= np.transpose(tmp)
            
    if SaveFlag == 1: 
        fig=plt.figure(figsize=(10,8))
        plt.rcParams.update({'font.size': 14})
        plt.pcolormesh(D0arrayEdge,D_D0Bins_edge,DmeasD0VsD0.transpose(), cmap='gist_earth') #The grid orientation follows the standard matrix convention: An array C with shape (nrows, ncolumns) is plotted with the column number as X and the row number as Y.
        plt.xscale('log')
        plt.xlabel('D$_{0}$, μm')
        plt.ylabel('D measured / D$_{0}$')
        cbar=plt.colorbar(orientation='vertical')
        cbar.set_label('Normalised frequency', rotation=270, labelpad=20)
        SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/DmeasD0VsD0_'+filena
        plt.savefig(SaveName,dpi=200)    
        plt.close(fig)
    
    DmeasVsD0_pdf, DmeasMidBins, DmeasEdgeBins = DmeasVsD0_habit(DmeasD0VsD0, D0arrayEdge, D0array,D_D0Bins_edge,D_D0Bins_mid,filena)
    
    return DmeasD0VsD0, D0arrayEdge, D0array,D_D0Bins_edge,D_D0Bins_mid


##_______________________________________________________________________________________
#  
    
def DmeasVsD0_habit(DmeasD0VsD0, D0arrayEdge, D0array,D_D0Bins_edge,D_D0Bins_mid,filena):    
#    
#
    DmeasMidBins=np.arange(10, 1280, 10) # new bins
    DmeasEdgeBins=np.arange(5,1285, 10)
    
#    DmeasMidBins=np.arange(10, 1280, 20) # new bins
#    DmeasEdgeBins=np.arange(0,1290, 20)
    
    DmeasVsD0_pdf=np.zeros([len(D0array),len(DmeasMidBins)])
    
    
    for i in range(len(D0array)) :
        DmeasVsD0_i= DmeasD0VsD0[i,:] # this is just frequncy need to change the bins not this
        D_Bins_mid_i = D_D0Bins_mid * D0array[i]
        DmeasVsD0_interp = np.interp(DmeasMidBins,D_Bins_mid_i,DmeasVsD0_i)
        DmeasVsD0_interp /= np.nansum(DmeasVsD0_interp)
        DmeasVsD0_interp[np.isnan(DmeasVsD0_interp)] = 0
        DmeasVsD0_pdf[i,:]= np.transpose(DmeasVsD0_interp)
    
    fig=plt.figure(figsize=(10,8))
    plt.rcParams.update({'font.size': 14})
    plt.pcolormesh(D0arrayEdge,DmeasEdgeBins,DmeasVsD0_pdf.transpose(), cmap='gist_earth',vmin=0,vmax=1) #The grid orientation follows the standard matrix convention: An array C with shape (nrows, ncolumns) is plotted with the column number as X and the row number as Y.
    plt.plot(D0arrayEdge,D0arrayEdge,'--',color='r') # 1 to 1 line
    plt.xlim([5,1285])
    plt.ylim([5,1285])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('D$_{0}$, μm')
    plt.ylabel('D measured, μm')
    cbar=plt.colorbar(orientation='vertical')
    cbar.set_label('Normalised frequency', rotation=270, labelpad=20)
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/DmeasVsD0_'+filena
    plt.savefig(SaveName,dpi=200)    
    plt.close(fig)
    
    
    return  DmeasVsD0_pdf, DmeasMidBins, DmeasEdgeBins


##_______________________________________________________________________________________
#  

#SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
#Prefix= 'Rosette_CPI'
#CPIImagePath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/Rosette/CP20080401_204059_387_1.jpg'

def ShapeVsZ_movie(SavePath,Prefix,CPI_flag,CPIImagePath):
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    #Prefix='Column20by200'
    #Prefix='Circle125'
    #Prefix= 'Rosette_CPI'

    pixel_size = 1 #um
    Lambda= 0.658 # laser wavelength
    #Zarray = np.arange(0, 100000, 2000)
    Zarray = np.arange(0, 100000, 500)
    #Zarray = np.arange(0, 50000, 1000)

    # Parameters    
    x_min = -1024*pixel_size # (µm)
    x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
    y_min = -1024*pixel_size # (µm)
    y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
    #x_min = -512*pixel_size # (µm)
    #x_max = 511*pixel_size # (µm) use 2**n for fast FFT
    #y_min = -512*pixel_size # (µm)
    #y_max = 511*pixel_size # (µm) use 2**n for fast FFT
    x = np.arange(x_min, x_max+.001, pixel_size)
    y = np.arange(y_min, y_max+.001, pixel_size)
    X, Y = np.meshgrid(x, y)

    # Set up shape
    
    if CPI_flag == 1 :
        # Use CPI image as mask 
        #CPIImagePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/Rosette/CP20080401_204059_387_1.jpg'
        M=LoadCpiImageMask(CPIImagePath)
    else :
        ### Circle
        #D= 150 # um
        #mask = X**2 + Y**2 < (D/2)**2 # Circle
        
        ###Rectangle
        #H=350
        #W=H/3
        #W=30
        #H=120
        W = 100
        H = 400 
        mask = (abs(X) < (W/2)) & (abs(Y) < (H/2)) # rectangle
        
        M= np.ones((x.size,y.size))
        M[mask] = 0
        D0 = 250
    
    Ims= np.zeros((100,100,len(Zarray)))
    Zd_array = np.zeros(len(Zarray))
    
    # These areas and diameter are for all pixels where the drop in light intestiy is greater than threshold.
    # Be aware of this when comparing with OShea AMT 2019. 
#    DiameterLevel0=np.zeros(len(Zarray))
#    DiameterLevel1=np.zeros(len(Zarray))
#    DiameterLevel2=np.zeros(len(Zarray))
#    Area0 =np.zeros(len(Zarray)) 
#    Area1 =np.zeros(len(Zarray))
#    Area2 =np.zeros(len(Zarray))
#    Circularity0 =np.zeros(len(Zarray)) 
#    Circularity1 =np.zeros(len(Zarray))
#    Circularity2 =np.zeros(len(Zarray))
#    BoxDiameterLevel0=np.zeros(len(Zarray))
#    BoxDiameterLevel1=np.zeros(len(Zarray))
#    BoxDiameterLevel2=np.zeros(len(Zarray))
    
    #Z= 10000 # distance from object plane um
    for i in range(len(Zarray)):
        Z= Zarray[i]
        I, A0, fx, fy = compute_diffraction(Z, Lambda, pixel_size, x, y, X, Y, M,0)
        xOAP, I_binnned, I_binned_75, I_binned_50, I_binned_25= Average2OAPpixels(I, x, y)
        #DiameterLevel0[i], Area0[i], Circularity0[i],BoxDiameterLevel0[i]= ImageParticleMoreStats(I_binned_25, 10)
        #DiameterLevel1[i], Area1[i], Circularity1[i],BoxDiameterLevel1[i]= ImageParticleMoreStats(I_binned_50, 10)
        #DiameterLevel2[i], Area2[i], Circularity2[i],BoxDiameterLevel2[i]= ImageParticleMoreStats(I_binned_75, 10)
         
            #Figurename=SavePath+Prefix+'at'+str(Z)+'.png.'
        #D0 = DiameterLevel1[0]
        
        Zd_array[i] = (4 * Lambda * Z) / ((D0)**2) 
            #plot_diffraction(xOAP, xOAP, I_binned_25, I_binned_50,I_binned_75,x,y, M,I,A0, fx, fy, Z,Zd_true, 0, Figurename)
            
        imageGreyscale = (I_binned_25 - I_binned_50 )*25 + (I_binned_50 - I_binned_75 )*50 + (I_binned_75 )*75
    
        Ims[:,:,i] = imageGreyscale
       
    ShapeVsZ_movie2(Ims, xOAP, Zarray, Zd_array,SavePath,Prefix)
    
    #return Ims, xOAP, Zarray, Zd_array
    
##_______________________________________________________________________________________
#
    
def ShapeVsZ_movie2(Ims, xOAP, Zarray, Zd_array,SavePath,Prefix): 

    
    
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.rcParams.update({'font.size': 12})
    
    #imageDiffractionThreshold= (np.where(imageDiffraction>0.5, 1, 0))
    
    #plt.title('50% Threshold')
    cax= ax.pcolormesh(xOAP,xOAP,Ims[:-1, :-1, 0], vmin=0, vmax=75, cmap='gist_stern_r')
    plt.ylabel('y, μm')
    plt.xlabel('x, μm')
    #plt.colorbar(orientation='vertical')
    plt.ylim([-300,300])
    plt.xlim([-300,300])
    
    
    
    
    #cax = ax.pcolormesh(x, y, G[:-1, :-1, 0], vmin=-1, vmax=1, cmap='Blues')
    #fig.colorbar(cax)
    
    def animate(i):
        cax.set_array(Ims[:-1, :-1, i].flatten())
        plt.title('z='+str(Zarray[i]) + ' μm, Zd = '+str(np.around(Zd_array[i],2)))

    # Save the animation
    anim = animation.FuncAnimation(fig, animate, frames=(len(Ims[0,0,:])-1), interval=1, repeat=True)
    #fig.show()
    
    # Save as movie 
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=2, bitrate=100000)
    anim.save(SavePath+Prefix+'_movie.mp4', writer=writer)
    
    return anim


##_______________________________________________________________________________________
#

def SaveDmeasVsD0():
    filena='Rosette'
    SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    DmeasD0VsD0, D0arrayEdge, D0array,D_D0Bins_edge,D_D0Bins_mid = DmeasD0VsD0_habit(filena,0)
    DmeasVsD0_pdf, DmeasMidBins, DmeasEdgeBins = DmeasVsD0_habit(DmeasD0VsD0, D0arrayEdge, D0array,D_D0Bins_edge,D_D0Bins_mid,filena)
    
    
            #FileName= 'Column50by200.h5'
            #FileName= 'Circle_125.h5'
    SaveFileName= 'DmeasD0VsD0_'+filena+'.h5'
    try:
        os.remove(SavePath+SaveFileName)
    except OSError:
        pass   
    h5f = h5py.File(SavePath+SaveFileName, 'w')
        
    h5f.create_dataset('DmeasVsD0_pdf',data=DmeasVsD0_pdf)
    h5f.create_dataset('D0arrayEdge',data=D0arrayEdge)
    h5f.create_dataset('D0array',data=D0array)
    h5f.create_dataset('DmeasMidBins',data=DmeasMidBins)
    h5f.create_dataset('DmeasEdgeBins',data=DmeasEdgeBins)
    h5f.close()
    
    
##_______________________________________________________________________________________
#
    
def TestInversion():
    
    #mododisperse distribution
    D0_test = 100 
    Counts = 100
    
    # load probe reponse file
    filena='Rosette'
    SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    SaveFileName= 'DmeasD0VsD0_'+filena+'.h5'
    Data_h5 = h5py.File(SavePath+SaveFileName, 'r')              
    DmeasVsD0_pdf=np.array(Data_h5['DmeasVsD0_pdf'])
    D0arrayEdge=np.array(Data_h5['D0arrayEdge'])
    D0array=np.array(Data_h5['D0array'])
    DmeasMidBins=np.array(Data_h5['DmeasMidBins'])
    DmeasEdgeBins=np.array(Data_h5['DmeasEdgeBins'])
    Data_h5.close()
    
    # probe reponse from monodisperse distribution
    D0_test_idx=bisect.bisect_left(D0array,D0_test) #
    Counts_probe = Counts * DmeasVsD0_pdf[D0_test_idx,:] # 
    

    print(np.nansum(Counts_probe))
    
    #InversionMatrix = np.linalg.inv(DmeasVsD0_pdf)
    #CountsRetrieved = np.matmul(InversionMatrix,Counts_probe)
    #plt.plot(D0array,CountsRetrieved)

    plt.plot(DmeasMidBins,Counts_probe)
    plt.xlabel('Diameter, um')
    plt.ylabel('Counts')


    #return D0array,CountsRetrieved, Counts_probe

###_______________________________________________________________________________________
##  
#    
##_______________________________________________________________________________________
##calculate stats (diameter, area fraction etc) for an image as a funciton of Z
#
#
##CP20080401_204059_387_1.jpg
## THE PLOTTING WILL BE WRONG IF NOT STARTED AT z=0
#
#
##ShapeVsZ(SavePath,'Rosette_CPI',1,CPIImagePath)
#
#
#
#def ShapeVsZ_Analog(BinaryImage,pixel_size):
#
#    SaveDataFlag = 1
#    PlotDataFlag = 0
#    Lambda= 0.658 # laser wavelength
#    #Zarray = np.arange(0, 100000, 2000)
#    Zarray = np.arange(0, 100000, 500)
#    #Zarray = np.arange(0, 50000, 1000)
#
#    # Parameters    
#    x_min = -1024*pixel_size # (µm)
#    x_max = 1023*pixel_size # (µm) use 2**n for fast FFT
#    y_min = -1024*pixel_size # (µm)
#    y_max = 1023*pixel_size # (µm) use 2**n for fast FFT
#    #x_min = -512*pixel_size # (µm)
#    #x_max = 511*pixel_size # (µm) use 2**n for fast FFT
#    #y_min = -512*pixel_size # (µm)
#    #y_max = 511*pixel_size # (µm) use 2**n for fast FFT
#    x = np.arange(x_min, x_max+.001, pixel_size)
#    y = np.arange(y_min, y_max+.001, pixel_size)
#    X, Y = np.meshgrid(x, y)
#
#    # Set up shape
#    
#    if CPI_flag == 1 :
#        # Use CPI image as mask 
#        #CPIImagePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/Rosette/CP20080401_204059_387_1.jpg'
#        M=LoadCpiImageMask(CPIImagePath)
#    else :
#        ### Circle
#        #D= 150 # um
#        #mask = X**2 + Y**2 < (D/2)**2 # Circle
#        
#        ###Rectangle
#        #H=350
#        #W=H/3
#        #W=30
#        #H=120
#        W = 50
#        H = 200 
#        mask = (abs(X) < (W/2)) & (abs(Y) < (H/2)) # rectangle
#        
#        M= np.ones((x.size,y.size))
#        M[mask] = 0
#    
#
#    
#    # These areas and diameter are for all pixels where the drop in light intestiy is greater than threshold.
#    # Be aware of this when comparing with OShea AMT 2019. 
#    DiameterLevel0=np.zeros(len(Zarray))
#    DiameterLevel1=np.zeros(len(Zarray))
#    DiameterLevel2=np.zeros(len(Zarray))
#    Area0 =np.zeros(len(Zarray)) 
#    Area1 =np.zeros(len(Zarray))
#    Area2 =np.zeros(len(Zarray))
#    Circularity0 =np.zeros(len(Zarray)) 
#    Circularity1 =np.zeros(len(Zarray))
#    Circularity2 =np.zeros(len(Zarray))
#    BoxDiameterLevel0=np.zeros(len(Zarray))
#    BoxDiameterLevel1=np.zeros(len(Zarray))
#    BoxDiameterLevel2=np.zeros(len(Zarray))
#    
#    #Z= 10000 # distance from object plane um
#    for i in range(len(Zarray)):
#        Z= Zarray[i]
#        I, A0, fx, fy = compute_diffraction(Z, Lambda, pixel_size, x, y, X, Y, M,0)
#        xOAP, I_binnned, I_binned_75, I_binned_50, I_binned_25= Average2OAPpixels(I, x, y)
#        DiameterLevel0[i], Area0[i], Circularity0[i],BoxDiameterLevel0[i]= ImageParticleMoreStats(I_binned_25, 10)
#        DiameterLevel1[i], Area1[i], Circularity1[i],BoxDiameterLevel1[i]= ImageParticleMoreStats(I_binned_50, 10)
#        DiameterLevel2[i], Area2[i], Circularity2[i],BoxDiameterLevel2[i]= ImageParticleMoreStats(I_binned_75, 10)
#        if PlotDataFlag == 1 : 
#            Figurename=SavePath+Prefix+'at'+str(Z)+'.png.'
#            D0 = DiameterLevel1[0]
#            Zd_true = (4 * Lambda * Z) / ((D0)**2) 
#            plot_diffraction(xOAP, xOAP, I_binned_25, I_binned_50,I_binned_75,x,y, M,I,A0, fx, fy, Z,Zd_true, 1, Figurename)
#
#    
#        #x_bins, I_binnned, I_binned_75, I_binned_50, I_binned_25= Average2OAPpixels(I, x,y)
#    
#    AreaFraction0=(Area0-Area1)/Area0 # Alow notation from OShea AMT 2019   
#    AreaFraction1=(Area1-Area2)/Area0 # Amid
#    AreaFraction2=Area2/Area0 # Ahigh   
#    
#    if SaveDataFlag == 1 : 
#        #FileName= 'Column50by200.h5'
#        #FileName= 'Circle_125.h5'
#        FileName= Prefix+'.h5'
#        try:
#            os.remove(SavePath+FileName)
#        except OSError:
#            pass   
#        h5f = h5py.File(SavePath+FileName, 'w')
#    
#        h5f.create_dataset('Zarray', data=Zarray)
#        h5f.create_dataset('DiameterLevel0', data=DiameterLevel0)
#        h5f.create_dataset('DiameterLevel1', data=DiameterLevel1)
#        h5f.create_dataset('DiameterLevel2', data=DiameterLevel2)
#        h5f.create_dataset('AreaFraction0', data=AreaFraction0)
#        h5f.create_dataset('AreaFraction1', data=AreaFraction1)
#        h5f.create_dataset('AreaFraction2', data=AreaFraction2)
#        h5f.create_dataset('Circularity0', data=Circularity0)
#        h5f.create_dataset('Circularity1', data=Circularity1)
#        h5f.create_dataset('Circularity2', data=Circularity2)
#        h5f.create_dataset('BoxDiameterLevel0',data=BoxDiameterLevel0)
#        h5f.create_dataset('BoxDiameterLevel1',data=BoxDiameterLevel1)
#        h5f.create_dataset('BoxDiameterLevel2',data=BoxDiameterLevel2)
#        h5f.close()
#
#    #return Zarray, DiameterLevel0, DiameterLevel1, DiameterLevel2,AreaFraction0, AreaFraction1, AreaFraction2

#_______________________________________________________________________________________

#plot stats vs Zd for several particles in list

def BatchPlotStatsVsZ():

    SavePath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/'
    
    #H5FileList = ['Column20by200','Column30by200','Column40by200','Column50by200','Column60by200','Column70by200','Column80by200','Column90by200','Column100by200','Column110by200']
    #D0List = [110, 115, 120, 125, 130,135,140,145,150,155]
    
    #H5FileList = ['Circle125','Rosette_CPI','Column20by200','Column30by200','Column40by200','Column50by200','Column60by200']
    #D0List = [125,190, 110, 115, 120, 125, 130]
    
    H5FileList = ['Column50by300','Column50by100']
    D0List = [175, 75]
    
    fig=plt.figure(figsize=(10,10)) 
    plt.rcParams.update({'font.size': 12})
    panel1= plt.subplot(2, 2, 1)    
    plt.ylabel('Diameter, x')
    plt.xlabel('Z')
    #plt.xlim([0,9])
    
    
    panel2=plt.subplot(2, 2, 2)
    plt.ylabel('Diameter, y')
    plt.xlabel('Z')
    #plt.xlim([0,9])

    panel3=plt.subplot(2, 2, 3)
    
    plt.ylabel('Diameter, x')
    plt.xlabel('Zd')
    plt.xlim([0,9])

    panel4=plt.subplot(2, 2, 4)
    
    plt.ylabel('Diameter, y')
    plt.xlabel('Zd')
    plt.xlim([0,9])
    
        
    for i in range(len(H5FileList)) :
        Filename = H5FileList[i]
        Data_h5 = h5py.File(SavePath+Filename+'.h5', 'r')              
        Zarray_mm = np.array(Data_h5['Zarray']) / 1000
        DiameterBGy1=np.array(Data_h5['DiameterBGy1'])
        DiameterBGx1=np.array(Data_h5['DiameterBGx1'])
#        DiameterLevel2=np.array(Data_h5['DiameterLevel2'])
#        AreaFraction0=np.array(Data_h5['AreaFraction0'])
#        AreaFraction1=np.array(Data_h5['AreaFraction1'])
#        AreaFraction2=np.array(Data_h5['AreaFraction2'])
#        #Circularity0=np.array(Data_h5['Circularity0'])
#        #Circularity1=np.array(Data_h5['Circularity1'])
#        #Circularity2=np.array(Data_h5['Circularity2'])
        Data_h5.close()
        D0 = D0List[i]
        #D0 = DiameterLevel1[0]
        Lambda= 0.658 # laser wavelength
        Lambda_mm = Lambda /1000
        Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
        
        #D_D0 = DiameterLevel1 / D0
        
        panel1.plot(Zarray_mm, DiameterBGx1,label=Filename)
        panel2.plot(Zarray_mm, DiameterBGy1,label=Filename)
        panel3.plot(Zd_true, DiameterBGx1,label=Filename)
        panel4.plot(Zd_true, DiameterBGy1,label=Filename)

        #plt.legend()

    #SaveName= SavePath+'ColumnAspect3GreyRatios.png'
    #plt.savefig(SaveName,dpi=200)
    #plt.close(fig)
    

#_______________________________________________________________________________________

#Plot K07 corrected diameter vs Zd for CPI images 
                    
def PlotALLCPI_Korolev():

    
    #fig=plt.figure(figsize=(12,6)) 
    #plt.rcParams.update({'font.size': 12})
    #panel1= plt.subplot(1, 2, 1)    
    panel1= plt.subplot(3,2, 5)
    plt.ylabel('D / D$_{0}$ (K07)' )
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,2.5]) 
    #panel2= plt.subplot(1, 2, 2)
    panel2= plt.subplot(3, 2, 6)
    plt.xlabel('Normalised Frequency')
    plt.xlim([5E-4,1])
    plt.ylim([0,2.5])
    plt.xscale('log')
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    ZdBins_min = 0
    ZdBins_max = 10
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Median = {}
    IQR={}
    Hist={}
    HistBinsEdge=np.arange(0,2.5,0.1)
    HistBinsMid = (HistBinsEdge[:-1] + HistBinsEdge[1:])/2
    
    AllBatch_D_KorolevCorr_D0 = []
    AllBatch_Zd =[]
    for filena in os.listdir(CPIPath): #loop habits
        if filena != 'Droplets'  :
            #if filena.endswith(".h5"):
            print(filena)
            LW = 0.1
            Batch_D_KorolevCorr_D0 = []
            Batch_Zd = []
            for ImageName in os.listdir(CPIPath+filena): #loop ice crystals
                if ImageName.endswith(".h5"): # 
                    print(ImageName)    
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
                    Area_BBox_FilledLevel1=np.array(Data_h5['Area_BBox_FilledLevel1'])
                    Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                    #Dcircle =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                    Data_h5.close()
                    D0 = BoxDiameterLevel1[0]
                    #D0=Dcircle[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    VoidArea = Area_BBox_FilledLevel1- Area_BBoxLevel1
                    #D_KorolevCorr=KorolevCorrectedD(Area_BBox_FilledLevel1, VoidArea,BoxDiameterLevel1)
                    D_KorolevCorr=KorolevCorrectedD(Area_BBox_FilledLevel1, VoidArea,BoxDiameterLevel1)
                    D_KorolevCorr_D0 = D_KorolevCorr / D0
                    D_KorolevCorr_D0_interp = np.interp(ZdBins_mid,Zd_true,D_KorolevCorr_D0) # interpolate diameter ratio to common zd base
                    Batch_D_KorolevCorr_D0= np.append(Batch_D_KorolevCorr_D0,D_KorolevCorr_D0_interp)
                    Batch_Zd = np.append(Batch_Zd,ZdBins_mid)
                    panel1.plot(Zd_true, D_KorolevCorr_D0,color='silver',linewidth=LW)             
                
            #Median and IQR for habit
            Output_D_KorolevCorr_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_KorolevCorr_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
            Median[filena] = np.nanmedian(Output_D_KorolevCorr_D0, axis=0) # median at a given Zd
            IQR[filena] = (np.nanpercentile(Output_D_KorolevCorr_D0,75, axis =0) - np.nanpercentile(Output_D_KorolevCorr_D0,25, axis =0))
            #Habit histogram
            D_D0hist,tmp = np.histogram(Batch_D_KorolevCorr_D0,bins=HistBinsEdge)
            D_D0hist = D_D0hist / np.nansum(D_D0hist)
            Hist[filena] = D_D0hist
                
    AllBatch_D_KorolevCorr_D0 = np.append(AllBatch_D_KorolevCorr_D0,Batch_D_KorolevCorr_D0)
    AllBatch_Zd =np.append(AllBatch_Zd,Batch_Zd)
            
    # Median and IQR for all ice crystals
    Output_D_KorolevCorr_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(AllBatch_D_KorolevCorr_D0,AllBatch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    Median['All'] = np.nanmedian(Output_D_KorolevCorr_D0, axis=0) # median at a given Zd
    IQR['All'] = (np.nanpercentile(Output_D_KorolevCorr_D0,75, axis =0) - np.nanpercentile(Output_D_KorolevCorr_D0,25, axis =0))

    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
            LW=4
            #Hist,ZdBinsMid,ZdBins_edge,D_D0Bins_edge, D_D0Bins_mid, D_D0Median,D_D0_IQR, AreaFraction0Median, AreaFraction1Median, AreaFraction2Median = CPI_HabitDmeasD0PDF(CPIPath,filena,0)
            panel1.plot(ZdBinsMid, Median[filena],linewidth=LW, label=filena)
            #panel1.legend()
            panel2.plot(Hist[filena],HistBinsMid,linewidth=LW, label=filena)
    
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_K07_habits.png'
    #plt.savefig(SaveName,dpi=200)
    #plt.close(fig) 
    #return Hist, HistBinsMid,HistBinsEdge 

#_______________________________________________________________________________________

# Return median and IQR Dof 75% threshold for all ice crystals
                   
def FindDoF_Zd_75():

        
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    ZdBins_min = 0
    ZdBins_max = 20
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Batch_Dof = np.empty(0)
    PixelSize =10
    
    for filena in os.listdir(CPIPath): #loop habits

        for ImageName in os.listdir(CPIPath+filena): # loop crystals
            if ImageName.endswith(".h5"):
                #load cpi h5 file
                Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                Area_BBoxLevel2=np.array(Data_h5['Area_BBoxLevel2'])
                Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                Data_h5.close()
                #D0 = BoxDiameterLevel1[0]
                DcircleLevel1 =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                D0 = DcircleLevel1[0]
                
                Lambda= 0.658 # laser wavelength
                Lambda_mm = Lambda /1000
                Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                
                Area_BBoxLevel2_interp = np.interp(ZdBins_mid,Zd_true,Area_BBoxLevel2) # interpolate diameter ratio to common zd base             
                i, Dof = ReturnXwhenY_0(Area_BBoxLevel2_interp, ZdBins_mid)
                Batch_Dof= np.append(Batch_Dof,Dof)
                
        print(filena)
        #print(len(Batch_Dof))
    print(np.round(np.nanmedian(Batch_Dof),1))
    print(np.round(np.nanpercentile(Batch_Dof,75) - np.nanpercentile(Batch_Dof,25),1))
        
    return Batch_Dof

#_________________________________________________________________________________


#_______________________________________________________________________________________

# Return median and IQR Dof of habit
                   
def FindDoF_Zd():

        
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    ZdBins_min = 0
    ZdBins_max = 20
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    #Batch_Dof = {}
    PixelSize =10
    
    for filena in os.listdir(CPIPath):
        #if filena.endswith(".h5"):
        #Batch_Dof = {}
        Batch_Dof=np.nan
        for ImageName in os.listdir(CPIPath+filena):
            if ImageName.endswith(".h5"):
                #load cpi h5 file
                Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
                Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                Data_h5.close()
                #D0 = BoxDiameterLevel1[0]
                DcircleLevel1 =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                D0 = DcircleLevel1[0]
                
                Lambda= 0.658 # laser wavelength
                Lambda_mm = Lambda /1000
                Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                
                DiameterLevel1_interp = np.interp(ZdBins_mid,Zd_true,DcircleLevel1) # interpolate diameter ratio to common zd base             
                i, Dof = ReturnXwhenY_0(DiameterLevel1_interp, ZdBins_mid)
                Batch_Dof= np.append(Batch_Dof,Dof)
                
        print(filena)
        #print(len(Batch_Dof))
        print(np.round(np.nanmedian(Batch_Dof),1))
        print(np.round(np.nanpercentile(Batch_Dof,75) - np.nanpercentile(Batch_Dof,25),1))
        
        #return Batch_Dof

#_________________________________________________________________________________
#returns x for first y that <=0

def ReturnXwhenY_0(yarray, xarray):

    i= 0
    while (yarray[i] > 0) and (i < len(xarray)-1): 
        i+=1
    
    return i-1, xarray[i-1]
        
            
#_________________________________________________________________________________
    
#output medians at given zd values

def Habit_D0_D_median_output(Median , IQR, ZdBinsMid): 
    
    HabitList= {'Bullets','ColumnAgg','Columns', 'Plate', 'PlateAgg', 'QuasiSpherical', 'Rosette', 'RosetteAgg'}
    Zd_idx = [10,20,30,40,50,60,70,80,90]
    for Habit in HabitList :
        print(Habit)
        for idx in Zd_idx :
            print(str(np.round(Median[Habit][idx],1))+' ('+str(np.round(IQR[Habit][idx],1))+')')
        
 
    
#_________________________________________________________________________________    

# plot diameter with and without K07

def PlotAllCPI_D_D0_6panel():
    fig=plt.figure(figsize=(14,18)) 
    plt.rcParams.update({'font.size': 16})
    PlotAllCPI_D_CircleEqiv()
    PlotAllCPI_D_BBox()
    PlotALLCPI_Korolev()
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_D_D0_habits.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig) 
            

#_______________________________________________________________________________________

#Plot bbox diameter for CPI images 
                    
def PlotAllCPI_D_BBox():

    #fig=plt.figure(figsize=(12,6)) 
    #plt.rcParams.update({'font.size': 12}) 
    #panel1= plt.subplot(1,2, 1)
    panel1= plt.subplot(3,2, 3)
    plt.ylabel('D / D$_{0}$ (Mean X-Y)' )
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,2.5]) 
    #panel2= plt.subplot(1, 2, 2)
    panel2= plt.subplot(3, 2, 4)
    plt.xlabel('Normalised Frequency')
    plt.xlim([5E-4,1])
    plt.ylim([0,2.5])
    plt.xscale('log')
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    ZdBins_min = 0
    ZdBins_max = 10
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Median = {}
    IQR={}
    Hist={}
    HistBinsEdge=np.arange(0,2.5,0.1)
    HistBinsMid = (HistBinsEdge[:-1] + HistBinsEdge[1:])/2
    
    AllBatch_D_D0 = []
    AllBatch_Zd =[]
    for filena in os.listdir(CPIPath): #loop habits
        if filena != 'Droplets'  :
            #if filena.endswith(".h5"):
            print(filena)
            LW = 0.1
            Batch_D_D0 = []
            Batch_Zd = []
            for ImageName in os.listdir(CPIPath+filena): #loop ice crystals
                if ImageName.endswith(".h5"): # 
                    print(ImageName)    
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
                    #Area_BBox_FilledLevel1=np.array(Data_h5['Area_BBox_FilledLevel1'])
                    #Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                    #Dcircle =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                    Data_h5.close()
                    D0 = BoxDiameterLevel1[0]
                    #D0=Dcircle[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    D_D0 = BoxDiameterLevel1 / D0
                    D_D0_interp = np.interp(ZdBins_mid,Zd_true,D_D0) # interpolate diameter ratio to common zd base
                    Batch_D_D0= np.append(Batch_D_D0,D_D0_interp)
                    Batch_Zd = np.append(Batch_Zd,ZdBins_mid)
                    panel1.plot(Zd_true, D_D0,color='silver',linewidth=LW)             
                
            #Median and IQR for habit
            Output_D_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
            Median[filena] = np.nanmedian(Output_D_D0, axis=0) # median at a given Zd
            IQR[filena] = (np.nanpercentile(Output_D_D0,75, axis =0) - np.nanpercentile(Output_D_D0,25, axis =0))
            #Habit histogram
            D_D0hist,tmp = np.histogram(Batch_D_D0,bins=HistBinsEdge)
            D_D0hist = D_D0hist / np.nansum(D_D0hist)
            Hist[filena] = D_D0hist
                
    AllBatch_D_D0 = np.append(AllBatch_D_D0,Batch_D_D0)
    AllBatch_Zd =np.append(AllBatch_Zd,Batch_Zd)
            
    # Median and IQR for all ice crystals
    Output_D_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(AllBatch_D_D0,AllBatch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    Median['All'] = np.nanmedian(Output_D_D0, axis=0) # median at a given Zd
    IQR['All'] = (np.nanpercentile(Output_D_D0,75, axis =0) - np.nanpercentile(Output_D_D0,25, axis =0))

    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
            LW=4
            panel1.plot(ZdBinsMid, Median[filena],linewidth=LW, label=filena)
            #panel1.legend()
            panel2.plot(Hist[filena],HistBinsMid,linewidth=LW, label=filena)
    
    #SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_D_D0_bbox_habits.png'
    #plt.savefig(SaveName,dpi=200)
    #plt.close(fig) 
    
    #return Median , IQR, ZdBinsMid
#_______________________________________________________________________________________


#Plot circle equivalent K07 corrected diameter for CPI images 
                    
def PlotAllCPI_D_CircleEqiv():
    PixelSize = 10
    
    #fig=plt.figure(figsize=(12,6)) 
    #plt.rcParams.update({'font.size': 12}) 
    #panel1= plt.subplot(1,2, 1)
    panel1= plt.subplot(3,2, 1)
    plt.ylabel('D / D$_{0}$ (Circle equivalent)' )
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,2.5]) 
    #panel2= plt.subplot(1, 2, 2)
    panel2= plt.subplot(3, 2, 2)
    plt.xlabel('Normalised Frequency')
    plt.xlim([5E-4,1])
    plt.ylim([0,2.5])
    plt.xscale('log')
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    i=0
    ZdBins_min = 0
    ZdBins_max = 10
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Median = {}
    IQR={}
    Hist={}
    HistBinsEdge=np.arange(0,2.5,0.1)
    HistBinsMid = (HistBinsEdge[:-1] + HistBinsEdge[1:])/2
    
    AllBatch_D_D0 = []
    AllBatch_Zd =[]
    for filena in os.listdir(CPIPath): #loop habits
        if filena != 'Droplets'  :
            #if filena.endswith(".h5"):
            print(filena)
            LW = 0.1
            Batch_D_D0 = []
            Batch_Zd = []
            for ImageName in os.listdir(CPIPath+filena): #loop ice crystals
                if ImageName.endswith(".h5"): # 
                    print(ImageName)    
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    #BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
                    #Area_BBox_FilledLevel1=np.array(Data_h5['Area_BBox_FilledLevel1'])
                    Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                    Dcircle =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                    Data_h5.close()
                    #D0 = BoxDiameterLevel1[0]
                    D0=Dcircle[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    D_D0 = Dcircle / D0
                    D_D0_interp = np.interp(ZdBins_mid,Zd_true,D_D0) # interpolate diameter ratio to common zd base
                    Batch_D_D0= np.append(Batch_D_D0,D_D0_interp)
                    Batch_Zd = np.append(Batch_Zd,ZdBins_mid)
                    panel1.plot(Zd_true, D_D0,color='silver',linewidth=LW)             
                
            #Median and IQR for habit
            Output_D_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
            Median[filena] = np.nanmedian(Output_D_D0, axis=0) # median at a given Zd
            IQR[filena] = (np.nanpercentile(Output_D_D0,75, axis =0) - np.nanpercentile(Output_D_D0,25, axis =0))
            #Habit histogram
            D_D0hist,tmp = np.histogram(Batch_D_D0,bins=HistBinsEdge)
            D_D0hist = D_D0hist / np.nansum(D_D0hist)
            Hist[filena] = D_D0hist
                
    AllBatch_D_D0 = np.append(AllBatch_D_D0,Batch_D_D0)
    AllBatch_Zd =np.append(AllBatch_Zd,Batch_Zd)
            
    # Median and IQR for all ice crystals
    Output_D_D0, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(AllBatch_D_D0,AllBatch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
    Median['All'] = np.nanmedian(Output_D_D0, axis=0) # median at a given Zd
    IQR['All'] = (np.nanpercentile(Output_D_D0,75, axis =0) - np.nanpercentile(Output_D_D0,25, axis =0))

    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
            LW=4
            panel1.plot(ZdBinsMid, Median[filena],linewidth=LW, label=filena)
            panel2.plot(Hist[filena],HistBinsMid,linewidth=LW, label=filena)
    panel2.legend(fontsize=14,ncol=2)
    #SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_D_D0_bbox_habits.png'
    #plt.savefig(SaveName,dpi=200)
    #plt.close(fig) 
    
    #return Output_D_D0, Batch_D_D0, Batch_Zd



#_______________________________________________________________________________________


#plot cpi grayscale ratios using circle equivalent and mean x-y

def PlotAllCPI_Greyscale_6panel():
    fig=plt.figure(figsize=(14,18)) 
    plt.rcParams.update({'font.size': 16})
    PlotAllCPI_ShapeVsZ_CircleEquiv()
    PlotAllCPI_ShapeVsZ_BBox()
    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPI_Greyscale_habits.png'
    plt.savefig(SaveName,dpi=200)
    plt.close(fig) 



#_______________________________________________________________________________________

#Plot the CPI greyscale ratios vs z using circle equivalent D
                    


def PlotAllCPI_ShapeVsZ_CircleEquiv():

    PixelSize = 10
    
    #fig=plt.figure(figsize=(10,10)) 
    #plt.rcParams.update({'font.size': 14})
    #fig.suptitle('Bounding box circle equivalent diameter')
    
#    panel1= plt.subplot(3, 2, 2)    
#    plt.title('Bounding box circle equivalent diameter')
#    plt.ylabel('D / D$_{0}$')
#    #plt.xlabel('Z$_{d}$')
#    plt.xlim([0,10])
#    plt.ylim([0,2.5])
    
    panel2=plt.subplot(3, 2, 2)
    plt.title('Circle equivalent diameter')
    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
    #plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])

    panel3=plt.subplot(3, 2, 4)
    
    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    panel4=plt.subplot(3, 2, 6)
    
    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    ZdBins_min = 0
    ZdBins_max = 20
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Median = {}
    
    
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #CPIPath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/'
    #SavePath=C:\Users\Admin TEMP\Dropbox (The University of Manchester)\Diffraction\CPI_output
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
        #if filena.endswith(".h5"):
            print(filena)
            LW = 0.1
            Batch_D_D0 = []
            Batch_Zd = []
            Batch_AreaFraction0 = []
            Batch_AreaFraction1 = []
            Batch_AreaFraction2 = []
    
            for ImageName in os.listdir(CPIPath+filena):
                if ImageName.endswith(".h5"):      
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    Area_BBoxLevel0=np.array(Data_h5['Area_BBoxLevel0'])
                    Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                    Area_BBoxLevel2=np.array(Data_h5['Area_BBoxLevel2'])
                    Data_h5.close()
                        
                    DcircleLevel1 =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                    D0 = DcircleLevel1[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    D_D0 = DcircleLevel1 / D0
                        
                    AreaFraction0=(Area_BBoxLevel0-Area_BBoxLevel1)/Area_BBoxLevel0 # Alow notation from OShea AMT 2019   
                    AreaFraction1=(Area_BBoxLevel1-Area_BBoxLevel2)/Area_BBoxLevel0 # Amid
                    AreaFraction2=Area_BBoxLevel2/Area_BBoxLevel0 # Ahigh   
                  
                    D_D0_interp = np.interp(ZdBins_mid,Zd_true,D_D0) # interpolate diameter ratio to common zd base
                    AreaFraction0_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction0) # interpolate to common zd base
                    AreaFraction1_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction1) # interpolate to common zd base
                    AreaFraction2_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction2) # interpolate to common zd base
    
                    Batch_Zd = np.append(Batch_Zd,ZdBins_mid)
                    Batch_D_D0= np.append(Batch_D_D0,D_D0_interp)
                    Batch_AreaFraction0 = np.append(Batch_AreaFraction0,AreaFraction0_interp)
                    Batch_AreaFraction1 = np.append(Batch_AreaFraction1,AreaFraction1_interp)
                    Batch_AreaFraction2 = np.append(Batch_AreaFraction2,AreaFraction2_interp)
                                             
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'D_D0'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction0'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction1,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction1'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction2,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction2'] = np.nanmedian(Output, axis =0)
                    
                    #panel1.plot(Zd_true, D_D0,color='silver',linewidth=LW)
                    panel2.plot(Zd_true, AreaFraction0,color='silver',linewidth=LW)
                    panel3.plot(Zd_true, AreaFraction1,color='silver',linewidth=LW)
                    panel4.plot(Zd_true, AreaFraction2,color='silver',linewidth=LW)               
                
    #Habit_D0_D_median={}
    #HabitList= {'Bullets','ColumnAgg','Columns', 'Plate', 'PlateAgg', 'QuasiSpherical', 'Rosette', 'RosetteAgg'}
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
            LW = 4
            #panel1.plot(ZdBinsMid, Median[filena,'D_D0'],linewidth=LW)
            panel2.plot(ZdBinsMid, Median[filena,'AreaFraction0'],linewidth=LW)
            panel3.plot(ZdBinsMid, Median[filena,'AreaFraction1'],linewidth=LW)
            panel4.plot(ZdBinsMid, Median[filena,'AreaFraction2'],linewidth=LW, label=filena)
            
    panel4.legend(fontsize=14,ncol=2)
    
#    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPIGreyscaleRatios2_CircleEquiv.png'
#    plt.savefig(SaveName,dpi=200)
#    plt.close(fig)



    
#_______________________________________________________________________________________

#Plot the CPI greyscale ratios vs z using circle equivalent D
                    


def PlotAllCPI_ShapeVsZ_BBox():

    PixelSize = 10
    
#    fig=plt.figure(figsize=(10,10)) 
#    plt.rcParams.update({'font.size': 14})

#    panel1= plt.subplot(2, 2, 1)    
#    plt.ylabel('D / D$_{0}$')
#    #plt.xlabel('Z$_{d}$')
#    plt.xlim([0,10])
#    plt.ylim([0,2.5])
    
    panel2=plt.subplot(3, 2, 1)
    plt.title('Mean X-Y')
    plt.ylabel('A$_{25-50}$ / A$_{25-100}$')
    #plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])

    panel3=plt.subplot(3, 2, 3)
    
    plt.ylabel('A$_{50-75}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    panel4=plt.subplot(3, 2, 5)
    
    plt.ylabel('A$_{75-100}$ / A$_{25-100}$')
    plt.xlabel('Z$_{d}$')
    plt.xlim([0,10])
    plt.ylim([0,1])
    
    ZdBins_min = 0
    ZdBins_max = 20
    ZdBins_delta = 0.1
    ZdBins_edge = np.arange(ZdBins_min, ZdBins_max, ZdBins_delta)
    ZdBins_mid= (ZdBins_edge[:-1:1] + ZdBins_edge[1::1]) /2
    Median = {}
    
    
    #SavePath = 'C:/Users/Admin TEMP/Documents/Diffraction/Plots/'
    #CPIPath = 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/CPI_images/'
    #SavePath=C:\Users\Admin TEMP\Dropbox (The University of Manchester)\Diffraction\CPI_output
    
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    
    ColourList = ['r', 'b', 'g','c','m','y','tab:orange', 'tab:pink','tab:olive']
    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
        #if filena.endswith(".h5"):
            print(filena)
            LW = 0.1
            Batch_D_D0 = []
            Batch_Zd = []
            Batch_AreaFraction0 = []
            Batch_AreaFraction1 = []
            Batch_AreaFraction2 = []
    
            for ImageName in os.listdir(CPIPath+filena):
                if ImageName.endswith(".h5"):      
                    #load cpi h5 file
                    Data_h5 = h5py.File(CPIPath+filena+'/'+ImageName, 'r')              
                    Zarray_mm = np.array(Data_h5['Zarray']) / 1000
                    BoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
                    Area_BBoxLevel0=np.array(Data_h5['Area_BBoxLevel0'])
                    Area_BBoxLevel1=np.array(Data_h5['Area_BBoxLevel1'])
                    Area_BBoxLevel2=np.array(Data_h5['Area_BBoxLevel2'])
                    Data_h5.close()
                        
                    #BoxDiameterLevel1 =  PixelSize *2* np.sqrt(Area_BBoxLevel1/np.pi)
                    D0 = BoxDiameterLevel1[0]
                    Lambda= 0.658 # laser wavelength
                    Lambda_mm = Lambda /1000
                    Zd_true = (4 * Lambda_mm * Zarray_mm) / ((D0/1000)**2)
                    D_D0 = BoxDiameterLevel1 / D0
                        
                    AreaFraction0=(Area_BBoxLevel0-Area_BBoxLevel1)/Area_BBoxLevel0 # Alow notation from OShea AMT 2019   
                    AreaFraction1=(Area_BBoxLevel1-Area_BBoxLevel2)/Area_BBoxLevel0 # Amid
                    AreaFraction2=Area_BBoxLevel2/Area_BBoxLevel0 # Ahigh   
                  
                    D_D0_interp = np.interp(ZdBins_mid,Zd_true,D_D0) # interpolate diameter ratio to common zd base
                    AreaFraction0_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction0) # interpolate to common zd base
                    AreaFraction1_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction1) # interpolate to common zd base
                    AreaFraction2_interp = np.interp(ZdBins_mid,Zd_true,AreaFraction2) # interpolate to common zd base
    
                    Batch_Zd = np.append(Batch_Zd,ZdBins_mid)
                    Batch_D_D0= np.append(Batch_D_D0,D_D0_interp)
                    Batch_AreaFraction0 = np.append(Batch_AreaFraction0,AreaFraction0_interp)
                    Batch_AreaFraction1 = np.append(Batch_AreaFraction1,AreaFraction1_interp)
                    Batch_AreaFraction2 = np.append(Batch_AreaFraction2,AreaFraction2_interp)
                                             
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_D_D0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'D_D0'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction0,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction0'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction1,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction1'] = np.nanmedian(Output, axis =0)
                    Output, ZdBinsLower, ZdBinsUpper, ZdBinsMid=BinLinearV2(Batch_AreaFraction2,Batch_Zd,ZdBins_min,ZdBins_max-ZdBins_delta,int((ZdBins_max-ZdBins_min-ZdBins_delta)/ ZdBins_delta))
                    Median[filena,'AreaFraction2'] = np.nanmedian(Output, axis =0)
                    
                    #panel1.plot(Zd_true, D_D0,color='silver',linewidth=LW)
                    panel2.plot(Zd_true, AreaFraction0,color='silver',linewidth=LW)
                    panel3.plot(Zd_true, AreaFraction1,color='silver',linewidth=LW)
                    panel4.plot(Zd_true, AreaFraction2,color='silver',linewidth=LW)               
                
    #Habit_D0_D_median={}
    #HabitList= {'Bullets','ColumnAgg','Columns', 'Plate', 'PlateAgg', 'QuasiSpherical', 'Rosette', 'RosetteAgg'}
    CPIPath = 'C:/Users/Admin TEMP/Documents/CPI_Habit_Images/'
    for filena in os.listdir(CPIPath):
        if filena != 'Droplets'  :
            LW = 4
            #panel1.plot(ZdBinsMid, Median[filena,'D_D0'],linewidth=LW)
            panel2.plot(ZdBinsMid, Median[filena,'AreaFraction0'],linewidth=LW)
            panel3.plot(ZdBinsMid, Median[filena,'AreaFraction1'],linewidth=LW)
            panel4.plot(ZdBinsMid, Median[filena,'AreaFraction2'],linewidth=LW, label=filena)
            
    #panel4.legend()
    
#    SaveName= 'C:/Users/Admin TEMP/Dropbox (The University of Manchester)/Diffraction/Plots/CPIGreyscaleRatios2_BBox.png'
#    plt.savefig(SaveName,dpi=200)
#    plt.close(fig)

    #return Habit_D0_D_median


#_______________________________________________________________________________________


    