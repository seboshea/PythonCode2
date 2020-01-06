# -*- coding: utf-8 -*-

import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import math
import bisect
import h5py
import scipy.io as sio
import matplotlib.pyplot as plt  
import matplotlib.dates as mdates
import matplotlib.colors as colors
#from matplotlib.mlab import bivariate_normal
from matplotlib.dates import DateFormatter
import os
import gc
from scipy import ndimage
from skimage import filters

from OAP import PlotImage, LoadImage, PlotImage_GS
from skimage.measure import label, regionprops
from AngularSpectrumTheory import SelectLargestParticle, ShapeVsZ, ImageParticleMoreStats
from GrayscaleFunctions import LoadGrayScaleStats, GrayscaleRatios 
from MyFunctions import BinLinear_median


#import GreyscaleFunctions
#____________________________________________________________________________
# Load stats from 2DS oasis h5 file 


def Load2DS_h5(ImagePath, FileName):
    PixelSize =10 
    Data_h5 = h5py.File(ImagePath+FileName, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    BufferNumber = ParticleTimesWv[:,6]
    Data_h5.close()

    #Get particle buffer time in seconds from midnight. 3 decimal places
    ParticleBufferTimeS= np.zeros(len(BufferNumber))*np.nan
    AnalogFlag = np.zeros(len(BufferNumber)) 
    
    for i in range(len(BufferNumber)):
        if BufferNumber[i] == 4294967295 :
            ParticleBufferTimeS[i] = np.nan
        else:   
            ParticleBufferTimeS[i] = 3600*HeaderMatrixWv[BufferNumber[i],3]+ 60*HeaderMatrixWv[BufferNumber[i],4] + HeaderMatrixWv[BufferNumber[i],5] + (HeaderMatrixWv[BufferNumber[i],6]/1000)

    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)
    DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    XMaxDiode = ParticleStatsWv[:,4]
    XMinDiode = ParticleStatsWv[:,3]
    YMinDiode = ParticleStatsWv[:,5]
    YMaxDiode = ParticleStatsWv[:,6]
    
    MaxDiameterBG= np.sqrt(DiameterBGx**2 + DiameterBGy**2)   
    AspectRatioBG = DiameterBGy / DiameterBGx
    AspectRatio = (ParticleStatsWv[:,6] - ParticleStatsWv[:,5]) / (ParticleStatsWv[:,4] - ParticleStatsWv[:,3])
    Channel = ParticleTimesWv[:,4] 
    #ImageTimes = np.round(ParticleTimeSeconds, decimals = 3)    
    ImageID = ParticleTimesWv[:,3]
  
    return ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode, MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID

#____________________________________________________________________________
# Load house keeping from 2DS oasis file  


def Load2DS_h5_HK(ImagePath, FileName):
    TAS = 0.2
    TAS_um = TAS *1E6
    ProbeResolution = 10 
    Data_h5 = h5py.File(ImagePath+FileName, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'],dtype =float)
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'],dtype =float)
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'],dtype =float)
    HouseKPcktMat = np.array(Data_h5['HouseKPcktMat'],dtype =float)
    BufferNumber = ParticleTimesWv[:,6]
    Data_h5.close()

    HouseKPcktMat[HouseKPcktMat==4294967295.] = np.nan

    
    HElement0V = HouseKPcktMat[:,2]*0.00244140625
    HElement64V = HouseKPcktMat[:,3]*0.00244140625
    HElement127V = HouseKPcktMat[:,4]*0.00244140625
    VElement0V = HouseKPcktMat[:,5]*0.00244140625
    VElement64V = HouseKPcktMat[:,6]*0.00244140625
    VElement127V = HouseKPcktMat[:,7]*0.00244140625
    HElement21V = HouseKPcktMat[:,26]*0.00244140625
    HElement42V = HouseKPcktMat[:,27]*0.00244140625
    HElement85V = HouseKPcktMat[:,28]*0.00244140625
    HElement106V = HouseKPcktMat[:,29]*0.00244140625
    VElement21V = HouseKPcktMat[:,30]*0.00244140625
    VElement42V = HouseKPcktMat[:,31]*0.00244140625
    VElement85V = HouseKPcktMat[:,32]*0.00244140625
    VElement106V = HouseKPcktMat[:,33]*0.00244140625
 
    ClockCounts = HouseKPcktMat[:,52]*(2**16) + HouseKPcktMat[:,53]
    ClockSeconds = ClockCounts * ProbeResolution / TAS_um # Probe clock time in seconds
    
    plt.plot(VElement0V,markersize=1)
    plt.plot(VElement64V,markersize=1)
    plt.plot(VElement127V,markersize=1)
    plt.plot(VElement21V,markersize=1)
    plt.plot(VElement42V,markersize=1)
    plt.plot(VElement85V,markersize=1)
    plt.plot(VElement106V,markersize=1)
    
    return ClockSeconds
#____________________________________________________________________________
# Load position/time data from high speed stage


def LoadStageLog(StagePath,StageName,RoundedInteger):
    
    StageLog = pd.read_csv(StagePath+StageName)
    dz = StageLog['dz']
    dz = dz[0] 
    StagePosition = StageLog['zmax']
    StagePosition = np.array(StagePosition[2:].astype(np.float))
    StageTime = StageLog['zmin']
    StageTime = np.array(StageTime[2:])
    StageTime = [datetime.datetime.strptime(StageTime[x], '%Y-%m-%d %H:%M:%S') for x in range(len(StageTime[2:]))]
    StageSecondsSinceMidnight = [(StageTime[x] - StageTime[x].replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds() for x in range(len(StageTime))]
    StageSecondsSinceMidnight = np.array(StageSecondsSinceMidnight)
    #StagePositon got rounded to integer. Need to add dz to odd indexes
    if RoundedInteger == 1 : 
        StagePosition[1::2] += dz
    
    return StageTime, StagePosition, StageSecondsSinceMidnight

#____________________________________________________________________________


def Ros300_vchannel():
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/FilteredImages/'
    ImageFileName = 'Export_base191007135317.h5'
    StatsFileName = 'base191007135317.h5'
    Z0 = 80
    DiffractionCalcs = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/Diffraction/'
    D0 = 270
    Lambda= 0.658 # laser wavelength
    PixelSize = 10
    ArrayWidth = 128
    
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/'
    StageName = 'Stagelog.txt' 
    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode,MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID = Load2DS_h5(ImagePath, StatsFileName)
    
    StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    count= 0
    AnalogFlag = np.zeros(len(ImageID))
    Squash_MeanXY_BG = np.zeros(len(ImageID)) * np.nan
    Squash_Area_BG = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGx = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGy = np.zeros(len(ImageID)) * np.nan
    Squash_Circularity_Filled = np.zeros(len(ImageID)) * np.nan
    Squash_MeanXY = np.zeros(len(ImageID)) * np.nan
    Squash_Area_Filled = np.zeros(len(ImageID)) * np.nan
    
    
    MeanX = (XMaxDiode + XMinDiode)/2
    MeanY = (YMaxDiode + YMinDiode)/2
    
        
    for i in range(len(ImageID)):
        if (DiameterBGy[i] > 300) and (DiameterBGy[i] < 1280) and (DiameterBGx[i] > 150) and (DiameterBGx[i] < 700) and (XMinDiode[i] > 60) and (XMaxDiode[i] < 110) and (Channel[i] == 1) and (AspectRatio[i]<5): 
            
            if (count % 1 == 0) :
                #print(count)
                #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
                
                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                Image[Image == 0 ]= 1
                Image[Image == 255 ] = 0
                
                if np.nansum(Image[int(MeanX[i]):,int(MeanY[i]):]) < np.nansum(Image[int(MeanX[i]):,:int(MeanY[i])]) : 
                    count+=1
                    #if (count % 1 == 0) :
                        #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
                    
                    SquashedImage = SquashImageY_2_5(Image)
                    Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                    #PlotOAPImage(SquashedImage,10, ImageSavePath,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]))
                    AnalogFlag[i] = 1
    
    print(count)
    DiameterBGx[AnalogFlag != 1] = np.nan
    AspectRatioBG[AnalogFlag != 1] = np.nan
    AspectRatio[AnalogFlag != 1] = np.nan
    ParticleTimeSeconds[AnalogFlag != 1] = np.nan

    SquashAspectRatioBG = Squash_DiameterBGy / Squash_DiameterBGx
    
    # Load analog data 
    if DiffractionCalcs == 1 : 
        
        PixelRatio = 0.57/10
        #PixelRatio = 1
        
        DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/Diffraction/'
        DiffractionFileName = 'Ros300_.h5'
        #DiffractionFileName = 'Ros300_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
#        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
#        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
#        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (PixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()      
    
    
    if 1==1: 
        #fig=plt.figure(figsize=(10,10))
        #plt.rcParams.update({'font.size': 12})
        plt.subplot(3, 2, 5)
        plt.title('ROS300')
        plt.plot(ImageStagePosition,Squash_MeanXY_BG,'o', markersize = 2, label = 'Analog')
        plt.plot(AST_Zarray_mm,AST_MeanXY_BG1, label = 'AST')
        #plt.legend()
        plt.xlabel('Z, mm')    
        plt.ylabel('Diameter, μm')
        plt.xlim([-35, 35])
        plt.ylim([200, 350])
        plt.subplot(3, 2, 6)
        plt.title('ROS300')
        plt.plot(ImageStagePosition,Squash_Area_BG,'o', markersize = 2)
        plt.plot(AST_Zarray_mm,ASTArea1_BG)
        #plt.xlabel('Z, mm')    
        plt.ylabel('Number of pixels')
        plt.xlabel('Z, mm')
        plt.xlim([-35, 35])
        plt.ylim([300, 550])
#        plt.subplot(3, 1, 3)
#        plt.plot(ImageStagePosition,Squash_DiameterBGy,'o', markersize = 2, label = 'x')
#        plt.plot(ImageStagePosition,Squash_DiameterBGx,'+', markersize = 2, label = 'y')
#        plt.legend()
#        plt.ylim([200, 400])
#        plt.xlabel('Z, mm')    
#        plt.ylabel('Diameter') 
#        plt.savefig(ImagePath+'DvsZ.png',dpi=200)
#        plt.close(fig)
    
     

#____________________________________________________________________________
# Load single 2DS image an use it as a mask for the diffraction calculations
def SingleImageMask(): 
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/'
    #ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/FilteredImages/'
    #ImageFileName = 'Export_base191007135317.h5'
    #ParticleBufferTimeS = 50490.208
    #ImageID = 61877
    #Prefix = 'Ros300_2DSimage'
       
#    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/'
#    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/Diffraction/'
#    ImageFileName = 'Export_base191007141713.h5'
#    ParticleBufferTimeS = 51919.229
#    ImageID = 22272
#    Prefix ='Ros250_2DSimage'
    
#    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/Output/'
#    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/Output/Diffraction/'
#    ImageFileName = 'Export_base191008092103.h5'
#    ParticleBufferTimeS = 34197.376
#    ImageID = 60956  
#    ParticleBufferTimeS = 34193.676
#    ImageID = 60133 
#    Prefix ='Ros85_2DSimage'
    
#    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/'
#    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/Diffraction/'
#    ImageFileName = 'Export_base191007100707.h5'
#    StatsFileName = 'base191007100707.h5'
#    #ParticleBufferTimeS = 36934.449
#    #ImageID = 39558
#    ParticleBufferTimeS = 36930.364
#    ImageID = 39208 
#    Prefix ='Ros250_2DSimageH'
    
    
#    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros300_run2/20191125162125/Output/'
#    FileName = 'Imagefile_1CIP Grayscale_20191125162126.h5'
#    ImageFileName = 'Export_Imagefile_1CIP Grayscale_20191125162126.h5'
#    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros300_run2/'
    #ParticleBufferTimeS = 59173.01 # backward
    #ImageID = 12909
    #Prefix ='Ros300_CIP15_1'
#    ParticleBufferTimeS = 59174.039 # forward
#    ImageID = 12955
#    Prefix ='Ros300_CIP15_2'
#    PixelSize = 15 
#    ArrayWidth = 64
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros250_run3/20191125163508/Output/'
    FileName = 'Imagefile_1CIP Grayscale_20191125163508.h5'
    ImageFileName = 'Export_Imagefile_1CIP Grayscale_20191125163508.h5'
    ImageSavePath = ImagePath
    PixelSize = 15 
    ArrayWidth = 64
    ParticleBufferTimeS = 59993.086
    ImageID = 20560
    Prefix ='Ros250_CIP15'
    #ParticleBufferTimeS = 59996.059
    #ImageID = 20642
    #Prefix ='Ros250_CIP15_2'
    #Lambda= 0.658 # laser wavelength
    
    
    
    
    
    Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS,ImageID, PixelSize, ArrayWidth)
    #Image[Image == 0 ]= 1
    #Image[Image == 255 ] = 0
#    Image[90:,:] = 0 # Just needed for 'Ros250_2DSimageH'
    
    SquashedImage = SquashImageY_4_GS(Image)
    SquashedImage[SquashedImage == 255] = 12.5 
    SquashedImage[SquashedImage == 1] = 37.5 
    SquashedImage[SquashedImage == 2] = 62.5
    SquashedImage[SquashedImage == 3] = 87.5
    
    #SquashedImage = SquashImageY_2_5(Image)
    #SquashedImage = SquashImageY_3_5(Image)
    #SquashedImage = SquashImageY_3_5(Image)
    #SquashedImage = SquashImageY_2(Image)
    
    #zoom image to 1 um pixel size 
    ZoomArray = ndimage.zoom(SquashedImage, PixelSize , order=3)    
    #BinaryImage = (np.where(ZoomArray>0.5, 0, 1))
    BinaryImage = (np.where(ZoomArray>=50, 0, 1))
    
    MaskArray = np.ones((2048,2048)) 
#    #Put shape on to common grid
    shape=np.shape(BinaryImage)
    xlen=shape[0]
    ylen=shape[1]   
    MaskArray[1024-int(xlen/2) : 1024+int(xlen/2+0.5), 1024-int(ylen/2+0.5) :1024+int(ylen/2) ]= BinaryImage

    #Do AST calcs
    PixelSizeInput = 1 
    PixelSizeOutput = 1 
    
    ShapeVsZ(ImageSavePath,Prefix,2,MaskArray,PixelSizeInput,PixelSizeOutput)
    
    
    
    #ShapeVsZ(DiffractionPath,'Ros300_',2,SourceImage)
#    return ZoomArray

#____________________________________________________________________________

#Use tif image for diffraction calculations

def AnalogAST():
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/ros300x10.tif'
    #DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros300_tas_02/Output/Diffraction/'
    #SubsetFlag =0
    
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/rosette118x4.tif'
    #DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/'
    #PixelSizeInput = 0.89 # um
    #PixelSizeOutput = 0.89 # um
    
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/ROS250.tif'
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/'
    PixelSizeInput = 1.25 # um
    PixelSizeOutput = 1.25 # um
    
    SubsetFlag = 1 

    
    SourceImage = LoadTiffImage(ImagePath,SubsetFlag)
    ShapeVsZ(DiffractionPath,'Ros250_',2,SourceImage, PixelSizeInput,PixelSizeOutput)

#____________________________________________________________________________

def LoadTiffImage(ImagePath,SubsetFlag):
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/ros300x10.tif'
    #ImagePath = 'C:/Users/Admin TEMP/Documents/Diffraction/CPI_images/Rosette/CP20080401_204059_387_1.jpg'
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/rosette118x4.tif'
    
    #ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/ROS250.tif'
    
    MaskArray = np.ones((2048,2048)) 
    
   
    import_image = (ndimage.imread(ImagePath,flatten=True).astype(float))
    if SubsetFlag == 1 :
        #import_image= import_image[600:900, 1050:1250] #Subset rosette118x4.tif
        import_image= import_image[600:900, 950:1200] #Subset ROS250.tif
        
        
        
    val = filters.threshold_otsu(import_image)
    BinaryImage = (np.where(import_image<val, 0, 1))
    

    
    
    

#    BoxStats = regionprops(BinaryImage, cache=False)
#    Boxbbox = [r.bbox for r in BoxStats]
#    Diametery = Boxbbox[0][2] - Boxbbox[0][0]
#    Diameterx = Boxbbox[0][3] - Boxbbox[0][1]
#    BoxMeanXY = (Boxbbox[0][2] - Boxbbox[0][0] + Boxbbox[0][3] - Boxbbox[0][1]) / 2
#    
#    #select largest particle in images
#    labels_max= SelectLargestParticle(BinaryImage)
#    stats = regionprops(labels_max, cache=False)
#    bbox = [r.bbox for r in stats]                     
#    DiameterBGy = Boxbbox[0][2] - Boxbbox[0][0]   
#    DiameterBGx = Boxbbox[0][3] - Boxbbox[0][1]
#    MeanXY_BG = (bbox[0][2] - bbox[0][0] + bbox[0][3] - bbox[0][1]) / 2   
#        
#    #fill internal voids largest particle in image
#    FilledImage = ndimage.morphology.binary_fill_holes(labels_max).astype(int) # fill internal voids for circularity 
#    FilledStats = regionprops(FilledImage, cache=False) 
#    bbox = [r.bbox for r in stats]                     
#    Diametery_BGFilled = Boxbbox[0][2] - Boxbbox[0][0]   
#    Diameterx_BGFilled = Boxbbox[0][3] - Boxbbox[0][1]
#    MeanXY_BGFilled = (bbox[0][2] - bbox[0][0] + bbox[0][3] - bbox[0][1]) / 2
#    Circularity_BGFilled = FilledStats[0].perimeter** 2 / (4 * np.pi * FilledStats[0].area )

#    shape=np.shape(BinaryImage)
#    xlen=shape[0]
#    ylen=shape[1]
#    COMx = 100
#    COMy = 100
#    
#    MaskArray[1024-int(xlen/2) : 1024+int(xlen/2+0.5), 1024-int(ylen/2+0.5) :1024+int(ylen/2) ]= BinaryImage
    
    MaskArray = np.ones((2048,2048)) 
    
    #Put shape on to common grid
    shape=np.shape(BinaryImage)
    xlen=shape[0]
    ylen=shape[1]   
    MaskArray[1024-int(xlen/2) : 1024+int(xlen/2+0.5), 1024-int(ylen/2+0.5) :1024+int(ylen/2) ]= BinaryImage
    

    #return BinaryImage, Diametery, Diameterx, BoxMeanXY, DiameterBGy, DiameterBGx, MeanXY_BG, Diametery_BGFilled, Diameterx_BGFilled, MeanXY_BGFilled, Circularity_BGFilled
    return MaskArray

#__________________________________________________________________________________
#BinaryImage = np.array(Image)
#BinaryImage[BinaryImage==0] = 1
#BinaryImage[BinaryImage==255] = 0

#Squash image by factor 2.5 


def SquashImageY_2_5(BinaryImage): 
    
    #Stretch image by duplicating each y pixel. 
    BinaryImage = np.array(BinaryImage, dtype=float)
    StretchImage= np.repeat(BinaryImage,2,axis=1)
    Pad= np.zeros([128,5])
    StretchImage = np.concatenate((Pad,StretchImage,Pad),axis = 1)
    #SquashImage = (BinaryImage[:,:-3:3]+BinaryImage[:,1:-2:3]+0.5*BinaryImage[:,2:-1:3])/2.5
    TmpSum = np.cumsum(StretchImage, axis =1)
    Tmp = TmpSum[:,4::5]
    SquashImage = Tmp[:,1::1] - Tmp[:,0:-1:1]
    SquashImage /= 5 
    
    SquashImage[SquashImage < 0.5] = int(0)
    SquashImage[SquashImage >= 0.5] = int(1)
    SquashImage= np.array(SquashImage, dtype=int)
    
    return SquashImage#, StretchImage, TmpSum, Tmp


#__________________________________________________________________________________
#BinaryImage = np.array(Image)
#BinaryImage[BinaryImage==0] = 1
#BinaryImage[BinaryImage==255] = 0

#Squash image by factor 2.5 


def SquashImageY_3_5(BinaryImage): 
    
    BinaryImage = np.array(BinaryImage, dtype=float)
    StretchImage= np.repeat(BinaryImage,2,axis=1)
    Pad= np.zeros([128,7])
    StretchImage = np.concatenate((Pad,StretchImage,Pad),axis = 1)
    #SquashImage = (BinaryImage[:,:-3:3]+BinaryImage[:,1:-2:3]+0.5*BinaryImage[:,2:-1:3])/2.5
    TmpSum = np.cumsum(StretchImage, axis =1)
    Tmp = TmpSum[:,6::7]
    SquashImage = Tmp[:,1::1] - Tmp[:,0:-1:1]
    SquashImage /= 7 
    
    SquashImage[SquashImage < 0.5] = int(0)
    SquashImage[SquashImage >= 0.5] = int(1)
    SquashImage= np.array(SquashImage, dtype=int)
    
    return SquashImage#, StretchImage, TmpSum, Tmp


#__________________________________________________________________________________
#BinaryImage = np.array(Image)
#BinaryImage[BinaryImage==0] = 1
#BinaryImage[BinaryImage==255] = 0

#Squash image by factor 2


def SquashImageY_2(BinaryImage): 
    BinaryImage = np.array(BinaryImage, dtype=float)   
    SquashImage = (BinaryImage[:,:-1:2] + BinaryImage[:,1::2])/ 2  
    SquashImage[SquashImage < 0.5] = int(0)
    SquashImage[SquashImage >= 0.5] = int(1)
    SquashImage= np.array(SquashImage, dtype=int)
    
    return SquashImage#, StretchImage, TmpSum, Tmp



#__________________________________________________________________________________

def PlotOAPImage(Image,PixelSize, SavePath, ImageTimes, ImageID, ImageZD, ImageZ):
    
    shape = np.shape(Image)
    xlen = shape[0]
    ylen = shape[1]
    
    x = np.arange(-PixelSize/2,xlen*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,ylen*PixelSize+PixelSize/2, PixelSize)
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    plt.title(str(ImageTimes)+'s, ID = ' + str(ImageID)+', Z = '+str(ImageZ)+'mm, Zd = '+str(np.round(ImageZD,decimals =3)))
    plt.pcolormesh(y,x,Image, cmap='Greys')
    plt.xlabel('x, μm')
    plt.ylabel('y, μm')
    plt.ylim([0,1280])
    plt.xlim([0,1280])
    plt.savefig(SavePath+str(ImageTimes)+'s_ID_' + str(ImageID)+'Squashed.png',dpi=200)
    plt.close(fig)
    
#__________________________________________________________________________________
# Plot a OAP greyscale image


def PlotOAPImageGS(Image,PixelSize, SavePath, ImageTimes, ImageID, ImageZD, ImageZ):
    
    shape = np.shape(Image)
    xlen = shape[0]
    ylen = shape[1]
    
    x = np.arange(-PixelSize/2,xlen*PixelSize+PixelSize/2, PixelSize)
    y = np.arange(-PixelSize/2,ylen*PixelSize+PixelSize/2, PixelSize)
    
    
    fig=plt.figure(figsize=(10,10))
    plt.rcParams.update({'font.size': 12})
    plt.title(str(ImageTimes)+'s, ID = ' + str(ImageID)+', Z = '+str(ImageZ)+'mm, Zd = '+str(np.round(ImageZD,decimals =3)))
    plt.pcolormesh(y,x,Image, cmap='gist_stern_r', vmin=0, vmax=3)
    plt.xlabel('x, μm')
    plt.ylabel('y, μm')
    plt.ylim([0,960])
    plt.xlim([0,960])
    plt.savefig(SavePath+str(ImageTimes)+'s_ID_' + str(ImageID)+'Squashed.png',dpi=200)
    plt.close(fig)




#__________________________________________________________________________________


def Ros250_vchannel():
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/FilteredImages/'
    ImageFileName = 'Export_base191007141713.h5'
    StatsFileName = 'base191007141713.h5'
    Z0 = 80
    DiffractionCalcs = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/Diffraction/'
    D0 = 270
    Lambda= 0.658 # laser wavelength
    PixelSize = 10
    ArrayWidth = 128
    
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/'
    StageName = 'Stagelog.txt' 
    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode,MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID = Load2DS_h5(ImagePath, StatsFileName)
    StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    count= 0
    AnalogFlag = np.zeros(len(ImageID))
    Squash_MeanXY_BG = np.zeros(len(ImageID)) * np.nan
    Squash_Area_BG = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGx = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGy = np.zeros(len(ImageID)) * np.nan
    Squash_Circularity_Filled = np.zeros(len(ImageID)) * np.nan
    Squash_MeanXY = np.zeros(len(ImageID)) * np.nan
    Squash_Area_Filled = np.zeros(len(ImageID)) * np.nan
    
    for i in range(len(ImageID)):
        if (DiameterBGy[i] > 400) and (DiameterBGy[i] < 2000) and (DiameterBGx[i] > 150) and (DiameterBGx[i] < 700) and (XMinDiode[i] > 60) and (XMaxDiode[i] < 110) and (Channel[i] == 1):# and (AspectRatio[i]<5): 
            count+=1
            if (count % 1 == 0) :
                
                #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                Image[Image == 0 ]= 1
                Image[Image == 255 ] = 0
                SquashedImage = SquashImageY_3_5(Image)
                Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                
#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] > 500)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/BigA/')
#                
#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] < 450)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/SmallA/')
            
                
                
                AnalogFlag[i] = 1
    print(count)
    DiameterBGy[AnalogFlag != 1] = np.nan
    DiameterBGx[AnalogFlag != 1] = np.nan
    AspectRatioBG[AnalogFlag != 1] = np.nan
    AspectRatio[AnalogFlag != 1] = np.nan
    SquashAspectRatioBG = Squash_DiameterBGy / Squash_DiameterBGx
    
    
    # Load analog data 
    if DiffractionCalcs == 1 : 
        
        #PixelRatio = 0.57/10
        PixelRatio = 1
        
        DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/Diffraction/'
        DiffractionFileName = 'Ros250_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
#        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
#        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
#        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (PixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()      
    
    
    if 1==1: 
        fig=plt.figure(figsize=(10,10))
        plt.rcParams.update({'font.size': 12})
        plt.subplot(3, 1, 1)
        plt.plot(ImageStagePosition,Squash_MeanXY_BG,'o', markersize = 2, label = 'Analog')
        plt.plot(AST_Zarray_mm,AST_MeanXY_BG1, label = 'AST')
        plt.legend()
        #plt.xlabel('Z, mm')    
        plt.ylabel('Diameter, μm')
        plt.xlim([-50, 50])
        plt.ylim([200, 400])
        plt.subplot(3, 1, 2)
        plt.plot(ImageStagePosition,Squash_Area_BG,'o', markersize = 2)
        plt.plot(AST_Zarray_mm,ASTArea1_BG)
        #plt.xlabel('Z, mm')    
        plt.ylabel('Number of pixels')
        plt.xlim([-50, 50])
        plt.ylim([300, 700])
        plt.subplot(3, 1, 3)
        plt.plot(ImageStagePosition,DiameterBGy,'o', markersize = 2, label = 'y')
        plt.plot(ImageStagePosition,Squash_DiameterBGy,'o', markersize = 2, label = 'y squashed')
        plt.plot(ImageStagePosition,Squash_DiameterBGx,'+', markersize = 2, label = 'x')
        plt.legend()
        plt.xlim([-50, 50])
        #plt.ylim([200, 400])
        #plt.plot(AST_Zarray_mm,ASTCircularity1_filled)
        #plt.plot(ImageStagePosition,SquashAspectRatioBG,'o')
        plt.xlabel('Z, mm')    
        plt.ylabel('Diameter') 
        plt.savefig(ImagePath+'DvsZ.png',dpi=200)
        plt.close(fig)



#__________________________________________________________________________________


def Ros85_vchannel():
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/Output/'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/Output/FilteredImages/'
    ImageFileName = 'Export_base191008092103.h5'
    StatsFileName = 'base191008092103.h5'
    Z0 = 80
    DiffractionCalcs = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191007_ros250_tas_02/Output/Diffraction/'
    D0 = 85
    Lambda= 0.658 # laser wavelength
    PixelSize = 10
    ArrayWidth = 128
    
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/'
    StageName = 'Stagelog.txt' 
    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode,MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID = Load2DS_h5(ImagePath, StatsFileName)
    
    StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    count= 0
    AnalogFlag = np.zeros(len(ImageID))
    Squash_MeanXY_BG = np.zeros(len(ImageID)) * np.nan
    Squash_Area_BG = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGx = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGy = np.zeros(len(ImageID)) * np.nan
    Squash_Circularity_Filled = np.zeros(len(ImageID)) * np.nan
    Squash_MeanXY = np.zeros(len(ImageID)) * np.nan
    Squash_Area_Filled = np.zeros(len(ImageID)) * np.nan
    
    for i in range(len(ImageID)):
        if (DiameterBGy[i] > 50) and (DiameterBGy[i] < 200) and (DiameterBGx[i] > 50) and (DiameterBGx[i] < 200) and (XMinDiode[i] > 53) and (XMaxDiode[i] < 67) and (Channel[i] == 1):# and (AspectRatio[i]<5): 
            count+=1
            if (count % 1 == 0) :
                
                #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                Image[Image == 0 ]= 1
                Image[Image == 255 ] = 0
                SquashedImage = SquashImageY_2(Image)
                Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                #PlotOAPImage(SquashedImage,10, ImageSavePath,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]))

#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] > 500)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/BigA/')
#                
#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] < 450)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/SmallA/')
                AnalogFlag[i] = 1
    print(count)
    DiameterBGy[AnalogFlag != 1] = np.nan
    DiameterBGx[AnalogFlag != 1] = np.nan
    AspectRatioBG[AnalogFlag != 1] = np.nan
    AspectRatio[AnalogFlag != 1] = np.nan
    SquashAspectRatioBG = Squash_DiameterBGy / Squash_DiameterBGx
    
    
    # Load analog data 
    if DiffractionCalcs == 1 : 
        
        #PixelRatio = 0.57/10
        PixelRatio = 1
        
        DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191008_ros85_tas02/Output/Diffraction/'
        DiffractionFileName = 'Ros85_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
#        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
#        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
#        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (PixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()      
    
    
    if 1==1: 
        #fig=plt.figure(figsize=(10,10))
        #plt.rcParams.update({'font.size': 12})
        plt.subplot(3, 2, 1)
        plt.title('ROS85')
        plt.plot(ImageStagePosition,Squash_MeanXY_BG,'o', markersize = 2, label = 'Analog')
        plt.plot(AST_Zarray_mm,AST_MeanXY_BG1, label = 'AST')
        
        #plt.xlabel('Z, mm')    
        plt.ylabel('Diameter, μm')
        plt.xlim([-35, 35])
        #plt.ylim([200, 400])
        plt.subplot(3, 2, 2)
        #plt.title('ROS85')
        plt.plot(ImageStagePosition,Squash_Area_BG,'o', markersize = 2, label = 'Analog')
        plt.plot(AST_Zarray_mm,ASTArea1_BG,label = 'AST')
        #plt.xlabel('Z, mm')    
        plt.ylabel('Number of pixels')
        #plt.ylabel('A$_{50 - 100}$')
        plt.xlim([-35, 35])
        plt.legend()
        #plt.ylim([300, 700])
#        plt.subplot(3, 1, 3)
#        plt.plot(ImageStagePosition,DiameterBGy,'o', markersize = 2, label = 'y')
#        plt.plot(ImageStagePosition,Squash_DiameterBGy,'o', markersize = 2, label = 'y squashed')
#        plt.plot(ImageStagePosition,Squash_DiameterBGx,'+', markersize = 2, label = 'x')
#        plt.legend()
#        plt.xlim([-50, 50])
#        plt.xlabel('Z, mm')    
#        plt.ylabel('Diameter') 
#        plt.savefig(ImagePath+'DvsZ.png',dpi=200)
#        plt.close(fig)
# plt.ylabel('A$_{50 - 100}$ / A$_{'+Low+' - 100}$')

#__________________________________________________________________________________


def Ros250_hchannel():
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/FilteredImages/'
    ImageFileName = 'Export_base191007100707.h5'
    StatsFileName = 'base191007100707.h5'
    Z0 = 80
    DiffractionCalcs = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/Diffraction/'
    D0 = 270
    Lambda= 0.658 # laser wavelength
    PixelSize = 10
    ArrayWidth = 128
    
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/'
    StageName = 'Stagelog.txt' 
    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode,MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID = Load2DS_h5(ImagePath, StatsFileName)

    StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    count= 0
    AnalogFlag = np.zeros(len(ImageID))
    Squash_MeanXY_BG = np.zeros(len(ImageID)) * np.nan
    Squash_Area_BG = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGx = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGy = np.zeros(len(ImageID)) * np.nan
    Squash_Circularity_Filled = np.zeros(len(ImageID)) * np.nan
    Squash_MeanXY = np.zeros(len(ImageID)) * np.nan
    Squash_Area_Filled = np.zeros(len(ImageID)) * np.nan
    
    for i in range(len(ImageID)):
        #if (DiameterBGy[i] > 300) and (DiameterBGy[i] < 600) and (DiameterBGx[i] > 150) and (DiameterBGx[i] < 500) and (XMinDiode[i] > 0) and (XMaxDiode[i] < 90) and (Channel[i] == 0):# and (AspectRatio[i]<5): 
        #if (DiameterBGx[i] > 150) and (XMinDiode[i] > 40) and (XMaxDiode[i] > 60) and (DiameterBGy[i] > 300) and (DiameterBGy[i] < 1000) and (Channel[i] == 0):# and (AspectRatio[i]<5): 
        if (DiameterBGx[i] > 150) and (XMinDiode[i] > 40) and (XMinDiode[i] < 65) and (DiameterBGy[i] > 300) and (DiameterBGy[i] < 1000) and (Channel[i] == 0):# and (AspectRatio[i]<5): 
       
            count+=1
            if (count % 1 == 0) :
                
                #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                Image[Image == 0 ]= 1
                Image[Image == 255 ] = 0
                SquashedImage = SquashImageY_2(Image)
                Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                
#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] > 500)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/BigA/')
#                
#                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] < 450)) :
#                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/SmallA/')
            
                
                
            AnalogFlag[i] = 1
    print(count)
    DiameterBGy[AnalogFlag != 1] = np.nan
    DiameterBGx[AnalogFlag != 1] = np.nan
    AspectRatioBG[AnalogFlag != 1] = np.nan
    AspectRatio[AnalogFlag != 1] = np.nan
    SquashAspectRatioBG = Squash_DiameterBGy / Squash_DiameterBGx
    
#    plt.plot(ImageStagePosition, AspectRatio,'o')
#    plt.plot(ImageStagePosition, AspectRatioBG,'o')
#    
#    
    # Load analog data 
    if DiffractionCalcs == 1 : 
        
        #PixelRatio = 0.57/10
        PixelRatio = 1
        
        DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/Diffraction/'
        #DiffractionFileName = 'Ros250_2DSimageH_36934.h5'       
        DiffractionFileName = 'Ros250_2DSimageH_39208.h5'
        
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
#        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
#        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
#        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (PixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()      
    
    
    if 1==1: 
        #fig=plt.figure(figsize=(10,10))
        #plt.rcParams.update({'font.size': 12})
        plt.subplot(3, 2, 3)
        plt.title('ROS250')
        plt.plot(ImageStagePosition,Squash_MeanXY_BG,'o', markersize = 2, label = 'Analog')
        plt.plot(AST_Zarray_mm,AST_MeanXY_BG1, label = 'AST')
        #plt.legend()
        #plt.xlabel('Z, mm')    
        plt.ylabel('Diameter, μm')
        plt.xlim([-50, 50])
        plt.ylim([150, 300])
        plt.subplot(3, 2, 4)
        plt.title('ROS250')
        plt.plot(ImageStagePosition,Squash_Area_BG,'o', markersize = 2)
        plt.plot(AST_Zarray_mm,ASTArea1_BG)
        #plt.xlabel('Z, mm')    
        plt.ylabel('Number of pixels')
        plt.xlim([-35, 35])
        plt.ylim([150, 450])
#        plt.subplot(3, 1, 3)
#        plt.plot(ImageStagePosition,DiameterBGy,'o', markersize = 2, label = 'y')
#        plt.plot(ImageStagePosition,Squash_DiameterBGy,'o', markersize = 2, label = 'y squashed')
#        plt.plot(ImageStagePosition,Squash_DiameterBGx,'+', markersize = 2, label = 'x')
#        plt.legend()
#        plt.xlim([-50, 50])
#        plt.xlabel('Z, mm')    
#        plt.ylabel('Diameter') 
#        plt.savefig(ImagePath+'DvsZ.png',dpi=200)
#        plt.close(fig)


#__________________________________________________________________________________
def Plot3Rosettes2DS(): 
    SavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/'
    fig=plt.figure(figsize=(12,14))
    plt.rcParams.update({'font.size': 12})
    Ros85_vchannel()
    Ros250_hchannel()
    Ros300_vchannel()
    plt.savefig(SavePath+'DvsZ.png',dpi=200)
    plt.close(fig)
     
#__________________________________________________________________________________

def Dots125_vchannel():
    
    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191011_125dots/Output/'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191011_125dots/Output/FilteredImages/'
    ImageFileName = 'Export_base191011102407.h5'
    StatsFileName = 'base191011102407.h5'
    Z0 = 80
    DiffractionCalcs = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/Diffraction/'
    D0 = 125
    Lambda= 0.658 # laser wavelength
    PixelSize = 10
    ArrayWidth = 128
    
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Vchannel/20191011_125dots/'
    StageName = 'Stagelog.txt' 
    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    ParticleBufferTimeS, ParticleTimeSeconds, MeanXY_BG, DiameterBGy, DiameterBGx, XMaxDiode,XMinDiode, YMaxDiode, YMinDiode,MaxDiameterBG, AspectRatioBG, AspectRatio, Channel, ImageID = Load2DS_h5(ImagePath, StatsFileName)

    StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    count= 0
    AnalogFlag = np.zeros(len(ImageID))
    Squash_MeanXY_BG = np.zeros(len(ImageID)) * np.nan
    Squash_Area_BG = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGx = np.zeros(len(ImageID)) * np.nan
    Squash_DiameterBGy = np.zeros(len(ImageID)) * np.nan
    Squash_Circularity_Filled = np.zeros(len(ImageID)) * np.nan
    Squash_MeanXY = np.zeros(len(ImageID)) * np.nan
    Squash_Area_Filled = np.zeros(len(ImageID)) * np.nan
    
    for i in range(len(ImageID)):
        #if (DiameterBGy[i] > 300) and (DiameterBGy[i] < 600) and (DiameterBGx[i] > 150) and (DiameterBGx[i] < 500) and (XMinDiode[i] > 0) and (XMaxDiode[i] < 90) and (Channel[i] == 0):# and (AspectRatio[i]<5): 
        #if (DiameterBGx[i] > 150) and (XMinDiode[i] > 40) and (XMaxDiode[i] > 60) and (DiameterBGy[i] > 300) and (DiameterBGy[i] < 1000) and (Channel[i] == 0):# and (AspectRatio[i]<5): 
        if (DiameterBGx[i] > 50) and (DiameterBGx[i] < 200) and (Channel[i] == 1):# and (AspectRatio[i]<5): 
       
            count+=1
#            if (count % 1 == 0) :
#                
#                #PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath)
#                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
#                Image[Image == 0 ]= 1
#                Image[Image == 255 ] = 0
#                #SquashedImage = SquashImageY_2(Image)
#                #Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
#                
##                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] > 500)) :
##                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/BigA/')
##                
##                if ((ImageStagePosition[i] > 15 ) and (Squash_Area_BG[i] < 450)) :
##                    PlotImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]), ImageSavePath+'/SmallA/')
#            
#                
                
            AnalogFlag[i] = 1
    print(count)
    DiameterBGy[AnalogFlag != 1] = np.nan
    DiameterBGx[AnalogFlag != 1] = np.nan
    AspectRatioBG[AnalogFlag != 1] = np.nan
    AspectRatio[AnalogFlag != 1] = np.nan
    
    if 1==1: 
        fig=plt.figure(figsize=(10,10))
        plt.rcParams.update({'font.size': 12})
        plt.subplot(2, 1, 1)
        plt.title('DOT150')
        plt.plot(ImageStagePosition,DiameterBGx,'o', markersize = 2, label = 'Analog')
        #plt.plot(AST_Zarray_mm,AST_MeanXY_BG1, label = 'AST')
        #plt.legend()
        #plt.xlabel('Z, mm')    
        plt.ylabel('Diameter, μm')
        plt.xlim([-50, 50])
        #plt.ylim([150, 300])
        plt.subplot(2, 1, 2)
        plt.title('DOT150')
        plt.plot(ImageStagePosition,AspectRatioBG,'o', markersize = 2)
        plt.ylabel('Aspect ratio')
        plt.xlim([-35, 35])
        #plt.ylim([150, 450])
#        plt.subplot(3, 1, 3)
#        plt.plot(ImageStagePosition,DiameterBGy,'o', markersize = 2, label = 'y')
#        plt.plot(ImageStagePosition,Squash_DiameterBGy,'o', markersize = 2, label = 'y squashed')
#        plt.plot(ImageStagePosition,Squash_DiameterBGx,'+', markersize = 2, label = 'x')
#        plt.legend()
#        plt.xlim([-50, 50])
#        plt.xlabel('Z, mm')    
#        plt.ylabel('Diameter') 
#        plt.savefig(ImagePath+'DvsZ.png',dpi=200)
#        plt.close(fig)
#    
#    
#    # Load analog data 
#    if DiffractionCalcs == 1 : 
#        
#        #PixelRatio = 0.57/10
#        PixelRatio = 1
#        
#        DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/2DS/Hchannel/20191007_ros250_tas02/Output/Diffraction/'
#        #DiffractionFileName = 'Ros250_2DSimageH_36934.h5'       
#        DiffractionFileName = 'Ros250_2DSimageH_39208.h5'
#        
#        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
#        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
#        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
#        
#        
#        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
#        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
##        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
##        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
##        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
##        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
#        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
#        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
##        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
#        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (PixelRatio**2)
#        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
#        
#        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
#        Data_h5.close()      
    
#__________________________________________________________________________________

def LoadCIP_h5(ImagePath, FileName) :
    PixelSize =15
    ArrayWidth = 64
    Data_h5 = h5py.File(ImagePath+FileName, 'r')              
    ParticleTimesWv=np.array(Data_h5['ParticleTimesWv'])
    HeaderMatrixWv=np.array(Data_h5['HeaderMatrixWv'])
    ParticleStatsWv=np.array(Data_h5['ParticleStatsWv'])
    BufferNumber = ParticleTimesWv[:,6]
    Data_h5.close()

    #Get particle buffer time in seconds from midnight. 3 decimal places
    ParticleBufferTimeS= np.zeros(len(BufferNumber))*np.nan
    #AnalogFlag = np.zeros(len(BufferNumber)) 
    
    for i in range(len(BufferNumber)):
        if BufferNumber[i] == 4294967295 :
            ParticleBufferTimeS[i] = np.nan
        else:   
            ParticleBufferTimeS[i] = 3600*HeaderMatrixWv[BufferNumber[i],3]+ 60*HeaderMatrixWv[BufferNumber[i],4] + HeaderMatrixWv[BufferNumber[i],5] + (HeaderMatrixWv[BufferNumber[i],6]/1000)

    
    ParticleTimeSeconds = ParticleTimesWv[:,0] +(ParticleTimesWv[:,1]/1E9)
    #MeanXY_BG = PixelSize+PixelSize*((ParticleStatsWv[:,12]-ParticleStatsWv[:,11] + ParticleStatsWv[:,10] - ParticleStatsWv[:,9]) /2)
    #DiameterBGy=PixelSize+PixelSize*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11])
    #DiameterBGx=PixelSize+PixelSize*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
    XMinDiode = np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31])
    XMaxDiode = np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32])
    #YMinDiode = ParticleStatsWv[:,5]
    #YMaxDiode = ParticleStatsWv[:,6]
    
    #MaxDiameterBG= np.sqrt(DiameterBGx**2 + DiameterBGy**2)   
    #AspectRatioBG = DiameterBGy / DiameterBGx
    #AspectRatio = (ParticleStatsWv[:,6] - ParticleStatsWv[:,5]) / (ParticleStatsWv[:,4] - ParticleStatsWv[:,3])
    #Channel = ParticleTimesWv[:,4] 
    #ImageTimes = np.round(ParticleTimeSeconds, decimals = 3)    
    ImageID = ParticleTimesWv[:,3]
    
    DiameterLevel0x=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32]),ParticleStatsWv[:,24])-np.minimum(np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]),ParticleStatsWv[:,23]))
    DiameterLevel1x=15+15*(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32])-np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]))
    DiameterLevel2x=15+15*(ParticleStatsWv[:,32]-ParticleStatsWv[:,31])    
    DiameterBGx=15+15*(ParticleStatsWv[:,10]-ParticleStatsWv[:,9]) #x diameter of largerst paricle in bounding box 
        
    DiameterLevel0y=15+15*(np.maximum(np.maximum(ParticleStatsWv[:,30],ParticleStatsWv[:,34]),ParticleStatsWv[:,26])-np.minimum(np.minimum(ParticleStatsWv[:,29],ParticleStatsWv[:,33]),ParticleStatsWv[:,25]))
    DiameterLevel1y=15+15*(np.maximum(ParticleStatsWv[:,30],ParticleStatsWv[:,34])-np.minimum(ParticleStatsWv[:,29],ParticleStatsWv[:,33]))
    DiameterLevel2y=15+15*(ParticleStatsWv[:,34]-ParticleStatsWv[:,33])    
    DiameterBGy=15+15*(ParticleStatsWv[:,12]-ParticleStatsWv[:,11]) #x diameter of largerst paricle in bounding box 

    DiameterLevel0= (DiameterLevel0x + DiameterLevel0y)/2
    DiameterLevel1= (DiameterLevel1x + DiameterLevel1y)/2
    DiameterLevel2= (DiameterLevel2x + DiameterLevel2y)/2  
    DiameterBG= (DiameterBGx + DiameterBGy)/2  
    
    AspectRatioLevel1 = DiameterLevel1y / DiameterLevel1x
    
    
    
    #at 25% 
    XMaxDiode = np.maximum(np.maximum(ParticleStatsWv[:,28],ParticleStatsWv[:,32]),ParticleStatsWv[:,24])
    XMinDiode = np.minimum(np.minimum(ParticleStatsWv[:,27],ParticleStatsWv[:,31]),ParticleStatsWv[:,23])
    YMinDiode = np.minimum(np.minimum(ParticleStatsWv[:,25],ParticleStatsWv[:,29]),ParticleStatsWv[:,33])
    YMaxDiode = np.maximum(np.maximum(ParticleStatsWv[:,26],ParticleStatsWv[:,30]),ParticleStatsWv[:,34])
    
    return ParticleBufferTimeS,ParticleTimeSeconds,ImageID,ParticleStatsWv, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode, AspectRatioLevel1
    
    #return ParticleBufferTimeS,ParticleTimeSeconds, MeanXY_BG, DiameterBGx, DiameterBGy, XMaxDiode, XMinDiode, YMinDiode,YMaxDiode, MaxDiameterBG,  AspectRatioBG, AspectRatio, ImageID
    
#__________________________________________________________________________________
    
#Analyse Cip15 experiment using ROS300 and high speed stage 

def CIP15_ros300():

    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros300_run2/20191125162125/Output/'
    FileName = 'Imagefile_1CIP Grayscale_20191125162126.h5'
    ImageFileName = 'Export_Imagefile_1CIP Grayscale_20191125162126.h5'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros300_run2/20191125162125/Output/FilteredImages/'
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros300_run2/'
    StageName = 'Stagelog.txt' 
    DiffractionCalcs  = 1
    DiffractionPath = StagePath
    DiffractionFileName = 'Ros300_CIP15_2.h5'
    DiffractionPixelRatio = 1/15
    
    Lambda= 0.658 # laser wavelength
    Z0 = 101
    D0 = 270
    PixelSize = 15
    ArrayWidth = 64
    
    ParticleBufferTimeS, ParticleTimeSeconds, ImageID, ParticleStatsWv, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode, AspectRatioLevel1  = LoadCIP_h5(ImagePath, FileName)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG, Zd_fromRatios = GrayscaleRatios(0, 0, ParticleStatsWv,0,'255075')

    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    #StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan
    
    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    
    AnalogFlag = np.zeros(len(ImageID))
    # Search for analog images 
        
    SquashedDiameterLevel0_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel1_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel2_BG=np.zeros(len(ImageID))*np.nan
    SquashedArea0_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea1_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea2_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea0_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea1_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea2_Filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity0_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity1_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity2_filled = np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel0=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel1=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel2=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGy1=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGx1=np.zeros(len(ImageID))*np.nan
    
    count = 0
    #for i in range(len(ImageID)):
        #if (DiameterLevel1[i] > 250) and (DiameterLevel1[i] < 600) :
        
            #count+=1
            #AnalogFlag[i] = 1
    #Select particles     
    AnalogFlag[(DiameterLevel1 > 250) & (DiameterLevel1 < 600)] = 1   
    ParticleTimeSeconds[AnalogFlag != 1] = np.nan
    #use differnce in particle element number to determine stage direction
    Index=np.where(ParticleTimeSeconds==ParticleTimeSeconds)
    Index=np.array(Index[0]).astype(float)
    IndexHigh = Index[1:] - Index[0:-1]
    IndexHigh = np.append(IndexHigh,[np.nan] )
    IndexLow = Index[1:] - Index[0:-1]
    IndexLow = np.append([np.nan],IndexLow )
    Direction = np.zeros(len(ParticleTimeSeconds))*np.nan
    Index = Index[IndexHigh > IndexLow] # > for forward, < for backward
    Index = np.array(Index).astype(int)
    Direction[Index] = 1
            
    for i in range(len(ImageID)):
        if Direction[i] == 1:
            count+=1
            if (count % 1 == 0) :
                #PlotImage_GS(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],ImageZd[i], ImageStagePosition[i], ImageSavePath)

                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                SquashedImage = SquashImageY_4_GS(Image)
                #Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                #PlotOAPImageGS(SquashedImage,PixelSize, ImageSavePath,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]))
                
                I_binned_75 = (np.where(SquashedImage>=3, 1, 0))
                I_binned_50 = (np.where(SquashedImage>=2, 1, 0))
                I_binned_25 = (np.where(SquashedImage>=1, 1, 0))

                SquashedDiameterLevel0_BG[i], SquashedArea0_BG[i], tmp, tmp, SquashedCircularity0_filled[i],SquashedBoxDiameterLevel0[i],tmp, SquashedArea0_Filled[i] = ImageParticleMoreStats(I_binned_25, PixelSize)
                SquashedDiameterLevel1_BG[i], SquashedArea1_BG[i], SquashedDiameterBGx1[i],SquashedDiameterBGy1[i], SquashedCircularity1_filled[i],SquashedBoxDiameterLevel1[i],tmp, SquashedArea1_Filled[i] = ImageParticleMoreStats(I_binned_50, PixelSize)
                SquashedDiameterLevel2_BG[i], SquashedArea2_BG[i], tmp, tmp,SquashedCircularity2_filled[i],SquashedBoxDiameterLevel2[i],tmp, SquashedArea2_Filled[i] = ImageParticleMoreStats(I_binned_75, PixelSize)
                
        
    SquashedAreaFraction0=(SquashedArea0_BG-SquashedArea1_BG)/SquashedArea0_BG # Alow notation from OShea AMT 2019   
    SquashedAreaFraction1=(SquashedArea1_BG-SquashedArea2_BG)/SquashedArea0_BG # Amid
    SquashedAreaFraction2=SquashedArea2_BG/SquashedArea0_BG # Ahigh   
            
    # Calculate medians 
    SquashedDiameterLevel1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedDiameterLevel1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedArea1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedArea1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction0_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction0,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction1_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction1,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction2_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction2,ImageStagePosition,-50.5,50.5, 101)

    
    
    print(count)
    
    
        # Load analog data 
    if DiffractionCalcs == 1 : 
        
        
        #PixelRatio = 1
        #DiffractionFileName = 'Ros300_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
        ASTAreaFraction0= np.concatenate((ASTAreaFraction0,[np.nan],ASTAreaFraction0), axis=0)
        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
        ASTAreaFraction1= np.concatenate((ASTAreaFraction1,[np.nan],ASTAreaFraction1), axis=0)
        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
        ASTAreaFraction2= np.concatenate((ASTAreaFraction2,[np.nan],ASTAreaFraction2), axis=0)
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (DiffractionPixelRatio **2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()
    
    
    
    fig=plt.figure(figsize=(12,14))
    plt.rcParams.update({'font.size': 14})
    plt.subplot(3,2,1)
    plt.plot(ImageStagePosition, SquashedDiameterLevel1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, AST_MeanXY_BG1, color = 'r')
    plt.plot( ImageStagePositionMid,SquashedDiameterLevel1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('Diameter, μm')
    plt.xlim([-50,50])
    plt.subplot(3,2,2)
    plt.plot(ImageStagePosition, SquashedArea1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTArea1_BG, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedArea1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 100}$, pixels')
    plt.xlim([-50,50])
    plt.subplot(3,2,3)
    plt.plot(ImageStagePosition, SquashedAreaFraction0,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction0, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction0_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{25 - 50}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.subplot(3,2,4)
    plt.plot(ImageStagePosition, SquashedAreaFraction1,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction1, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction1_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 75}$ / A$_{25 - 100}$')
    plt.xlabel('Z, mm')
    plt.xlim([-50,50])
    plt.subplot(3,2,5)
    plt.plot(ImageStagePosition, SquashedAreaFraction2,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction2, color = 'r', label = 'Angular spectrum theory')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction2_Median,'o', markersize =5, color = 'k', label = 'CIP15')
    plt.legend(bbox_to_anchor=(1.3, 0.5), loc='center left', borderaxespad=0.)
    plt.ylabel('A$_{75 - 100}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.xlabel('Z, mm')
    plt.savefig(ImagePath+'GreyscaleVsZ.png',dpi=200)
    plt.close(fig)
#    i =6601
#    Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], 15, 64)
#    print(ParticleBufferTimeS[i])
#    print(ImageID[i])


    #return ParticleTimeSeconds, Index, IndexHigh, IndexLow, Direction
#__________________________________________________________________________________

def SquashImageY_4_GS(BinaryImage): 
    BinaryImage = np.array(BinaryImage, dtype=float)   
    BinaryImage[BinaryImage == 255] = 12.5 
    BinaryImage[BinaryImage == 1] = 37.5 
    BinaryImage[BinaryImage == 2] = 62.5
    BinaryImage[BinaryImage == 3] = 87.5
    
    TmpSum = np.cumsum(BinaryImage, axis =1)
    Tmp = TmpSum[:,3::4]
    SquashImage = Tmp[:,1::1] - Tmp[:,0:-1:1]
    SquashImage /= 4 
    
#    SquashImage = (BinaryImage[:,:-1:2] + BinaryImage[:,1::2])/ 2  
#    SquashImage[SquashImage < 0.5] = int(0)
#    SquashImage[SquashImage >= 0.5] = int(1)
#    SquashImage= np.array(SquashImage, dtype=int)
    
    SquashImage[SquashImage < 25] = 0 
    #SquashImage  = (np.where(np.logical_and(SquashImage >= 25,SquashImage < 50), 1,SquashImage))
    SquashImage[(SquashImage >= 25) & (SquashImage < 50)] = 1
    #SquashImage  = (np.where(np.logical_and(SquashImage >= 50, SquashImage < 75), 2,SquashImage))
    SquashImage[(SquashImage >= 50) & (SquashImage < 75)] = 2
    SquashImage[SquashImage >= 75] = 3 
    Pad= np.zeros([64,1])
    SquashImage = np.concatenate((Pad,SquashImage,Pad),axis = 1)
    
    
    
    return SquashImage#, StretchImage, TmpSum, Tmp
#__________________________________________________________________________________

#Analyse Cip15 experiment using ROS118 and high speed stage
   
def CIP15_ros118():

    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/261119/ros118/20191126134140/Output/'
    FileName = 'Imagefile_1CIP Grayscale_20191126134140.h5'
    ImageFileName = 'Export_Imagefile_1CIP Grayscale_20191126134140.h5'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/261119/ros118/20191126134140/Output/FilteredImages/'
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/261119/ros118/'
    StageName = 'Stagelog.txt' 
    DiffractionCalcs  = 1
    DiffractionPath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/AnalogImages/'
    DiffractionFileName = 'Ros118_.h5'
    DiffractionPixelRatio = 0.8/15
    
    Lambda= 0.658 # laser wavelength
    Z0 = 101
    D0 = 118
    PixelSize = 15
    ArrayWidth = 64
    
    ParticleBufferTimeS, ParticleTimeSeconds, ImageID, ParticleStatsWv, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode, AspectRatioLevel1  = LoadCIP_h5(ImagePath, FileName)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG, Zd_fromRatios = GrayscaleRatios(0, 0, ParticleStatsWv,0,'255075')

    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    #StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan

    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    
    AspectRatioBBOX = (YMaxDiode + 1 - YMinDiode) / (XMaxDiode + 1 - XMinDiode)
    
    #plt.plot(ImageStagePosition, AspectRatioBBOX)
    
    AnalogFlag = np.zeros(len(ImageID))
    # Search for analog images 
        
    SquashedDiameterLevel0_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel1_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel2_BG=np.zeros(len(ImageID))*np.nan
    SquashedArea0_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea1_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea2_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea0_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea1_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea2_Filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity0_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity1_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity2_filled = np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel0=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel1=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel2=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGy1=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGx1=np.zeros(len(ImageID))*np.nan
    
    count = 0    
    AnalogFlag[(np.absolute(ImageStagePosition)>7 ) & (AspectRatioBBOX>2) & (XMinDiode > 15) & (XMaxDiode < 50) & ((DiameterLevel0 > 100)|(DiameterLevel1 >= 40)|(DiameterLevel2 >= 40))] = 1
    AnalogFlag[(np.absolute(ImageStagePosition)<7 ) & (AreaFraction0 <0.9)& ((DiameterLevel0 < 1000)|(DiameterLevel1 < 1000)|(DiameterLevel2 < 1000))& ((DiameterLevel0 > 75)|(DiameterLevel1 >75)|(DiameterLevel2 >75))] = 1
            
    for i in range(len(ImageID)):
        if AnalogFlag[i] == 1:
            count+=1
            if (count % 1 == 0) :
#                PlotImage_GS(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],ImageZd[i], ImageStagePosition[i], ImageSavePath)

                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                SquashedImage = SquashImageY_4_GS(Image)
                #Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                PlotOAPImageGS(SquashedImage,PixelSize, ImageSavePath,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]))
                
                I_binned_75 = (np.where(SquashedImage>=3, 1, 0))
                I_binned_50 = (np.where(SquashedImage>=2, 1, 0))
                I_binned_25 = (np.where(SquashedImage>=1, 1, 0))

                SquashedDiameterLevel0_BG[i], SquashedArea0_BG[i], tmp, tmp, SquashedCircularity0_filled[i],SquashedBoxDiameterLevel0[i],tmp, SquashedArea0_Filled[i] = ImageParticleMoreStats(I_binned_25, PixelSize)
                SquashedDiameterLevel1_BG[i], SquashedArea1_BG[i], SquashedDiameterBGx1[i],SquashedDiameterBGy1[i], SquashedCircularity1_filled[i],SquashedBoxDiameterLevel1[i],tmp, SquashedArea1_Filled[i] = ImageParticleMoreStats(I_binned_50, PixelSize)
                SquashedDiameterLevel2_BG[i], SquashedArea2_BG[i], tmp, tmp,SquashedCircularity2_filled[i],SquashedBoxDiameterLevel2[i],tmp, SquashedArea2_Filled[i] = ImageParticleMoreStats(I_binned_75, PixelSize)
                
        
    SquashedAreaFraction0=(SquashedArea0_BG-SquashedArea1_BG)/SquashedArea0_BG # Alow notation from OShea AMT 2019   
    SquashedAreaFraction1=(SquashedArea1_BG-SquashedArea2_BG)/SquashedArea0_BG # Amid
    SquashedAreaFraction2=SquashedArea2_BG/SquashedArea0_BG # Ahigh   
            
    # Calculate medians 
    SquashedDiameterLevel1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedDiameterLevel1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedArea1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedArea1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction0_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction0,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction1_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction1,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction2_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction2,ImageStagePosition,-50.5,50.5, 101)
   
#    AspectRatioBBOX[AnalogFlag == 0] = np.nan
#    plt.plot(ImageStagePosition, AspectRatioBBOX,'o')
    print(count)
#    
#    
#        # Load analog data 
    if DiffractionCalcs == 1 : 
        
        
        #PixelRatio = 1
        #DiffractionFileName = 'Ros300_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
        ASTAreaFraction0= np.concatenate((ASTAreaFraction0,[np.nan],ASTAreaFraction0), axis=0)
        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
        ASTAreaFraction1= np.concatenate((ASTAreaFraction1,[np.nan],ASTAreaFraction1), axis=0)
        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
        ASTAreaFraction2= np.concatenate((ASTAreaFraction2,[np.nan],ASTAreaFraction2), axis=0)
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (DiffractionPixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()
#    
#    
#    
#    #DiameterLevel1[AnalogFlag==0] =np.nan
#    #AspectRatioLevel1[AnalogFlag==0] =np.nan
    fig=plt.figure(figsize=(12,14))
    plt.rcParams.update({'font.size': 14})
    plt.subplot(3,2,1)
    plt.plot(ImageStagePosition, SquashedDiameterLevel1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, AST_MeanXY_BG1, color = 'r')
    plt.plot( ImageStagePositionMid,SquashedDiameterLevel1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('Diameter, μm')
    plt.xlim([-50,50])
    plt.subplot(3,2,2)
    plt.plot(ImageStagePosition, SquashedArea1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTArea1_BG, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedArea1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 100}$, pixels')
    plt.xlim([-50,50])
    plt.subplot(3,2,3)
    plt.plot(ImageStagePosition, SquashedAreaFraction0,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction0, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction0_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{25 - 50}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.subplot(3,2,4)
    plt.plot(ImageStagePosition, SquashedAreaFraction1,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction1, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction1_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 75}$ / A$_{25 - 100}$')
    plt.xlabel('Z, mm')
    plt.xlim([-50,50])
    plt.subplot(3,2,5)
    plt.plot(ImageStagePosition, SquashedAreaFraction2,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction2, color = 'r', label = 'Angular spectrum theory')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction2_Median,'o', markersize =5, color = 'k', label = 'CIP15')
    plt.legend(bbox_to_anchor=(1.3, 0.5), loc='center left', borderaxespad=0.)
    plt.ylabel('A$_{75 - 100}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.xlabel('Z, mm')
    plt.savefig(ImagePath+'GreyscaleVsZ.png',dpi=200)
    plt.close(fig)
#    
#    return XMinDiode, XMaxDiode
#    return AnalogFlag, AspectRatioBBOX, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode
    
    
    
    
#__________________________________________________________________________________
  
#Analyse Cip15 experiment using ROS250 and high speed stage

def CIP15_ros250():

    ImagePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros250_run3/20191125163508/Output/'
    FileName = 'Imagefile_1CIP Grayscale_20191125163508.h5'
    ImageFileName = 'Export_Imagefile_1CIP Grayscale_20191125163508.h5'
    ImageSavePath = 'C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros250_run3/20191125163508/Output/FilteredImages/'
    StagePath ='C:/Users/Admin TEMP/Documents/HighSpeedStage/CIP15/251119/ros250_run3/'
    StageName = 'Stagelog.txt' 
    DiffractionCalcs  = 1
    DiffractionPath =ImagePath
    DiffractionFileName = 'Ros250_CIP15.h5'
    DiffractionPixelRatio = 1/15
    #DiffractionFileName = 'Ros250_Microscope.h5'
    #DiffractionPixelRatio = 1.25/15
       
    Lambda= 0.658 # laser wavelength
    Z0 = 101
    D0 = 250
    PixelSize = 15
    ArrayWidth = 64
    
    ParticleBufferTimeS, ParticleTimeSeconds, ImageID, ParticleStatsWv, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode, AspectRatioLevel1  = LoadCIP_h5(ImagePath, FileName)
    AreaFraction0,AreaFraction1,AreaFraction2,AreaRatio_2_1,AreaRatio_1_0,DiameterLevel0,DiameterLevel1,DiameterLevel2, FilledArea, VoidArea,DiameterBG, Zd_fromRatios = GrayscaleRatios(0, 0, ParticleStatsWv,0,'255075')

    StageTime, StagePosition, StageSecondsSinceMidnight = LoadStageLog(StagePath,StageName,1)        
    #StageSecondsSinceMidnight -= 3600 # The stage is an hour ahead of the 2DS. Stage BST
    Idx = np.searchsorted(ParticleBufferTimeS,StageSecondsSinceMidnight)
    ImageStagePosition  = np.zeros(len(ParticleBufferTimeS))*np.nan

    for i in range(len(Idx)-1) :     
        ImageStagePosition[Idx[i]:Idx[i+1]] = StagePosition[i]
    ImageStagePosition -= Z0
    Z_um = ImageStagePosition * 1000
    ImageZd = (4 * Lambda * Z_um) / ((D0)**2)

    
    AspectRatioBBOX = (YMaxDiode + 1 - YMinDiode) / (XMaxDiode + 1 - XMinDiode)
    
    #plt.plot(ImageStagePosition, AspectRatioBBOX)
    
    AnalogFlag = np.zeros(len(ImageID))
    # Search for analog images 
        
    SquashedDiameterLevel0_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel1_BG=np.zeros(len(ImageID))*np.nan
    SquashedDiameterLevel2_BG=np.zeros(len(ImageID))*np.nan
    SquashedArea0_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea1_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea2_BG = np.zeros(len(ImageID))*np.nan
    SquashedArea0_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea1_Filled = np.zeros(len(ImageID))*np.nan
    SquashedArea2_Filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity0_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity1_filled = np.zeros(len(ImageID))*np.nan
    SquashedCircularity2_filled = np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel0=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel1=np.zeros(len(ImageID))*np.nan
    SquashedBoxDiameterLevel2=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGy1=np.zeros(len(ImageID))*np.nan
    SquashedDiameterBGx1=np.zeros(len(ImageID))*np.nan
    
    count = 0

    #Select particles    
    AnalogFlag[(XMinDiode > 10) & (XMaxDiode < 50) & (DiameterLevel0 > 200) & (DiameterLevel0 < 750) & (AreaFraction0 <0.9)] = 1   
    #AnalogFlag[(XMinDiode > 15) & (XMaxDiode < 50) & (AspectRatioBBOX>3) & (AspectRatioBBOX<5) & ((DiameterLevel0 > 100)|(DiameterLevel1 >= 40)|(DiameterLevel2 >= 40))] = 1 #


#    ParticleTimeSeconds[AnalogFlag != 1] = np.nan
    #use differnce in particle element number to determine stage direction
#    Index=np.where(ParticleTimeSeconds==ParticleTimeSeconds)
#    Index=np.array(Index[0]).astype(float)
#    IndexHigh = Index[1:] - Index[0:-1]
#    IndexHigh = np.append(IndexHigh,[np.nan] )
#    IndexLow = Index[1:] - Index[0:-1]
#    IndexLow = np.append([np.nan],IndexLow )
#    Direction = np.zeros(len(ParticleTimeSeconds))*np.nan
#    Index = Index[IndexHigh < IndexLow] # > for forward, < for backward
#    Index = np.array(Index).astype(int)
#    Direction[Index] = 1
            
    for i in range(len(ImageID)):
        if AnalogFlag[i] == 1:
            count+=1
            if (count % 1 == 0) :
                #PlotImage_GS(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i],ImageZd[i], ImageStagePosition[i], ImageSavePath)

                Image, x, y = LoadImage(ImagePath,ImageFileName,ParticleBufferTimeS[i],ImageID[i], PixelSize, ArrayWidth)
                SquashedImage = SquashImageY_4_GS(Image)
                #Squash_MeanXY_BG[i], Squash_Area_BG[i], Squash_DiameterBGx[i],Squash_DiameterBGy[i], Squash_Circularity_Filled[i], Squash_MeanXY[i],tmp,Squash_Area_Filled[i] = ImageParticleMoreStats(SquashedImage, 10)
                #PlotOAPImageGS(SquashedImage,PixelSize, ImageSavePath,ParticleBufferTimeS[i],ImageID[i],np.absolute(ImageZd[i]), np.absolute(ImageStagePosition[i]))
                
                I_binned_75 = (np.where(SquashedImage>=3, 1, 0))
                I_binned_50 = (np.where(SquashedImage>=2, 1, 0))
                I_binned_25 = (np.where(SquashedImage>=1, 1, 0))

                SquashedDiameterLevel0_BG[i], SquashedArea0_BG[i], tmp, tmp, SquashedCircularity0_filled[i],SquashedBoxDiameterLevel0[i],tmp, SquashedArea0_Filled[i] = ImageParticleMoreStats(I_binned_25, PixelSize)
                SquashedDiameterLevel1_BG[i], SquashedArea1_BG[i], SquashedDiameterBGx1[i],SquashedDiameterBGy1[i], SquashedCircularity1_filled[i],SquashedBoxDiameterLevel1[i],tmp, SquashedArea1_Filled[i] = ImageParticleMoreStats(I_binned_50, PixelSize)
                SquashedDiameterLevel2_BG[i], SquashedArea2_BG[i], tmp, tmp,SquashedCircularity2_filled[i],SquashedBoxDiameterLevel2[i],tmp, SquashedArea2_Filled[i] = ImageParticleMoreStats(I_binned_75, PixelSize)
                
        
    SquashedAreaFraction0=(SquashedArea0_BG-SquashedArea1_BG)/SquashedArea0_BG # Alow notation from OShea AMT 2019   
    SquashedAreaFraction1=(SquashedArea1_BG-SquashedArea2_BG)/SquashedArea0_BG # Amid
    SquashedAreaFraction2=SquashedArea2_BG/SquashedArea0_BG # Ahigh   

# Calculate medians 
    SquashedDiameterLevel1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedDiameterLevel1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedArea1_BG_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedArea1_BG,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction0_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction0,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction1_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction1,ImageStagePosition,-50.5,50.5, 101)
    SquashedAreaFraction2_Median, BinsLower, BinsUpper, ImageStagePositionMid = BinLinear_median(SquashedAreaFraction2,ImageStagePosition,-50.5,50.5, 101)


#            
#    
#    AspectRatioBBOX[AnalogFlag == 0] = np.nan
#    plt.plot(ImageStagePosition, AspectRatioBBOX,'o')
#    print(count)
##    
##    
##        # Load analog data 
    if DiffractionCalcs == 1 : 
        
        
        #PixelRatio = 1
        #DiffractionFileName = 'Ros300_2DSimage.h5'
        Data_h5 = h5py.File(DiffractionPath+DiffractionFileName, 'r')              
        AST_Zarray_mm = np.array(Data_h5['Zarray']) / 1000        
        AST_Zarray_mm = np.concatenate((-AST_Zarray_mm,[np.nan],AST_Zarray_mm), axis=0)
        
        
        AST_MeanXY_BG1=np.array(Data_h5['DiameterLevel1'])
        AST_MeanXY_BG1 = np.concatenate((AST_MeanXY_BG1,[np.nan],AST_MeanXY_BG1), axis=0)
        ASTAreaFraction0=np.array(Data_h5['AreaFraction0'])
        ASTAreaFraction0= np.concatenate((ASTAreaFraction0,[np.nan],ASTAreaFraction0), axis=0)
        ASTAreaFraction1=np.array(Data_h5['AreaFraction1'])
        ASTAreaFraction1= np.concatenate((ASTAreaFraction1,[np.nan],ASTAreaFraction1), axis=0)
        ASTAreaFraction2=np.array(Data_h5['AreaFraction2'])
        ASTAreaFraction2= np.concatenate((ASTAreaFraction2,[np.nan],ASTAreaFraction2), axis=0)
#        ASTCircularity0_filled=np.array(Data_h5['Circularity0_filled'])
        ASTCircularity1_filled=np.array(Data_h5['Circularity1_filled'])
        ASTCircularity1_filled = np.concatenate((ASTCircularity1_filled,[np.nan],ASTCircularity1_filled), axis=0)
#        ASTCircularity2_filled=np.array(Data_h5['Circularity2_filled'])
        ASTArea1_BG=np.array(Data_h5['Area1_BG']) * (DiffractionPixelRatio**2)
        ASTArea1_BG = np.concatenate((ASTArea1_BG,[np.nan],ASTArea1_BG), axis=0)
        
        ASTBoxDiameterLevel1=np.array(Data_h5['BoxDiameterLevel1'])
        Data_h5.close()
##    
##    
##    
##    #DiameterLevel1[AnalogFlag==0] =np.nan
##    #AspectRatioLevel1[AnalogFlag==0] =np.nan
    fig=plt.figure(figsize=(12,14))
    plt.rcParams.update({'font.size': 14})
    plt.subplot(3,2,1)
    plt.plot(ImageStagePosition, SquashedDiameterLevel1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, AST_MeanXY_BG1, color = 'r')
    plt.plot( ImageStagePositionMid,SquashedDiameterLevel1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('Diameter, μm')
    plt.xlim([-50,50])
    plt.subplot(3,2,2)
    plt.plot(ImageStagePosition, SquashedArea1_BG,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTArea1_BG, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedArea1_BG_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 100}$, pixels')
    plt.xlim([-50,50])
    plt.subplot(3,2,3)
    plt.plot(ImageStagePosition, SquashedAreaFraction0,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction0, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction0_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{25 - 50}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.subplot(3,2,4)
    plt.plot(ImageStagePosition, SquashedAreaFraction1,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction1, color = 'r')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction1_Median,'o', markersize =5, color = 'k')
    plt.ylabel('A$_{50 - 75}$ / A$_{25 - 100}$')
    plt.xlabel('Z, mm')
    plt.xlim([-50,50])
    plt.subplot(3,2,5)
    plt.plot(ImageStagePosition, SquashedAreaFraction2,'o', markersize =1, color = 'silver')
    plt.plot(AST_Zarray_mm, ASTAreaFraction2, color = 'r', label = 'Angular spectrum theory')
    plt.plot(ImageStagePositionMid,SquashedAreaFraction2_Median,'o', markersize =5, color = 'k', label = 'CIP15')
    plt.legend(bbox_to_anchor=(1.3, 0.5), loc='center left', borderaxespad=0.)
    plt.ylabel('A$_{75 - 100}$ / A$_{25 - 100}$')
    plt.xlim([-50,50])
    plt.xlabel('Z, mm')
    plt.savefig(ImagePath+'GreyscaleVsZ.png',dpi=200)
    plt.close(fig)
##    
##    return XMinDiode, XMaxDiode
##    return AnalogFlag, AspectRatioBBOX, XMaxDiode, XMinDiode,YMaxDiode, YMinDiode    
    