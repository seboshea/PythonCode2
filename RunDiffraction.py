
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
from MyFunctions import BinLinearV2
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,AutoMinorLocator
import bisect
from AngularSpectrumTheory import DmeasVsD0_habit, DmeasD0VsD0_habit


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
h5f.create_dataset('D0arrayEdge',data=D0arrayEdge)
h5f.create_dataset('D0arrayEdge',data=D0arrayEdge)
h5f.close()