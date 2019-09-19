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
from GrayscaleFunctions import GetIAT_TimeInS_vector, BinHaloHolo_v2
from MyFunctions import datetime_2_s_sincemidnight
from PICASSO_functions import Avg_CDP_dndDp


