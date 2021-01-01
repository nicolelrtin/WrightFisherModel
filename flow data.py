%pylab inline

pip install FlowCytometryTools

import fcsparser, os, FlowCytometryTools
import pandas as pd
import scipy as sp 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import FlowCytometryTools
from FlowCytometryTools import FCMeasurement
sns.set(style="whitegrid")

C9 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G1_C9.fcs')
C10 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G1_C10.fcs')
C11 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G1_C11.fcs')
C5 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G2_C5.fcs')
B3 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G2_B3.fcs')
B7 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G2_B7.fcs')
B10 = FCMeasurement(ID='Test Sample', datafile = '1018_Well Mixed_G2_B10.fcs')
samples = [C9, C10, C11, C5, B3, B7, B10]
print(C9.data)