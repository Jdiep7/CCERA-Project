# Template for Sun Scan lab

import numpy as np
import matplotlib.pyplot as plt
import json 
from math import factorial, pi 
import time 
from datetime import datetime 

# Airy function .   We will use this in Lab03_02.py 

def airy(mean_time, base_temp, peak_temp, width) : 
    times = np.linspace(mean_time-1000.,mean_time+1000.,200)
    airy_function = [] 
    for tt in times :
        r = (tt-mean_time)/width 
        t, I = 0.5*pi*r, 1.
        for k in range(1,20) :
            dI = (-1)**k * t**(2*k) / (factorial(k)*factorial(k+1))
            I = I + dI 
        I *= I
        airy_function.append(I) 
    return times, base_temp + peak_temp*np.array(airy_function)

# begin execution

# this line specifies the data set to be used 
base_name = "./Lab03_data/2024-07-03-1340"

# get the metadata, which provides information about the data file we are reading 
with open(base_name + ".json") as json_file : metadata = json.load(json_file)
print("file={0:s} \nmetadata={1:s}".format(base_name,str(metadata)))

# read in the data time series
data_file = base_name + "_2.sum" 
power = np.fromfile(data_file, dtype=np.float32)
nVals = len(power)
print("Number of data values read={0:d}".format(nVals))

