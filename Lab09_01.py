# analyze pulsar data

import numpy as np
import matplotlib.pyplot as plt 
import json 

# get the JSON file 
base_name = './data/2024-07-26-1804'
with open(base_name + ".json") as json_file : metadata = json.load(json_file)

# read or calculate various run parameters 
f_sample = metadata['srate']
fft_size = metadata['fft_size']
n_decimate = metadata['decimation_factor']
c_rate = f_sample/fft_size/n_decimate
t_fft = 1./c_rate 

# for now, we will look at just a single file (there are two, 
# one for each polarization) 
file = base_name + "_1.sum"
power_time_series = 1000*np.fromfile(file, dtype=np.float32)
nSamples = len(power_time_series)
print("nSamples={0:d}".format(nSamples))
times = np.linspace(0.,nSamples*t_fft,nSamples) 