import matplotlib.pyplot as plt
import numpy as np
import time
import socket


def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

# fit spectrum to Chebyshev polynomial 
# restrict range of fit to |vDoppler| > vSignal 
def fitBackground(vDoppler,power,n,vSignal) :
    weights = np.ones_like(vDoppler)
    for i in range(len(vDoppler)) :
        if abs(vDoppler[i]) < vSignal : weights[i] = 1.e-6 
    series = np.polynomial.chebyshev.Chebyshev.fit(vDoppler, power, n, w=weights)
    background = series(vDoppler) 
    return background

# Begin execution here

# read in the metadata and the data 
base_name = './Lab05_data/2024-07-09-2244'
metadata = getMetaData(base_name + ".json")
fft_size = metadata['fft_size']

# we will use channel 1 
chan = 1
data_file = base_name + "_{0:d}.avg".format(chan)
power, rows, cols = getData(data_file,fft_size)

fCenter = 1.0e-6*metadata['freq']
f_sample = 1.0e-6*metadata["srate"]

fMin = fCenter - (f_sample/2)
fMax = fCenter + (f_sample/2)

freqs = np.linspace(fMin,fMax,metadata['fft_size'])

power *= 1.3e5

vDoppler = ((freqs - 1420.41)/1420.41)*(3e5)

plt.plot(vDoppler, power, 'b.')
plt.title("Antenna Temperature based on Doppler Velocity")
plt.xlabel("Doppler Velocity")
plt.ylabel("Antenna Temperature")
plt.show()

print("rows={0:d} cols={1:d}".format(rows,cols))