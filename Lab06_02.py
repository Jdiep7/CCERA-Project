import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, sin
import glob
import os

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
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


def anaSpectrum(base_name):
    # read in the metadata and the data 
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

    v1, v2 = -300., 300.
    i1 = np.searchsorted(vDoppler,v1)
    i2 = np.searchsorted(vDoppler,v2)

    vDoppler = vDoppler[i1:i2]
    power = power[i1:i2]

    background = fitBackground(vDoppler, power, 5, 200)

    power -= background
    
    return vDoppler, power
    
# Begin execution here


files = glob.glob("./Lab06_data/*.json")
files.sort()

base_name = os.path.splitext(files[0])[0]
vDoppler, power = anaSpectrum(base_name)
nRows, nCols = len(files), len(vDoppler)
mapData = np.zeros((nRows,nCols))

for row, file in enumerate(files):
    base_name = file.removesuffix(".json")
    vDoppler, power = anaSpectrum(base_name)
    mapData[row] = power


fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(
    mapData,
    extent=[vDoppler[0], vDoppler[-1], 0, len(files)],
    aspect='auto',
    origin='lower',
    cmap='viridis'
)

ax.set_title("HI Spectrum Time Series")
ax.set_xlabel("Doppler Velocity (km/s)")
ax.set_ylabel("Time (index)")

plt.colorbar(im, ax=ax, label="Antenna Temperature (K)")
plt.show()



    