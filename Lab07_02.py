import matplotlib.pyplot as plt
import numpy as np
import glob

def getMetaData(file) :
    import json
    with open(file) as json_file:
        dict = json.load(json_file)
        #print("From file {0:s} \nread dictionary={1:s}".format(file,str(dict)))
    return dict 

def getData(file,fft_size) :
    vals = np.fromfile(file, dtype=np.float32)
    cols = fft_size
    rows = int(len(vals)/fft_size) 
    return vals, rows, cols

def getFreqs(metadata) :
    dF = 1.e-6*metadata['srate']
    fMin = 1.e-6*metadata['freq'] - 0.5*dF
    fMax = 1.e-6*metadata['freq'] + 0.5*dF
    #print("getFreqs: fMin={0:e} fMax={1:e}".format(fMin,fMax))
    freqs = np.linspace(fMin,fMax,metadata['fft_size'])
    return freqs

def getVelocities(f) :
    f0, c = 1420.41, 3.0e5
    v = c*(f/f0 - 1.)
    return v

# fit spectrum to Chebyshev polynomial 
# restrict range of fit to |vDoppler| > vSignal 
def fitBackground(vDoppler,power,n,vSignal) :
    weights = np.ones_like(vDoppler)
    for i in range(len(vDoppler)) :
        if abs(vDoppler[i]) < vSignal : weights[i] = 1.e-6 
    series = np.polynomial.chebyshev.Chebyshev.fit(vDoppler, power, n, w=weights)
    background = series(vDoppler) 
    return background

def anaSpectrum(base_name) :
    #print("In anaSpectrum() base_name={0:s}".format(base_name))
    metadata = getMetaData(base_name + ".json")
    fft_size = metadata['fft_size']

    # we will use channel 1 
    chan = 1
    data_file = base_name + "_{0:d}.avg".format(chan)
    power, rows, cols = getData(data_file,fft_size)

    freqs = getFreqs(metadata)
    vDoppler = getVelocities(freqs)
    power *= 1.10e5

    vMin, vMax = -300., 300.
    i1 = np.searchsorted(vDoppler,vMin)
    i2 = np.searchsorted(vDoppler,vMax)
    #print("i1={0:d} i2={1:d}".format(i1,i2))
    freqs = freqs[i1:i2]
    vDoppler = vDoppler[i1:i2]
    power = power[i1:i2]

    background = fitBackground(vDoppler,power,5,200.)
  
    return vDoppler, power-background 


# begin execution here

# Get a list of files.  Make sure that they are ordered by time
files = glob.glob("./Lab07_data/*.json")
files.sort() 
print("files={0:s}".format(str(files)))

# analyse the first file to establish the parameters of the plot
base_name = "." + files[0].strip(".json")
vDoppler, power = anaSpectrum(base_name)
nRows, nCols = len(files), len(vDoppler)
print("nRows={0:d} nCols={1:d}".format(nRows,nCols))
mapData = np.zeros((nRows,nCols))
calculatedvLSR=[]
glonList = []
vPrime = []

for i, file in enumerate(files) :
    metadata = getMetaData(file)
    base_name =  "." + file.strip(".json")
    vDoppler, power = anaSpectrum(base_name)
    gLon = float(metadata['gLon'])
    row = int(gLon/3 + 0.5)
    mapData[row] = np.maximum(power,0.)
    vLSR = metadata['vlsr']
    calculatedvLSR.append(vLSR)
    glonList.append(gLon)

    for i, p in enumerate(power):
        if p > 5.:
            vPrime.append(vDoppler[i])
            break

fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111)
ax.set_title("Galactic Scan")
ax.set_xlabel("Approach velocity (km/s)")
ax.set_ylabel("Galactic Longitude (deg)")

ax.patch.set_facecolor('white')

im = ax.imshow(mapData,extent=[-300.,300.,90.,0.],aspect='auto')
im.set_cmap('jet')
plt.colorbar(im, use_gridspec=True)

ax.plot(calculatedvLSR, glonList, color='white', linestyle='-', linewidth=2, markersize=5, label='vLSR')
ax.plot(vPrime, glonList, color='pink', linestyle='-', linewidth=2, markersize=5, label='vPrime')
plt.show()