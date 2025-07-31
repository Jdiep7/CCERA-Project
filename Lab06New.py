import matplotlib.pyplot as plt
import numpy as np
from math import sqrt, sin
import glob
import os
from datetime import datetime

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

    vDoppler = ((freqs - 1420.41)/1420.41)*(3.0e5)

    v1, v2 = -300., 300.
    i1 = np.searchsorted(vDoppler,v1)
    i2 = np.searchsorted(vDoppler,v2)

    freqs = freqs[i1:i2]
    vDoppler = vDoppler[i1:i2]
    power = power[i1:i2]

    background = fitBackground(vDoppler, power, 5, 200)

    power -= background
    
    return vDoppler, power


def getBaseNames(t_start=None, t_stop=None):
    files = glob.glob("./AL045/*.json")
    files.sort()
    time_section = []    

    if t_start is None and t_stop is None:
        time_section = files
    elif t_start is None:
        end_file =  files.index("./AL045\\" + t_stop + '.json') + 1
        time_section = files[:end_file]
    elif t_stop is None:
        start_file = files.index("./AL045\\" + t_start + '.json')
        time_section = files[start_file:]
    else:  
        try:
            start_file = files.index("./AL045\\" + t_start + '.json')
            end_file =  files.index("./AL045\\" + t_stop + '.json') + 1
        except Exception:
            print("Time Stamp Not Found")
            return
        
        time_section = files[start_file: end_file]
        
    base_names = []
    
    for i in time_section:
        base_names.append(i.removesuffix(".json"))
        
    return base_names
       
def listMetaData(base_names):
    metadata = []
    for i in base_names:
        metadata.append(getMetaData(i + ".json"))
    return metadata
    
def getSpectra(base_names):
    v_dopplers = []
    powers = []

    for name in base_names:
        vDoppler, power = anaSpectrum(name)
        power = np.maximum(0., power)
        v_dopplers.append(vDoppler)
        powers.append(power)
    
    return v_dopplers, powers

def plotHeatMap(v_dopplers, powers, metadata):
    vDoppler = v_dopplers[0]  
    nRows = len(powers)
    nCols = len(vDoppler)

    mapData = np.zeros((nRows, nCols))
    times = []
    calculatedvLSR = []

    start_time = metadata[0]['t_start']
    end_time = metadata[-1]['t_start'] + metadata[-1]['run_time']
    total_duration_hours = (end_time - start_time) / 3600.0
    start_time_str = datetime.fromtimestamp(start_time).strftime('%Y-%m-%d %H:%M:%S')

    for row, (power, meta) in enumerate(zip(powers, metadata)):
        mapData[row] = power
        rel_time = (meta['t_start'] - start_time) / 3600.0
        times.append(rel_time)
        vLSR = meta['vlsr']
        calculatedvLSR.append(vLSR)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        mapData,
        extent=[vDoppler[0], vDoppler[-1], 0, total_duration_hours],
        aspect='auto',
        origin='lower',
        cmap='jet'
    )

    data_series_name = os.path.basename(os.path.dirname("./AL045/"))
    plot_title = f"HI Spectrum Time Series: {data_series_name}\nStart Time: {start_time_str}"
    ax.plot(calculatedvLSR, times, color='black', marker='o', linestyle='-', linewidth=1, markersize=1, label='vLSR')
    ax.set_title(plot_title)
    ax.set_xlabel("Doppler Velocity (km/s)")
    ax.set_ylabel("Time (hours since start)")
    plt.colorbar(im, use_gridspec=True, label='Power')

    plt.show()
    
def plotAnimation(v_dopplers, powers, base_names):        
    # set up for animation
    fig = plt.figure(1)
    #fig.canvas.set_window_title('21cm Spectrum') 
    ax = fig.add_subplot(111)
    vMin, vMax = -200., 200.
    ax.set_xlim([vMin,vMax])
    li, = ax.plot([], [], 'b.')
    ax.set_ylim([-5.,50.])
    ax.set_title("PSD vs Approach Velocity")
    ax.set_xlabel("v (km/s)")
    ax.set_ylabel("PSD (K)")
    timeText = ax.text(vMin+0.5*(vMax-vMin),40.," ",fontsize=14)
    fig.canvas.draw()
    plt.show(block=False)
    
    for vDoppler, power, name in zip(v_dopplers, powers, base_names):
        li.set_xdata(vDoppler)
        li.set_ydata(power)

        # extract timestamp from base name
        timeString = os.path.basename(name)
        timeText.set_text(timeString)

        plt.pause(0.5)
        #time.sleep(0.5)

    
def plotMountainRange(v_dopplers, powers):
    fig, ax = plt.subplots(figsize=(10, 8))

    offset = 10.0
    vDoppler = v_dopplers[0]

    for i, power in enumerate(powers):
        ax.plot(vDoppler, power + i * offset, color='blue', linewidth=1)

    ax.set_xlabel("Doppler Velocity (km/s)")
    ax.set_ylabel("Power + offset")
    ax.set_title("2D Waterfall Plot of HI Spectra")

    plt.tight_layout()
    plt.show()
    


t_start, t_stop = '2024-07-10-1001', '2024-07-10-1306'
base_names = getBaseNames()
metadata = listMetaData(base_names)
v_dopplers, powers = getSpectra(base_names)
plotHeatMap(v_dopplers, powers, metadata)
plotMountainRange(v_dopplers, powers)
plotAnimation(v_dopplers, powers, base_names)





    