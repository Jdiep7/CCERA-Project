import matplotlib.pyplot as plt
import numpy as np
import json

# Function to read metadata
def readMetadata(base_name):
    with open(base_name + ".json") as json_file:
        metadata = json.load(json_file)
        print("From file {0:s} \nread dictionary={1:s}".format(base_name, str(metadata)))
    return metadata

# Function to compute frequency array
def getFreqs(metadata):
    fCenter = 1.0e-6 * metadata['freq']   # MHz
    f_sample = 1.0e-6 * metadata["srate"] # MHz
    fMin = fCenter - (f_sample / 2)
    fMax = fCenter + (f_sample / 2)
    freqs = np.linspace(fMin, fMax, metadata['fft_size'])
    return freqs

# Function to compute doppler velocities
def getVdoppler(freqs):
    vDoppler = ((freqs - 1420.41) / 1420.41) * (3e5)  # km/s
    return vDoppler


def getPower(base_name):
    metadata = readMetadata(base_name)
    chan = 2
    data_file = base_name + "_{0:d}.avg".format(chan)
    vals = np.fromfile(data_file, dtype=np.float32)
    power = vals * 1.3e5
    return power

def limitRange(v_doppler, power, v_min, v_max):
    i1 = np.searchsorted(v_doppler, v_min)
    i2 = np.searchsorted(v_doppler, v_max)
    return v_doppler[i1:i2], power[i1:i2]

def getBackground(v_doppler, power, n=5, vSignal=200):
    weights = np.ones_like(v_doppler)
    for i in range(len(v_doppler)):
        if abs(v_doppler[i]) < vSignal:
            weights[i] = 1.e-6
    series = np.polynomial.chebyshev.Chebyshev.fit(v_doppler, power, n, w=weights)
    background = series(v_doppler)
    return background

base_name = './Lab05_data/2024-07-10-0514'

metadata = readMetadata(base_name)
freqs = getFreqs(metadata)
v_doppler = getVdoppler(freqs)
power = getPower(base_name)


plt.plot(freqs, power, 'b.')
plt.title('Power vs Frequency')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Antenna Temperature")
plt.show()

plt.plot(v_doppler, power, 'b.')
plt.title('Power vs Doppler Velocity')
plt.xlabel("Doppler Velocity (km/s)")
plt.ylabel("Antenna Temperature")
plt.show()


v_min, v_max = -300, 300
v_doppler, power = limitRange(v_doppler, power, v_min, v_max)
background = getBackground(v_doppler, power)

plt.plot(v_doppler, power, 'b.', label='Power')
plt.plot(v_doppler, background, 'r-', label='Background fit')
plt.title('Power and Background vs Doppler Velocity')
plt.xlabel("Doppler Velocity (km/s)")
plt.ylabel("Antenna Temperature")
plt.legend()
plt.show()

power -= background

plt.plot(v_doppler, power, 'b.')
plt.title('Doppler spectrum for {0:s}'.format(base_name.split("/")[-1]))
plt.xlabel("Doppler Velocity")
plt.ylabel("Antenna Temperature")
plt.show()