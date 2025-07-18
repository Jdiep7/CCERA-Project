import matplotlib.pyplot as plt
import numpy as np
import glob
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
import os

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

#def vlsr(t,loc,psrc,verbose=False):
def vlsr(metadata,loc,verbose=False):
    """Compute the line of sight radial velocity

    psrc: SkyCoord object or source
    loc: EarthLocation object of observer
    t: Time object
    """
    tra, tdec = metadata['RA'], metadata['dec']
    psrc = SkyCoord(ra = tra, dec = tdec ,frame = "icrs", unit = (u.deg,u.deg)) 
    t = Time(metadata['t_start'],scale="utc",format="unix")

    # Direction of motion of the Sun. Info from
    psun = SkyCoord(ra = "18:03:50.29", dec = "+30:00:16.8",frame = "icrs",unit = (u.hourangle,u.deg))
    vsun = -20.0*u.km/u.s

    # Radial velocity correction to solar system barycenter
    vsrc = psrc.radial_velocity_correction(obstime = t, location = loc)

    # Projection of solar velocity towards the source
    vsun_proj = psrc.cartesian.dot(psun.cartesian)*vsun

    if verbose:
        print("Barycentric radial velocity: {0:+8.3f}".format(vsrc.to(u.km/u.s)))
        print("Projected solar velocity:    {0:+8.3f}".format(vsun_proj.to(u.km/u.s)))
    
    return vsun_proj-vsrc

# Get a list of files.  Make sure that they are ordered by time
files = glob.glob("./AL045/*.json")
files.sort() 
print("files={0:s}".format(str(files)))

# analyse the first file to establish the parameters of the plot
base_name = os.path.splitext(files[0])[0]
vDoppler, power = anaSpectrum(base_name)

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

for row, file in enumerate(files) :
    base_name = "./" + file.strip(".json")
    vDoppler, power = anaSpectrum(base_name)

    li.set_xdata(vDoppler)  
    li.set_ydata(power)

    timeString = os.path.basename(base_name)
    timeText.set_text(timeString) 
    plt.pause (0.5)
    #time.sleep(0.5)

#plt.plot(vals,times,'w.')
plt.show()

