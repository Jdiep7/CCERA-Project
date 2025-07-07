#
# This script contains functions designed to simplify the plotting
# of the various theoretical rotation curves.

# compute various curves associated with the Exponential Disk Model
# See Freeman, Ap. J. Vol 160 (1970) 811.
#

import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt, atan
import scipy.special as sps

# define variations basic parameters 
c = 3.0e8
year = 3.154e7 
lightYear = c*year
R_D = 7000.*lightYear
G = 6.67e-11
mSun = 1.98e30
rSun = 26000.*lightYear
fMax = 10.                # determines the range of the plots 

lightYear_km = lightYear/1000


radius = [np.float64(2.4599965847805472e+17), np.float64(2.4568084907759565e+17), np.float64(2.4468091126821142e+17), np.float64(2.4301412581287254e+17), np.float64(2.406804404809573e+17), np.float64(2.376828367896667e+17), np.float64(2.340205722576638e+17), np.float64(2.297454315139858e+17), np.float64(2.2478392723988848e+17), np.float64(2.1919439670035142e+17), np.float64(2.1302923322324234e+17), np.float64(2.06418792627035e+17), np.float64(1.9915048473344138e+17), np.float64(1.9130664625210038e+17), np.float64(1.829742203654664e+17), np.float64(1.741116455350067e+17), np.float64(1.6476959358791392e+17), np.float64(1.5505736906331286e+17), np.float64(1.4493624742706442e+17), np.float64(1.3428680553769192e+17), np.float64(1.2324494688501624e+17), np.float64(1.120133948118845e+17), np.float64(1.0035536621491885e+17), np.float64(8.847661686939755e+16), np.float64(7.64642044143589e+16), np.float64(6.400322706493743e+16), np.float64(5.141660533470054e+16), np.float64(3.8857114110769384e+16), np.float64(2.6029279804145464e+16), np.float64(1.3216731290515486e+16), np.float64(310138103676114.75)]
vRotation = [np.float64(212.50548666471434), np.float64(214.81319650374465), np.float64(217.37682805737325), np.float64(221.04599629772477), np.float64(218.97074364854583), np.float64(211.99989050392503), np.float64(213.4576971258418), np.float64(212.22376296339428), np.float64(211.66873573512098), np.float64(216.99986040086247), np.float64(217.960221200087), np.float64(217.63123745360454), np.float64(221.0470606912392), np.float64(225.25274139513277), np.float64(222.5783543083386), np.float64(217.301138345224), np.float64(213.7481577706122), np.float64(215.03637646400057), np.float64(212.53183106286764), np.float64(216.0206135160327), np.float64(216.59800227499838), np.float64(204.54812681764645), np.float64(200.74247106506874), np.float64(191.16027245461208), np.float64(189.64435547307227), np.float64(125.84289567248427), np.float64(110.95896154097153), np.float64(99.13403060474387), np.float64(65.58854375608243), np.float64(50.56630523822377), np.float64(35.051932607273066)]

min_radius_ly = 10000 * lightYear_km

combined = list(zip(radius, vRotation))

filtered = [(r,v) for r,v in combined if r >= min_radius_ly]

radius, vRotation = zip(*filtered)
radius = list(radius)
vRotation = list(vRotation)

print(radius)

def getKepler(Md) :
    M = Md*mSun*1.0e11
    R = np.linspace(0.1*R_D,fMax*R_D,100)
    vKepler = np.sqrt(np.reciprocal(R))
    vKepler *= sqrt(G*M)
    print("M={0:e}".format(M))
    return R, vKepler 

def getExponentialDisk(Md) :
    M = Md*mSun*1.0e11
    r = np.linspace(0., fMax*R_D, 100)
    arg = 0.5*r/R_D
    I0 = sps.iv(0,arg)
    K0 = sps.kn(0,arg)
    I1 = sps.iv(1,arg)
    K1 = sps.kn(1,arg)

    rsq = np.multiply(r,r)
    factor = 0.5*M*G/(R_D**3)
    vsq = np.multiply(I0,K0) - np.multiply(I1,K1)
    vsq = factor*np.multiply(rsq,vsq)
    v = np.sqrt(vsq)
    return r, v 

def getIsothermalSphere(Md,mIso,aIso) :

    # mIso is the mass of the isothermal sphere in untis of 10^11 solar masses 
    # aIso is the characteristic radius in units of kilo-light years 
    r = np.linspace(0.1*R_D,fMax*R_D,100)
    aIso = 1000.*lightYear*aIso
    mIso = mIso*1.e11*mSun
    print("mIso ={0:e}".format(mIso))
    print("rSun ={0:e}".format(rSun))
    print("aIso ={0:e}".format(aIso))
    rhoIso = mIso/(4.*pi*(rSun**2+aIso**2)*(rSun - aIso*atan(rSun/aIso)))
    print("rhoIso={0:e}".format(rhoIso)) 
    vIso = 4.*pi*rhoIso*G*(rSun**2 + aIso**2)*(1. - np.multiply(aIso/r,np.arctan(r/aIso)))
    dummy, vExpoDisk = getExponentialDisk(Md)
    vsq = np.multiply(vExpoDisk,vExpoDisk) 
    vIso += vsq
    vIso = np.sqrt(vIso)
    mMW = 4*pi*rhoIso*(rSun**2 + aIso**2)*(390000*lightYear)
    print("mMW/mSun={0:.3e}".format(mMW/mSun))  
    return r, vIso 

# begin execution here 

# parameters for the models 
Md = 0.5           # Disk mass in units of mSun*1e11
mIso = 0.3      # M inside r=R in units of mSun*1e11   This may require adjustment. 
aIso = 20.         # a (in k.l-y) parameter in isothermal sphere model
# Add your code here.

rKepler, vKepler = getKepler(Md)
r,v = getExponentialDisk(Md)
rIso, vIso = getIsothermalSphere(Md, mIso, aIso)

plt.plot((np.array(radius)*1000/lightYear)/1000,np.array(vRotation),'ro',label='Data')

plt.plot(0.001*rKepler/lightYear,0.001*vKepler,'b--',label='Kepler')
plt.plot(0.001*r/lightYear, 0.001*v, 'g--', label="Exponential Disk")
plt.plot(0.001*rIso/lightYear, 0.001*vIso, 'r--', label="Isothermal Sphere")
plt.legend()
plt.show()